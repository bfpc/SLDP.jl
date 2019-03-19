import JuMP
import SDDP

columns(m) = [m[:,i] for i in 1:size(m,2)]

function reject_stdout(f)
    _stdout = STDOUT
    stdout_rd, stdout_wr = redirect_stdout()

    ret = try
        f()
    catch e
        println("ERROR in reject_stdout(): ", e)
    finally
        # restore
        redirect_stdout(_stdout)
        close(stdout_wr)
        close(stdout_rd)
    end
    ret
end


# Auxiliary functions to save and load logs
function writelogline!(line, file)
    v = (line.iteration, line.bound, line.lower_statistical_bound, line.upper_statistical_bound,
         line.timecuts, line.simulations, line.timesimulations, line.timetotal)
    print(file, v[1])
    for vi in v[2:end]
        print(file, ",", vi)
    end
    print(file, "\n")
end

function writelog!(m::SDDP.SDDPModel, filename)
    open(filename, "w") do file
        for line in m.log
            writelogline!(line, file)
        end
    end
end

function readlog!(m::SDDP.SDDPModel, filename)
    open(filename, "r") do file
        for li in readlines(file)
            ls = split(li, ",")
            vals = [parse(fieldtype(SDDP.SolutionLog,fn), String(v)) for (v,fn) in zip(ls,fieldnames(SDDP.SolutionLog))]
            push!(m.log, SDDP.SolutionLog(vals...))
        end
    end
end

#
# Functions to calculate the Future cost functions
#

# First, Qfrak, the current estimate, from cuts
function cutvalue(c::ALDCut, t)
    c.v + dot(c.l, t - c.xhat) - c.rho*norm(t - c.xhat, 1)
end

function estimatevf(ac::Vector{ALDCut}, states...)
    vf = Vector{Float64}[]
    if length(states) == 1
        return [maximum(cutvalue(c,ti) for c in ac) for ti in states[1]]
    else
        ts = columns(SDDP.mesh(states...)::Array{Float64,2})
        return [maximum(cutvalue(c,ti) for c in ac) for ti in ts]
    end
end

function Qfrak(m,t,ms, states::Union{Float64, AbstractVector{Float64}}...; ncuts=nothing, nald=ncuts)
    sp = SDDP.getsubproblem(m,t,ms)
    vf = SDDP.valueoracle(sp)

    ac = aldcuts(sp)
    cuts = SDDP.valid_cuts(SDDP.cutoracle(vf))
    if ncuts != nothing
      ac = ac[1:nald]
      cuts = cuts[1:ncuts]
    end

    _, y_cont = SDDP.processvaluefunctiondata(cuts, true, states...)
    if length(ac) == 0
        return y_cont
    else
        y_ald = estimatevf(ac, states...)
        return maximum(hcat(y_cont, y_ald), 2)[:,1]
    end
end

# Then, the value function using exact values from next stage
function setstate!(sp::JuMP.Model, vs::AbstractVector{Float64})
    for (st, v) in zip(SDDP.states(sp), vs)
        lb = JuMP.getlowerbound(st.variable)
        ub = JuMP.getupperbound(st.variable)
        if v < lb
            SDDP.setvalue!(st, lb)
        elseif v > ub
            SDDP.setvalue!(st, ub)
        else
            SDDP.setvalue!(st, v)
        end
    end
end

# For one sample of the noise
function Q(m::SDDP.SDDPModel,t::Int,ms::Int, inoise::Int, states::Union{Float64, AbstractVector{Float64}}...; debug::Bool=false, relaxation::Bool=false)
    # Get subproblem and set solver accordingly
    sp = SDDP.getsubproblem(m,t+1,ms)
    solvers = sp.ext[:solvers]
    if relaxation
        JuMP.setsolver(sp, solvers.LP)
    else
        JuMP.setsolver(sp, solvers.MIP)
    end

    # Set noise
    if SDDP.hasnoises(sp)
        ex = SDDP.ext(sp)
        SDDP.setnoise!(sp, ex.noises[inoise])
    else
        print("  /!\\ WARNING Using Q with a deterministic problem.\n")
    end

    # Build states vector to loop over
    if length(states) == 2
        m = SDDP.mesh(states...)
    else
        m = states[1]'
    end

    # Solve for all states
    ans = []
    for i in 1:size(m,2)
        if debug
            print(i, ", ")
        end
        setstate!(sp,m[:,i])
        JuMP.solve(sp, relaxation=relaxation, ignore_solve_hook=true)
        push!(ans, JuMP.getobjectivevalue(sp))
    end
    if debug
        println()
    end
    return ans
end

# Average over all noises
function Qtilde(sp::JuMP.Model, state::AbstractVector{Float64}; relaxation::Bool=false)
    setstate!(sp,state)
    if SDDP.hasnoises(sp)
        ex = SDDP.ext(sp)
        vs = []
        for i in 1:length(ex.noiseprobability)
            SDDP.setnoise!(sp, ex.noises[i])
            JuMP.solve(sp, relaxation=relaxation, ignore_solve_hook=true)
            push!(vs, JuMP.getobjectivevalue(sp))
        end
        return mean(vs)
    else
        JuMP.solve(sp, relaxation=relaxation, ignore_solve_hook=true)
        return JuMP.getobjectivevalue(sp)
    end
end

function Qtilde(m::SDDP.SDDPModel,t::Int,ms::Int, states::Union{Float64, AbstractVector{Float64}}...; debug::Bool=false, relaxation::Bool=false)
    sp = SDDP.getsubproblem(m,t+1,ms)
    solvers = sp.ext[:solvers]
    if relaxation
        JuMP.setsolver(sp, solvers.LP)
    else
        JuMP.setsolver(sp, solvers.MIP)
    end
    ans = []
    if length(states) == 2
        m = SDDP.mesh(states...)
    else
        m = states[1]'
    end
    for i in 1:size(m,2)
        if debug
            print(i, ", ")
        end
        push!(ans, Qtilde(sp,m[:,i], relaxation=relaxation))
    end
    if debug
        println()
    end
    return ans
end

function make_cut(m, t, x, rho)
    settings = SDDP.Settings()
    close(settings.cut_output_file)
    prev_state = SDDP.getstage(m,t-1).state
    SDDP.padvec!(prev_state, length(x))
    for (i,xi) in enumerate(x)
        prev_state[i] = xi
    end
    cut_it(m, t, settings, rho=rho)
end
