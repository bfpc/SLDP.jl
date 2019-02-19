import JuMP
import SDDP

columns(m) = [m[:,i] for i in 1:size(m,2)]

function reject_stdout(f)
    _stdout = STDOUT
    stdout_rd, stdout_wr = redirect_stdout()

    ret = try
        f()
    catch e
        println("ERROR in capture_streams(): ", e)
    finally
        # restore
        redirect_stdout(_stdout)
        close(stdout_wr)
        close(stdout_rd)
    end
    ret
end

function cutvalue(c::ALDCut, ti)
    c.v + dot(c.l, ti - c.xi) - c.rho*norm(ti - c.xi, 1)
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

function Qtilde(sp::JuMP.Model, state::AbstractVector{Float64})
    setstate!(sp,state)
    if SDDP.hasnoises(sp)
        ex = SDDP.ext(sp)
        vs = []
        for i in 1:length(ex.noiseprobability)
            SDDP.setnoise!(sp, ex.noises[i])
            reject_stdout( JuMP.solve(sp) )
            push!(vs, JuMP.getobjectivevalue(sp))
        end
        return mean(vs)
    else
        reject_stdout( JuMP.solve(sp) )
        return JuMP.getobjectivevalue(sp)
    end
end

function Qtilde(m,t,ms, states::Union{Float64, AbstractVector{Float64}}...; debug=false)
    sp = SDDP.getsubproblem(m,t+1,ms)
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
        push!(ans, Qtilde(sp,m[:,i]))
    end
    if debug
        println()
    end
    return ans
end
