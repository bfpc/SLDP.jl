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

#
# Functions to calculate the Future cost functions
#

# Qfrak, the current estimate, from both linear and ALD cuts
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

function Qfrak(m,t,ms, states::Union{Float64, AbstractVector{Float64}}...; ncuts=nothing, nald=nothing)
    sp = SDDP.getsubproblem(m,t,ms)
    vf = SDDP.valueoracle(sp)

    ac = aldcuts(sp)
    cuts = SDDP.valid_cuts(SDDP.cutoracle(vf))
    if ncuts != nothing
      cuts = cuts[1:ncuts]
    end
    if nald != nothing
      ac = ac[1:nald]
    end

    y_cont = Float64[]
    y_ald  = Float64[]
    if length(cuts) > 0
      _, y_cont = SDDP.processvaluefunctiondata(cuts, true, states...)
    end
    if length(ac) > 0
        y_ald = estimatevf(ac, states...)
    end

    if length(cuts) == 0
        return y_ald
    else
        if length(ac) == 0
            return y_cont
        else
            return maximum(hcat(y_cont, y_ald), 2)[:,1]
        end
    end
end

#
# Other ALD-related utils
#

# Solve all instances in stage t to produce a cut for the previous stage
function make_cut(m, t, x, rho)
    settings = SDDP.Settings()
    close(settings.cut_output_file)
    prev_state = SDDP.getstage(m,t-1).state
    SDDP.padvec!(prev_state, length(x))
    for (i,xi) in enumerate(x)
        prev_state[i] = xi
    end
    SDDP.reset!(m.storage)
    cut_it(m, t, settings, rho=rho)
end
