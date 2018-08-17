import JuMP
import SDDP

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

function estimatevf(sp::JuMP.Model, states...)
    vf = []
    ac = ASDDiP.aldcuts(sp)
    if length(ac) == 0
        b = SDDP.ext(sp).problembound
        if length(states) == 1
            npts = length(states[1])
        else
            npts = length(states[1])*length(states[2])
        end
        return b * ones(npts)
    end
    for c in ac
        if length(states) == 1
           ts = states[1]
        else
            ts = SDDP.mesh(states...)
        end
        vf_i = [c.v + dot(c.l, ti - c.xi) - c.rho*norm(ti - c.xi, 1) for ti in ts]
        push!(vf, vf_i)
    end
    maximum(hcat(vf...), 2)[:,1]
end

function Qfrak(m,t,ms, states::Union{Float64, AbstractVector{Float64}}...)
    sp = SDDP.getsubproblem(m,t,ms)
    vf = SDDP.valueoracle(sp)
    _, y_cont = SDDP.processvaluefunctiondata(vf, true, states...)
    y_ald = estimatevf(sp, states...)
    maximum(hcat(y_cont, y_ald), 2)[:,1]
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
