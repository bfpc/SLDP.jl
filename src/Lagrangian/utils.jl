function sense(c::ConstraintRef)
    constr = c.m.linconstr[c.idx]
    if constr.lb == constr.ub
        return :eq
    elseif constr.lb == -Inf
        return :le
    elseif constr.ub == Inf
        return :ge
    end
    error("Range constraint not supported.")
end

function getslack(c::ConstraintRef)
    constr = c.m.linconstr[c.idx]
    if constr.ub == Inf
        return constr.terms - constr.lb
    else
        return constr.terms - constr.ub
    end
end

getdualsense(m::JuMP.Model) = getobjectivesense(m)==:Min? :Max : :Min

function issatisfied(slack::Float64, sense::Symbol, tol=1e-9)
    if sense == :eq
        return abs(slack) < tol
    elseif sense == :le
        return slack < tol
    elseif sense == :ge
        return slack > -tol
    end
    error("Unknown sense $(sense).")
end
