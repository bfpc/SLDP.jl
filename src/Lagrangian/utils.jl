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

getdualsense(m::JuMP.Model) = getobjectivesense(m)==:Min? :Max : :Min
