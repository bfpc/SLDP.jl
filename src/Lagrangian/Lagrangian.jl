# This file includes modified source code from https://github.com/lkapelevich/SDDiP.jl
# commit a59e55ba227533a399a6b535270f6ce1ba519862
#############################################################################

module Lagrangian
using JuMP

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

# Each method needs to know about the problem being solved
type LinearProgramData
    obj::QuadExpr                               # objective
    constraints::Vector{<:ConstraintRef}        # constraints being relaxed
    senses::Vector{Symbol}                      # we will cache the sense of constraints
    old_bound::Vector{Float64}                  # cache before relaxing constraints
end

"""
    LinearProgramData(obj::QuadExpr, constraints::Vector{<:ConstraintRef})

Creates a `LinearProgramData` object for calling `lagrangiansolve!`.

# Arguments
* obj:               Objective of the linear program without relaxation.
* constraints:       A vector of type `JuMP.ConstraintRef` with the contraints to be relaxed.
"""
function LinearProgramData(obj::QuadExpr, constraints::Vector{<:ConstraintRef})
    LinearProgramData(obj, constraints,
        sense.(constraints),
        zeros(Float64, length(constraints)),
    )
end

function relaxandcache!(l::LinearProgramData, m::JuMP.Model)
    ext = m.ext[:ALD]
    for (i, c) in enumerate(l.constraints)
        con = m.linconstr[c.idx]
        if l.senses[i] == :le
            l.old_bound[i] = con.ub
            con.ub = ext.xin_ub[i]
        elseif l.senses[i] == :ge
            l.old_bound[i] = con.lb
            con.lb = ext.xin_lb[i]
        else
            l.old_bound[i] = con.lb
            # Set to infinity because of Gurobi not handling RangeConstraints
            con.ub = Inf
            con.lb = ext.xin_lb[i]
        end
    end
end

function recover!(l::LinearProgramData, m::JuMP.Model, π::Vector{Float64})
    for (i, c) in enumerate(l.constraints)
        con = m.linconstr[c.idx]
        m.linconstrDuals[c.idx] = π[i]
        if l.senses[i] == :le
            con.ub = l.old_bound[i]
        elseif l.senses[i] == :ge
            con.lb = l.old_bound[i]
        else
            con.lb = con.ub = l.old_bound[i]
        end
    end
end

end # module
