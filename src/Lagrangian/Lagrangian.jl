module Lagrangian
using JuMP

export
# Each method needs to know about the problem being solved
LinearProgramData

type LinearProgramData
    obj::QuadExpr                               # objective
    constraints::Vector{<:ConstraintRef}        # constraints being relaxed
    senses::Vector{Symbol}                      # we will cache the sense of constraints
    slacks::Vector{AffExpr}                     # also cache Ax-b
    old_bound::Vector{Float64}                  # cache before relaxing constraints
end

include("utils.jl")

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
        getslack.(constraints),
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
