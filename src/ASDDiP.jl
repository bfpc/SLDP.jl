module ASDDiP

export setASDDiPsolver!, prepareALD!

# Our Lagrangian solver
include(joinpath(dirname(@__FILE__),"Lagrangian","Lagrangian.jl"))

using SDDP, JuMP, Reexport
@reexport using .Lagrangian

include("typedefs.jl")
include("solver.jl")
include("state.jl")

include("utils.jl")
end
