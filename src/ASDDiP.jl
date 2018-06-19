module ASDDiP

export setASDDiPsolver!, prepareALD!

# Our Lagrangian solver
include(joinpath(dirname(@__FILE__),"Lagrangian","Lagrangian.jl"))

using SDDP, JuMP, Reexport
import TimerOutputs: @timeit

@reexport using .Lagrangian

include("typedefs.jl")
include("solver.jl")
include("state.jl")
include("valuefunction.jl")

include("utils.jl")
end
