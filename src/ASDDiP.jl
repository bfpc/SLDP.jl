module ASDDiP

export setALDsolver!, prepareALD!

# Our Lagrangian solver
include(joinpath(dirname(@__FILE__),"Lagrangian","Lagrangian.jl"))

import SDDP
using JuMP, Reexport
using TimerOutputs: @timeit

@reexport using .Lagrangian

include("typedefs.jl")
include("solver.jl")
include("state.jl")
include("valuefunction.jl")
include("setup.jl")

include("utils.jl")
end
