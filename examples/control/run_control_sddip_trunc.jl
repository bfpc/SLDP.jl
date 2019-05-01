include("control_sddip.jl")

# Auxiliary functions to save logs
function writelogline!(line, file)
    v = (line.iteration, line.bound, line.lower_statistical_bound, line.upper_statistical_bound,
         line.timecuts, line.simulations, line.timesimulations, line.timetotal)
    print(file, v[1])
    for vi in v[2:end]
        print(file, ",", vi)
    end
    print(file, "\n")
end

function writelog!(m::SDDP.SDDPModel, filename)
    open(filename, "w") do file
        for line in m.log
            writelogline!(line, file)
        end
    end
end

#
# Simulation configs
# noise & discount; names, ramp_modes & #iter
include("constants.jl")
using Const


# Start simulation
niter = 100
for eps in [0.1, 0.01]
  for lagtol in [1e-2, 1e-4]
    name = "SDDiP_" * string(eps) * "_" * string(lagtol)
    println("Simulating SDDiP, eps = $eps, lagrangian tolerance = $lagtol")
    epsnoise = eps*round.(noise/eps)
    m = controlmodel(nstages=8, discount=discount, noise=epsnoise,
                     eps=eps, lagtol=lagtol)
    controlsolve(m,niter)

    SDDP.writecuts!(joinpath("cuts", name*"_benders.cuts"), m)
    writelog!(m, joinpath("logs", name*".log"))
  end
  println()
end
