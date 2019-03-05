include("control.jl")

#
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
            # print(file, line)
        end
    end
end

#
# Simulation configs
# noise & discount; names, ramp_modes & #iter
include("constants.jl")
using Const

#
# Start simulation
for (name,ramp,niter) in zip(run_names, ramp_modes, iters)
  println("Simulating $name")
  m = controlmodel(nstages=8, discount=discount, ramp_mode=ramp, noise=noise)
  controlsolve(m,niter)

  SDDP.writecuts!(joinpath("cuts", name*"_benders.cuts"), m)
  ASDDiP.save_aldcuts!(m, joinpath("cuts", name*"_ald.cuts"))
  writelog!(m, joinpath("logs", name*".log"))
end

