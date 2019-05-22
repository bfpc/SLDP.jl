include("control.jl")

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
  SLDP.save_aldcuts!(m, joinpath("cuts", name*"_ald.cuts"))
  SLDP.writelog!(m, joinpath("logs", name*".log"))
end

