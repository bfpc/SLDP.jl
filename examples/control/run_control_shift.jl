include("control.jl")

#
# Simulation configs
# noise & discount; names, ramp_modes & #iter
include("constants.jl")
using Const

#
# Start simulation
for (basename,ramp) in [("SB", :None)]
  niter = 100
  for sb in (true)
    name = basename * "_shift"
    println("Simulating $basename, using SB cuts")
    m = controlmodel(nstages=8, discount=discount, ramp_mode=ramp, noise=noise, shift=shift,
                    doSB=sb)
    controlsolve(m,niter)

    SDDP.writecuts!(joinpath("cuts", name*"_benders.cuts"), m)
    ASDDiP.save_aldcuts!(m, joinpath("cuts", name*"_ald.cuts"))
    ASDDiP.writelog!(m, joinpath("logs", name*".log"))
  end
end
for (basename,ramp) in [("ALD simple", :simple), ("ALD parallel", :parallel)]
  niter = 200
  for tent in (true, false)
    for sb in (true, false)
      name = basename * (tent ? "_tent" : "_l1") * (sb ? "_SB" : "_B") * "_shift"
      println("Simulating $basename, $(tent ? "with tents" : "with general cuts"), using $(sb? "SB" : "Benders") cuts")
      m = controlmodel(nstages=8, discount=discount, ramp_mode=ramp, noise=noise, shift=shift,
                      doSB=sb, tents=tent)
      controlsolve(m,niter)

      SDDP.writecuts!(joinpath("cuts", name*"_benders.cuts"), m)
      ASDDiP.save_aldcuts!(m, joinpath("cuts", name*"_ald.cuts"))
      ASDDiP.writelog!(m, joinpath("logs", name*".log"))
    end
  end
end

