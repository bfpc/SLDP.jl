include("caroe.jl")

# Simulation configs
run_names  = ["SB", "ALD parallel"]
ramp_modes = [:None, :parallel]
iters      = [100, 200]

#
# Start simulation
for (basename,ramp,niter) in zip(run_names, ramp_modes, iters)
  for noise in [2,3,6] #,11,21]
    for isint in [true,false]
      println("Simulating $basename with $noise discretization, int=$isint")
      m = caroe_model(ramp_mode=ramp, num_noise=noise, doSB=true, int=isint)
      caroe_solve(m,niter)

      name = basename * (isint ? "_int_" : "_cont_") * string(noise)
      SDDP.writecuts!(joinpath("cuts", name*"_benders.cuts"), m)
      ASDDiP.save_aldcuts!(m, joinpath("cuts", name*"_ald.cuts"))
      ASDDiP.writelog!(m, joinpath("logs", name*".log"))
    end
  end
end

