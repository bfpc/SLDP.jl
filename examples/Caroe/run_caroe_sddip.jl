include("caroe_sddip.jl")

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

# Warm-up
m = caroe_model(num_noise=2, int=false, eps=0.1, lagtol=1e-4)
caroe_solve(m,10)

# Start simulation
niter = 100
for noise in [2,3,6,11,21]
  for isint in [true,false]
    for eps in [0.1, 0.01]
      for lagtol in [1e-2, 1e-4]
        println("Simulating SDDiP with $noise-step discretization, int=$isint")
        println("                      eps = $eps, lagrangian tolerance = $lagtol")

        m = caroe_model(num_noise=noise, int=isint,
                        eps=eps, lagtol=lagtol)
        caroe_solve(m,niter)

        name = "SDDiP_" * string(eps) * "_" * string(lagtol)
        name = name * (isint ? "_int_" : "_cont_") * string(noise)
        SDDP.writecuts!(joinpath("cuts", name*"_benders.cuts"), m)
        writelog!(m, joinpath("logs", name*".log"))
        println()
      end
    end
  end
end

