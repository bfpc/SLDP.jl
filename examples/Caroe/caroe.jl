import JuMP: @variable, @constraint
import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import SLDP: prepareALD!, setALDsolver!, ALDsolve_optrho!

import Gurobi: GurobiSolver


# =============
# \rho policies
function rho_ramp(Lip::Float64)
  Lip/15, 0.
end

function rho_ramp_parallel(Lip::Float64)
  1/15, 0.
end

function rho_ramp_parallel2(Lip::Float64)
  1/15, 0.
end

function rho_zero(Lip::Float64)
  0.,0.
end

# And dictionnary
rho_fun = Dict(
  :None => rho_zero,
  :simple => rho_ramp,
  :parallel => rho_ramp_parallel,
  :parallel2 => rho_ramp_parallel2,
  :opt_rho => rho_zero
 )


module MyData
cost1 = [3/2 4]
cost2 = [16 19 23 28]
end

# ==============
# Stage problems
function stage1(sp, int)
  x = sp[:x]
  if int
    @variable(sp, z[1:2], Int)
    @constraint(sp, z .== x)
  end
  @stageobjective(sp, -sum(MyData.cost1*x))
end

function stage2(sp, noise)
  x0 = sp[:x0]
  @variable(sp, y[1:4], Bin)
  @rhsnoise(sp, w=noise, sum([2 3 4 5]*y) <= w[1] - x0[1])
  @rhsnoise(sp, w=noise, sum([6 1 3 2]*y) <= w[2] - x0[2])
  @stageobjective(sp, -sum(MyData.cost2*y))
end


function caroe_model(;nstages=2, ramp_mode=:None, num_noise=2, int=false,
                      doSB=true, tents=false)

# =====
# Noise
# uniform discrete distribution over [5,15]^2
noisestep = 10/(num_noise - 1)
noisegen = 1:num_noise
noise = 5 + noisestep*[[i-1,j-1] for i = noisegen for j = noisegen]

# ==================
# Lipschitz constant
function Lip(t)
    return 30.
end

rho_line = try rho_fun[ramp_mode]
catch
  warn("Invalid ramp mode $(ramp_mode).  Not using ALD cuts, only Strenghtened Benders.")
  rho_zero
end


# ======
# Solver
env = Gurobi.Env()

# ==========
# SDDP Model
m = SDDPModel(
        sense             = :Min,
        stages            = nstages,
        objective_bound   = -(sum(MyData.cost1)*5 + sum(MyData.cost2)*(nstages-1)),
        solver            = GurobiSolver(env,OutputFlag=0,TimeLimit=60)
       ) do sp, stage

    # ========
    # State: x
    @state(sp, 0 <= x[1:2] <= 5, x0 == 0)
    if stage == 1
      stage1(sp, int)
    else
      stage2(sp, noise)
    end

    # =================================
    # The solver hook for backwards ALD
    lip = Lip(stage)
    setALDsolver!(sp, lip, rho_line(lip); doSB=doSB, tents=tents)
end

# =======
# Prepare
if ramp_mode == :opt_rho
  prepareALD!(m)
  for stage in SDDP.stages(m)
    for sp in SDDP.subproblems(stage)
      JuMP.setsolvehook(sp, ALDsolve_optrho!)
    end
  end
else
  prepareALD!(m)
end

return m
end

sim = SDDP.MonteCarloSimulation(
     frequency = 50,
     min  = 100, step = 50, max  = 200,
     terminate = false
    )

function caroe_solve(m, niters)
    @time solvestatus = SDDP.solve(m,
        iteration_limit = niters,
        simulation = sim
       )
end
