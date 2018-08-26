import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import JuMP: @variable, @constraint
import Gurobi: GurobiSolver

using ASDDiP


# =============
# \rho policies
Lip(t) = 30
function rho_ramp(niter, t, i=1)
  Lip(t) * clamp((niter-15)/15, 0, 1)
end

function rho_ramp_parallel(niter, t, i=1)
  clamp((niter-15)/15, 0, Lip(t))
end

function rho_ramp_parallel2(niter, t, i=1)
  clamp((niter-15)/15, 0, 2*Lip(t))
end

function rho_zero(niter, t, i=1)
  0
end


# ==========
# Parameters
if !isdefined(:ramp_mode)
  ramp_mode = :None
end
if !isdefined(:niters)
  niters = 0
end
if !isdefined(:num_noise)
  num_noise = 2
end

# =====
# Noise
# uniform discrete distribution over [5,15]^2
noisestep = 10/(num_noise - 1)
noisegen = 1:num_noise
noise = 5 + noisestep*[[i-1,j-1] for i = noisegen for j = noisegen]

# ==============
# Stage problems
function stage1(sp)
  x = sp[:x]
  @stageobjective(sp, -sum([3/2 4]*x))
end

function stage2(sp)
  x0 = sp[:x0]
  @variable(sp, y[1:4], Bin)
  @rhsnoise(sp, w=noise, sum([2 3 4 5]*y) <= w[1] - x0[1])
  @rhsnoise(sp, w=noise, sum([6 1 3 2]*y) <= w[2] - x0[2])
  @stageobjective(sp, -sum([16 19 23 28]*y))
end


# ==========
# SDDP Model
m = SDDPModel(
        sense             = :Min,
        stages            = 2,
        objective_bound   = -100,
        solver            = GurobiSolver(OutputFlag=0,TimeLimit=60)
       ) do sp, stage

    # ========
    # State: x
    @state(sp, 0 <= x[1:2] <= 5, x0 == 0)
    if stage == 1
      stage1(sp)
    else
      stage2(sp)
    end

    # =================================
    # The solver hook for backwards ALD
    setASDDiPsolver!(sp)
end

# ===============
# Prepare & solve
if ramp_mode == :None
  prepareALD!(m, rho_zero)
elseif ramp_mode == :simple
  prepareALD!(m, rho_ramp)
elseif ramp_mode == :parallel
  prepareALD!(m, rho_ramp_parallel)
elseif ramp_mode == :parallel2
  prepareALD!(m, rho_ramp_parallel2)
elseif ramp_mode == :opt_rho
  prepareALD!(m, rho_zero)
  for stage in SDDP.stages(m)
    for sp in SDDP.subproblems(stage)
      JuMP.setsolvehook(sp, ASDDiP.ASDDiPsolve_optrho!)
    end
  end
else
  warn("Invalid ramp mode $(ramp_mode).  Not using ALD cuts, only Strenghtened Benders.")
  prepareALD!(m, rho_zero)
end

sim = SDDP.MonteCarloSimulation(
     frequency = 10,
     min  = 100, step = 50, max  = 200,
     terminate = false
    )

if niters > 0
    @time solvestatus = SDDP.solve(m,
        iteration_limit = niters,
        simulation = sim
       )
end
