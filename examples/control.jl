import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import JuMP: @variable, @constraint
import Gurobi: GurobiSolver

using ASDDiP


# =============
# \rho policies
Lip(t) = (1 - discount^(8+1-t))/(1 - discount)*discount^(t-1)
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
if !isdefined(:discount)
  discount = 0.9
end
if !isdefined(:ramp_mode)
  ramp_mode = :None
end
if !isdefined(:niters)
  niters = 0
end

# random noise
srand(11111)
A_noise   = 0.4
num_noise = 10
noise = A_noise * randn(num_noise)


# ==========
# SDDP Model
m = SDDPModel(
        sense             = :Min,
        stages            = 8,
        objective_bound   = 0,
        solver            = GurobiSolver(OutputFlag=0)
                                                ) do sp, stage

    # ======================
    # State: x-axis position
    @state(sp, -10 <= x <= 10, x0 == 2)

    # ========
    # Dynamics
    @variable(sp, 0 <= move <= 1, Int)  # Control: left or right
    @rhsnoise(sp, drift=noise, x == x0 + drift + 2*move - 1)

    # ==========================================
    # Cost: absolute difference from equilibrium
    @variable(sp, xpos >= 0)
    @variable(sp, xneg >= 0)
    @constraint(sp, x == xpos - xneg)
    @stageobjective(sp, discount^(stage-1)*(xpos + xneg))

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
else
  warn("Invalid ramp mode $(ramp_mode).  Not using ALD cuts, only Strenghtened Benders.")
  prepareALD!(m, rho_zero)
end

sim = SDDP.MonteCarloSimulation(
     frequency = 50,
     min  = 100, step = 50, max  = 200,
     terminate = false
    )

if niters == 0
  return
end

@time solvestatus = SDDP.solve(m,
    iteration_limit = niters,
    simulation = sim
   )
