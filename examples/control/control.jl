import JuMP: @variable, @constraint
import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import ASDDiP: prepareALD!, setALDsolver!

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


function controlmodel(;nstages=8, discount=0.9, ramp_mode=:None, noise=noise, shift=nothing,
                      doSB=true, tents=false)

# ==================
# Lipschitz constant
function Lip(t)
    l = (1.0 - discount^(nstages+1-t))/(1 - discount)*discount^(t-1)
    if ramp_mode == :parallel2
      return 2*l
    else
      return l
    end
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
        objective_bound   = 0,
        solver            = GurobiSolver(env,OutputFlag=0)
       ) do sp, stage

    # ======================
    # State: x-axis position
    @state(sp, -10 <= x <= 10, x0 == 2)

    # ========
    # Dynamics
    @variable(sp, 0 <= move <= 1, Int)  # Control: left or right
    # Correct to not have first stage noise
    if stage == 1
      @constraint(sp, x == x0 + 2*move - 1)
    else
      delta = (shift == nothing ? 0.0 : shift[stage])
      @rhsnoise(sp, drift=noise+delta, x == x0 + drift + 2*move - 1)
    end

    # ==========================================
    # Cost: absolute difference from equilibrium
    @variable(sp, xpos >= 0)
    @variable(sp, xneg >= 0)
    @constraint(sp, x == xpos - xneg)
    @stageobjective(sp, discount^(stage-1)*(xpos + xneg))

    # =================================
    # The solver hook for backwards ALD
    lip = Lip(stage)
    setALDsolver!(sp, lip, rho_line(lip); doSB=doSB, tents=tents)
end


# =======
# Prepare
if ramp_mode == :opt_rho
  prepareALD!(m, Lip, rho_zero)
  for stage in SDDP.stages(m)
    for sp in SDDP.subproblems(stage)
      JuMP.setsolvehook(sp, ASDDiP.ASDDiPsolve_optrho!)
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

function controlsolve(m, niters)
    @time solvestatus = SDDP.solve(m,
        iteration_limit = niters,
        simulation = sim
       )
end

