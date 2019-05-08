import JuMP: @variable, @constraint
import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import SDDiP: @binarystate, setSDDiPsolver!, Pattern, SubgradientMethod, LevelMethod

import Gurobi: GurobiSolver


function controlmodel(;nstages=8, discount=0.9, noise=noise, eps=0.1, pattern=(0,1,1), lagtol=1e-4)

# ==================
# Lipschitz constant
function Lip(t)
    return (1.0 - discount^(nstages+1-t))/(1 - discount)*discount^(t-1)
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
    @binarystate(sp, 0 <= x <= 20, x0 == 12, Cont, eps)

    # ========
    # Dynamics
    @variable(sp, 0 <= move <= 1, Int)  # Control: left or right
    @variable(sp, s) # slack
    # Correct to not have first stage noise
    if stage == 1
      @constraint(sp, x == x0 + 2*move - 1 + s)
    else
      @rhsnoise(sp, drift=noise, x == x0 + drift + 2*move - 1 + s)
    end

    # ====================================================================
    # Cost: absolute difference from equilibrium, plus integrality penalty
    @variable(sp, xpos >= 0)
    @variable(sp, xneg >= 0)
    @constraint(sp, x - 10 == xpos - xneg)

    @variable(sp, norm_s)
    @constraint(sp, norm_s >= s)
    @constraint(sp, norm_s >= -s)
    @stageobjective(sp, discount^(stage-1)*(xpos + xneg) + (1+Lip(stage))*norm_s)

    # =================================
    # The solver hook for backwards SDDiP
    pat  = Pattern(benders=pattern[1], strengthened_benders=pattern[2], lagrangian=pattern[3])
    meth = LevelMethod(initialbound=0.0, tol=SDDiP.Unit(lagtol), quadsolver=GurobiSolver(env,OutputFlag=0), maxit=1_000)
    setSDDiPsolver!(sp, method=meth, pattern=pat)
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
