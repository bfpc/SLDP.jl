import JuMP: @variable, @constraint
import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import SDDiP: @binarystate, setSDDiPsolver!, Pattern, SubgradientMethod, LevelMethod

import Gurobi: GurobiSolver


function controlmodel(;nstages=8, discount=0.9, noise=noise, eps=0.1, lagtol=1e-4)

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
    # Correct to not have first stage noise
    if stage == 1
      @constraint(sp, x == x0 + 2*move - 1)
    else
      @rhsnoise(sp, drift=noise, x == x0 + drift + 2*move - 1)
    end

    # ==========================================
    # Cost: absolute difference from equilibrium
    @variable(sp, xpos >= 0)
    @variable(sp, xneg >= 0)
    @constraint(sp, x - 10 == xpos - xneg)
    @stageobjective(sp, discount^(stage-1)*(xpos + xneg))

    # =================================
    # The solver hook for backwards SDDiP
    pat  = Pattern(benders=0, strengthened_benders=1, lagrangian=1)
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
