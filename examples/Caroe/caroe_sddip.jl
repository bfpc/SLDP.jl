import JuMP: @variable, @constraint
import SDDP: SDDPModel, @state, @rhsnoise, @stageobjective
import SDDiP: @binarystate, setSDDiPsolver!, Pattern, SubgradientMethod, LevelMethod

import Gurobi: GurobiSolver


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


function caroe_model(;nstages=2, num_noise=2, int=false, eps=0.1, pattern=(0,1,1), lagtol=1e-4)

# =====
# Noise
# uniform discrete distribution over [5,15]^2
noisestep = 10/(num_noise - 1)
noisegen = 1:num_noise
noise = 5 + noisestep*[[i-1,j-1] for i = noisegen for j = noisegen]

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
    @binarystate(sp, 0 <= x[1:2] <= 5, x0 == 0, Cont, eps)
    if stage == 1
      stage1(sp, int)
    else
      stage2(sp, noise)
    end

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

function caroe_solve(m, niters)
    @time solvestatus = SDDP.solve(m,
        iteration_limit = niters,
        simulation = sim
       )
end
