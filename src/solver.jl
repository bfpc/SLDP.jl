function default_rho_policy(niter :: Int, t :: Int, i=1 :: Int)
  return clamp((niter-10)/10, 0, 10)
end

function prepareALD!(m :: SDDP.SDDPModel, Lip::Function, rho_policy=default_rho_policy)
  for stage in SDDP.stages(m)
    sp = stage.subproblems[1]
    stage.ext[:ALDbounds] = [[st.lb, st.ub] for st in aldstates(sp)]
    stage.ext[:Lip]  = Lip(stage.t)
    stage.ext[:rhos] = rho_policy
  end
end

function SDDP.postsolve!(:: Type{SDDP.ForwardPass}, m :: SDDP.Model, sp :: JuMP.Model)
  t = SDDP.ext(sp).stage
  stage = SDDP.getstage(m, t)
  for (st,b) in zip(aldstates(sp), stage.ext[:ALDbounds])
    b[1] = st.lb
    b[2] = st.ub
  end
end

# Adjusts the inequalities to calculate |z - xin| during ALD
function set_xhat_av!(sp::JuMP.Model, xhat)
  for (st,xi) in zip(aldstates(sp), xhat)
    JuMP.setRHS(st.av_c1,  xi)
    JuMP.setRHS(st.av_c2, -xi)
  end
end

function saveInfeasible(sp :: JuMP.Model)
  filepath = joinpath(pwd(), "infeasible_subproblem.lp")
  JuMP.writeLP(sp, filepath, genericnames=false)
  error("""Model in stage $(SDDP.ext(sp).stage) and markov state $(SDDP.ext(sp).markovstate)
        was not solved to Optimality. I wrote the offending LP file to
        $(filepath).

        This is most commonly caused by numerical issues with the solver.
        Consider reformulating the model or try different solver parameters.
        """)
end

# Overload how stage problems are solved in the backward pass
# Copied from SDDiP/SDDP
function SDDP.JuMPsolve(::Type{SDDP.BackwardPass}, m::SDDP.SDDPModel, sp::JuMP.Model)
    direction = SDDP.BackwardPass
    SDDP.presolve!(direction, m, sp)

    TT = STDOUT
    _rd,_wr = redirect_stdout()
    status = try
      @timeit SDDP.TIMER "Backwards jumpsolve w/ duals" begin
        SDDP.jumpsolve(sp, require_duals=true)
      end
    finally
      redirect_stdout(TT)
      close(_rd); close(_wr)
    end

    if status != :Optimal
      saveInfeasible(sp)
    end
    SDDP.postsolve!(direction, m, sp)
end


function set_ald_objective!(sp::JuMP.Model, π)
  aff_expr = 0
  xhat = sp.ext[:ALD].xin_v
  rho  = sp.ext[:ALD].rho[1]
  set_xhat_av!(sp, xhat)
  for (st,xi,li) in zip(aldstates(sp),xhat, π)
    aff_expr += li * (xi - st.xin) + rho*st.av
  end
  sp.obj += aff_expr
end

function get_info_prevstage!(sp::JuMP.Model, m::SDDP.SDDPModel, t)
  ext = sp.ext[:ALD]
  prev_st = SDDP.getstage(m,t)
  for (i,b) in enumerate(prev_st.ext[:ALDbounds])
    ext.xin_lb[i] = b[1]
    ext.xin_ub[i] = b[2]
  end
  for (i,xi) in enumerate(prev_st.state)
    ext.xin_v[i] = xi
  end
end

function make_cut(m, t, x, rho)
    settings = SDDP.Settings()
    close(settings.cut_output_file)
    prev_state = SDDP.getstage(m,t-1).state
    SDDP.padvec!(prev_state, length(x))
    for (i,xi) in enumerate(x)
        prev_state[i] = xi
    end
    cut_it(m, t, rho, settings)
end

function cut_it(m::SDDP.SDDPModel, t::Int, rho::Real, settings::SDDP.Settings)
    # solve all stage t problems
    SDDP.reset!(m.storage)
    for sp in SDDP.subproblems(m, t)
        SDDP.setstates!(m, sp)
        sp.ext[:ALD].rho[1] = rho
        # Hack, should be done in forward pass
        ASDDiP.get_info_prevstage!(sp, m, t-1)
        empty!(sp.ext[:ALD].vstore)
        empty!(sp.ext[:ALD].lstore)
        empty!(sp.ext[:ALD].rhostore)
        SDDP.solvesubproblem!(SDDP.BackwardPass, m, sp)
    end
    # add appropriate cuts
    for sp in SDDP.subproblems(m, t-1)
        @timeit SDDP.TIMER "Cut addition" begin
            SDDP.modifyvaluefunction!(m, settings, sp)
            ASDDiP.modify_ald_valuefunction!(m, settings, sp)
        end
    end
end

# Modify SDDP.backwardpass! because we must empty! all vstore/lstore for ALD,
# over _all_ subproblems in current stage (not only the current stage)
import SDDP: backwardpass!
function backwardpass!(m::SDDP.SDDPModel, settings::SDDP.Settings)
    niter = length(m.log)
    # walk backward through the stages
    for t in SDDP.nstages(m):-1:2
        stage = SDDP.getstage(m,t)
        Lip = stage.ext[:Lip]
        rho = stage.ext[:rhos](niter, Lip)
        cut_it(m, t, rho, settings)
    end
    SDDP.reset!(m.storage)
    for sp in SDDP.subproblems(m, 1)
        SDDP.solvesubproblem!(SDDP.BackwardPass, m, sp)
    end
    return dot(m.storage.objective, m.storage.probability)
end

# The solvehook
function ASDDiPsolve!(sp::JuMP.Model; require_duals::Bool=false, kwargs...)
    solvers = sp.ext[:solvers]
    if require_duals && SDDP.ext(sp).stage > 1
        # Update the objective we cache in case the objective has noises
        l = lagrangian(sp)
        l.obj = JuMP.getobjective(sp)

        # Strengthened Augmented Benders cut:
        # Solve relaxation, then lift as possible.

        # Get the LP duals
        JuMP.setsolver(sp, solvers.LP)
        @assert JuMP.solve(sp, ignore_solve_hook=true, relaxation=true) == :Optimal
        π  = SDDP.getdual.(SDDP.states(sp))
        # Update slacks because RHSs change each iteration
        # l.slacks = getslack.(l.constraints)
        # Relax bounds to formulate Lagrangian
        Lagrangian.relaxandcache!(l, sp)
        # Change the MIP objective
        set_ald_objective!(sp, π)
        # Solve the Lagrangian, with LP πs chosen and fixed
        JuMP.setsolver(sp, solvers.MIP)
        status = JuMP.solve(sp, ignore_solve_hook=true)
        push!(sp.ext[:ALD].vstore, JuMP.getobjectivevalue(sp))
        push!(sp.ext[:ALD].lstore, π)
        push!(sp.ext[:ALD].rhostore, sp.ext[:ALD].rho[1])
        # Undo changes
        sp.obj = l.obj
        Lagrangian.recover!(l, sp, π)

        # Solve original problem again to have correct level for automatic SDDP Benders cut
        if sp.ext[:ALD].rho[1] > 0
          JuMP.setsolver(sp, solvers.LP)
          JuMP.solve(sp, ignore_solve_hook=true, relaxation=true)
        end
    else
        # We are in the forward pass, or we are in stage 1
        JuMP.setsolver(sp, solvers.MIP)
        status = JuMP.solve(sp, ignore_solve_hook=true)
        #println("Objective at stage ", SDDP.ext(sp).stage, " going forward: ", JuMP.getobjectivevalue(sp))
    end
    status
end

function bissect_rho(sp, π, mipvalue, minrho, maxrho, valuetol=1e-3, rhotol=1e-3)
  candidate = (minrho + maxrho)/2
  sp.ext[:ALD].rho[1] = candidate
  # Change the MIP objective
  set_ald_objective!(sp, π)
  #print("Solving for rho = $candidate, ")
  status = JuMP.solve(sp, ignore_solve_hook=true)
  curvalue = JuMP.getobjectivevalue(sp)
  # Go back to previous objective
  sp.obj = lagrangian(sp).obj
  #println("got $curvalue, while MIP = $mipvalue")
  if maxrho - minrho > rhotol
    if mipvalue - curvalue > valuetol
      #println("Bissecting $candidate, $maxrho")
      return bissect_rho(sp, π, mipvalue, candidate, maxrho, valuetol, rhotol)
    else
      #println("Bissecting $minrho, $candidate")
      return bissect_rho(sp, π, mipvalue, minrho, candidate, valuetol, rhotol)
    end
  else
    #println("Done: $minrho, $candidate, $maxrho, $mipvalue, $curvalue")
    return candidate, status
  end
end

# Optimal rho solvehook
function ASDDiPsolve_optrho!(sp::JuMP.Model; require_duals::Bool=false, kwargs...)
    solvers = sp.ext[:solvers]
    if require_duals && SDDP.ext(sp).stage > 1
        # Update the objective we cache in case the objective has noises
        l = lagrangian(sp)
        l.obj = JuMP.getobjective(sp)

        # Strengthened Augmented Benders cut:
        # Solve relaxation, then lift to optimal \rho

        # Get optimal MIP value
        println("Solving stage $(SDDP.ext(sp).stage) backwards")
        @assert JuMP.solve(sp, ignore_solve_hook=true) == :Optimal
        mipvalue = JuMP.getobjectivevalue(sp)
        #println("Got MIP value = $mipvalue")

        # Get the LP duals
        JuMP.setsolver(sp, solvers.LP)
        @assert JuMP.solve(sp, ignore_solve_hook=true, relaxation=true) == :Optimal
        curvalue = JuMP.getobjectivevalue(sp)
        #println("Got LP value = $curvalue")
        π  = SDDP.getdual.(SDDP.states(sp))

        if curvalue == mipvalue
          optrho = 0
          status = :Optimal
          println(optrho, " ", status)
        else
          # Update slacks because RHSs change each iteration
          # l.slacks = getslack.(l.constraints)
          # Relax bounds to formulate Lagrangian
          Lagrangian.relaxandcache!(l, sp)
          JuMP.setsolver(sp, solvers.MIP)
          optrho,status = bissect_rho(sp, π, mipvalue, 0, sp.ext[:ALD].Lip)
          println(optrho, " ", status)
          push!(sp.ext[:ALD].vstore, JuMP.getobjectivevalue(sp))
          push!(sp.ext[:ALD].lstore, π)
          push!(sp.ext[:ALD].rhostore, optrho)

          # Undo changes
          Lagrangian.recover!(l, sp, π)
          # Solve original problem again to have correct level for automatic SDDP Benders cut
          JuMP.setsolver(sp, solvers.LP)
          JuMP.solve(sp, ignore_solve_hook=true, relaxation=true)
        end
    else
        # We are in the forward pass, or we are in stage 1
        JuMP.setsolver(sp, solvers.MIP)
        status = JuMP.solve(sp, ignore_solve_hook=true)
        #println("Objective at stage ", SDDP.ext(sp).stage, " going forward: ", JuMP.getobjectivevalue(sp))
    end
    status
end

"""
    setSDDiPsolver!(sp::JuMP.Model; method=Subgradient(0.0), pattern=Pattern(), MIPsolver=sp.solver, LPsolver=MIPsolver)

Sets a JuMP solvehook for integer SDDP to stage problem `sp` that will call a
a Lagrangian solver of type `method.` Argument `pattern` can be used to specify
a pattern of different cut types.

You should specify an LP/MIP solver if you are using different cut types in a cut
pattern, and you are not using a solver that can solve both MIPs and LPs.
"""
function setASDDiPsolver!(sp::JuMP.Model; MIPsolver=sp.solver, LPsolver=MIPsolver)
    n = length(SDDP.states(sp))
    sp.ext[:ALD] = ALDExtension([],10,zeros(n),zeros(n),zeros(n),[],[],[],[1.5],[])

    constraints = SDDP.LinearConstraint[]
    for s in SDDP.states(sp)
        # xinₜ = xoutₜ₋₁ is being relaxed
        push!(constraints, s.constraint)
        s0 = sp.linconstr[s.constraint.idx].terms.vars[1]
        ald_statevariable!(sp, s0, s.variable, s.constraint)
    end


    sp.ext[:Lagrangian] = Lagrangian.LinearProgramData(JuMP.QuadExpr(),         # objective
                                           constraints)         # relaxed constraints
    sp.ext[:solvers] = MixedSolvers(LPsolver, MIPsolver)

    JuMP.setsolvehook(sp, ASDDiPsolve!)
end
