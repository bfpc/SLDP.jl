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

function get_info_prevstage!(sp::JuMP.Model, m::SDDP.SDDPModel)
  ext = sp.ext[:ALD]
  t = sp.ext[:SDDP].stage
  prev_st = SDDP.getstage(m,t-1)
  for (i,b) in enumerate(prev_st.ext[:ALDbounds])
    ext.xin_lb[i] = b[1]
    ext.xin_ub[i] = b[2]
  end
  for (i,xi) in enumerate(prev_st.state)
    ext.xin_v[i] = xi
  end
end

function cut_it(m::SDDP.SDDPModel, t::Int, settings::SDDP.Settings; rho=nothing::Union{Void,Real})
    # solve all stage t problems
    niter = length(m.log)
    for sp in SDDP.subproblems(m, t)
        # Set xin
        SDDP.setstates!(m, sp)

        # Set rho for ALD
        if rho == nothing
          lip = aldparams(sp).Lip
          a,b = aldparams(sp).rho_line
          rho = clamp(a*niter + b, 0.0, lip)
        end
        sp.ext[:ALD].rho[1] = rho

        # Update ALD equations for |z - xin|
        get_info_prevstage!(sp, m)

        # Solve all noises, and store each cut information
        empty!(sp.ext[:ALD].vstore)
        empty!(sp.ext[:ALD].lstore)
        empty!(sp.ext[:ALD].rhostore)
        SDDP.solvesubproblem!(SDDP.BackwardPass, m, sp)
    end
    # add appropriate cuts
    for sp in SDDP.subproblems(m, t-1)
        @timeit SDDP.TIMER "Cut addition" begin
            SDDP.modifyvaluefunction!(m, settings, sp)
            modify_ald_valuefunction!(m, settings, sp)
        end
    end
end

# Modify SDDP.backwardpass! because we must empty! all vstore/lstore for ALD,
# over _all_ subproblems in current stage (not only the current markov state)
import SDDP: backwardpass!
function backwardpass!(m::SDDP.SDDPModel, settings::SDDP.Settings)
    # walk backward through the stages
    for t in SDDP.nstages(m):-1:2
        SDDP.reset!(m.storage)
        cut_it(m, t, settings)
    end
    SDDP.reset!(m.storage)
    for sp in SDDP.subproblems(m, 1)
        SDDP.solvesubproblem!(SDDP.BackwardPass, m, sp)
    end
    return dot(m.storage.objective, m.storage.probability)
end

function getpi_SB(sp::JuMP.Model)
  # If using tents, just return zeroes
  if aldparams(sp).tents
    return zeros(len(SDDP.states(sp)))
  else
    # π via relaxation
    lpsolver = sp.ext[:solvers].LP
    JuMP.setsolver(sp, lpsolver)
    @assert JuMP.solve(sp, ignore_solve_hook=true, relaxation=true) == :Optimal
    return SDDP.getdual.(SDDP.states(sp))
  end
end

# The solvehook
function ALDsolve!(sp::JuMP.Model; require_duals::Bool=false, kwargs...)
    solvers = sp.ext[:solvers]
    if require_duals && SDDP.ext(sp).stage > 1
        # Update the objective we cache in case the objective has noises
        l = lagrangian(sp)
        l.obj = JuMP.getobjective(sp)

        # Strengthened Augmented Benders cut: lift as much as possible given π and rho
        π = getpi_SB(sp)
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
          # Maybe insert here a SB cut, not just Benders cut
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

