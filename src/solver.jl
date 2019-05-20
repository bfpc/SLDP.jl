# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
function SDDP.JuMPsolve(direction::Type{SDDP.BackwardPass}, m::SDDP.SDDPModel, sp::JuMP.Model)
    SDDP.presolve!(direction, m, sp)

    @timeit SDDP.TIMER "Backwards jumpsolve w/ duals" begin
      status = SDDP.jumpsolve(sp, require_duals=true)
    end

    if status != :Optimal
      saveInfeasible(sp)
    end
    SDDP.postsolve!(direction, m, sp)
end


function set_sb_objective!(sp::JuMP.Model, π)
  aff_expr = 0
  xhat = sp.ext[:ALD].xin_v
  for (st,xi,li) in zip(aldstates(sp),xhat, π)
    aff_expr += li * (xi - st.xin)
  end
  sp.obj += aff_expr
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
          if aldparams(sp).tents
            rho = lip
          else
            a,b = aldparams(sp).rho_line
            rho = clamp(a*niter + b, 0.0, lip)
          end
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
            if !aldparams(sp).tents
                SDDP.modifyvaluefunction!(m, settings, sp)
            end
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
    return zeros(length(SDDP.states(sp)))
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
        # tents are "easy"
        if aldparams(sp).tents
            π = zeros(length(SDDP.states(sp)))
            JuMP.setsolver(sp, solvers.MIP)
            status = JuMP.solve(sp, ignore_solve_hook=true)
            push!(sp.ext[:ALD].vstore, JuMP.getobjectivevalue(sp))
            push!(sp.ext[:ALD].lstore, π)
            push!(sp.ext[:ALD].rhostore, sp.ext[:ALD].rho[1])
            return status
        end

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

        # Solve problem again to have correct RHS in the SDDP linear cut (Benders / Strenghtened Benders)
        if sp.ext[:ALD].rho[1] > 0
          sp.obj = l.obj
          if aldparams(sp).useSB
            set_sb_objective!(sp, π)
            JuMP.solve(sp, ignore_solve_hook=true)
          else
            JuMP.setsolver(sp, solvers.LP)
            JuMP.solve(sp, ignore_solve_hook=true, relaxation=true)
          end
        end

        # Undo changes
        sp.obj = l.obj
        Lagrangian.recover!(l, sp, π)
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
  status = JuMP.solve(sp, ignore_solve_hook=true)
  curvalue = JuMP.getobjectivevalue(sp)
  # Go back to previous objective
  sp.obj = lagrangian(sp).obj
  if maxrho - minrho > rhotol
    if mipvalue - curvalue > valuetol
      return bissect_rho(sp, π, mipvalue, candidate, maxrho, valuetol, rhotol)
    else
      return bissect_rho(sp, π, mipvalue, minrho, candidate, valuetol, rhotol)
    end
  else
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
        JuMP.setsolver(sp, solvers.MIP)
        @assert JuMP.solve(sp, ignore_solve_hook=true) == :Optimal
        mipvalue = JuMP.getobjectivevalue(sp)

        # Get the LP duals
        JuMP.setsolver(sp, solvers.LP)
        @assert JuMP.solve(sp, ignore_solve_hook=true, relaxation=true) == :Optimal
        curvalue = JuMP.getobjectivevalue(sp)
        π  = SDDP.getdual.(SDDP.states(sp))

        # Add valuetol as parameter to use here and in bissection
        if curvalue == mipvalue
          optrho = 0
          push!(sp.ext[:ALD].vstore, JuMP.getobjectivevalue(sp))
          push!(sp.ext[:ALD].lstore, π)
          push!(sp.ext[:ALD].rhostore, optrho)
          status = :Optimal
        else
          # Update slacks because RHSs change each iteration
          # l.slacks = getslack.(l.constraints)
          # Relax bounds to formulate Lagrangian
          Lagrangian.relaxandcache!(l, sp)
          JuMP.setsolver(sp, solvers.MIP)
          lip = aldparams(sp).Lip
          optrho,status = bissect_rho(sp, π, mipvalue, 0, lip)
          push!(sp.ext[:ALD].vstore, JuMP.getobjectivevalue(sp))
          push!(sp.ext[:ALD].lstore, π)
          push!(sp.ext[:ALD].rhostore, optrho)

          # Solve problem again to have correct RHS in the SDDP linear cut (Benders / Strenghtened Benders)
          if aldparams(sp).useSB
            set_sb_objective!(sp, π)
            JuMP.solve(sp, ignore_solve_hook=true)
            # Undo changes
            sp.obj = l.obj
            Lagrangian.recover!(l, sp, π)
          else
            JuMP.setsolver(sp, solvers.LP)
            JuMP.solve(sp, ignore_solve_hook=true, relaxation=true)
          end
        end
    else
        # We are in the forward pass, or we are in stage 1
        JuMP.setsolver(sp, solvers.MIP)
        status = JuMP.solve(sp, ignore_solve_hook=true)
        #println("Objective at stage ", SDDP.ext(sp).stage, " going forward: ", JuMP.getobjectivevalue(sp))
    end
    status
end

