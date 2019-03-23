function prepareALD!(m :: SDDP.SDDPModel)
  for stage in SDDP.stages(m)
    sp = stage.subproblems[1]
    stage.ext[:ALDbounds] = [[st.lb, st.ub] for st in aldstates(sp)]
  end
end


"""
    setALDsolver!(sp::JuMP.Model, lip::Float64, rho_line::Tuple{Float64,Float64};
                  tents::Bool=false, maxcuts=0, drop=0, MIPsolver=sp.solver, LPsolver=MIPsolver)

Sets a JuMP solvehook for integer SDDP to stage problem `sp` that will use ALD cuts,
with maximum lipschitz constant `lip`.
The value of the lagrangian augmentation parameter rho is given by a*niter + b,
where (a,b) = `rho_line`.

The solver will store at most `maxcuts` ALD cuts (0 = infinite),
and will remove ALD cuts according to `drop`.

You should specify an LP/MIP solver if you are using different cut types in a cut
pattern, and you are not using a solver that can solve both MIPs and LPs.
"""
function setALDsolver!(sp::JuMP.Model, lip::Float64, rho_line::Tuple{Float64,Float64};
                       doSB::Bool=true, tents::Bool=false, maxcuts=0, drop=0,
                       MIPsolver=sp.solver, LPsolver=MIPsolver)
    n = length(SDDP.states(sp))
    params = ALDparams(lip,rho_line,doSB,tents,maxcuts,drop)
    sp.ext[:ALD] = ALDExtension([],params,
                                zeros(n),zeros(n),zeros(n),[0.0],
                                [],[],[],
                                [])

    # Build list of "dynamic constraints", and setup ALD equations for |z - xin|
    constraints = SDDP.LinearConstraint[]
    for s in SDDP.states(sp)
        # xinₜ = xoutₜ₋₁ is being relaxed
        push!(constraints, s.constraint)
        s0 = sp.linconstr[s.constraint.idx].terms.vars[1]
        ald_statevariable!(sp, s0, s.variable, s.constraint)
    end

    # Save a pointer to the dynamic constraints
    # The objective could change in different scenarios, so don't store it just now
    sp.ext[:Lagrangian] = Lagrangian.LinearProgramData(JuMP.QuadExpr(), constraints)

    # Setup solver for the backwardpass
    sp.ext[:solvers] = MixedSolvers(LPsolver, MIPsolver)
    JuMP.setsolvehook(sp, ALDsolve!)
end
