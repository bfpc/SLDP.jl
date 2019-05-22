# Copyright Bernardo Freitas Paulo da Costa 2019
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function prepareALD!(m :: SDDP.SDDPModel)
  for stage in SDDP.stages(m)
    sp = stage.subproblems[1]
    stage.ext[:ALDbounds] = [[st.lb, st.ub] for st in aldstates(sp)]
  end
end


"""
    setALDsolver!(sp::JuMP.Model, lip::Float64, rho_line::Tuple{Float64,Float64};
                  tents::Bool=false, doSB::Bool=true, maxcuts=0, drop=0,
                  MIPsolver=sp.solver, LPsolver=MIPsolver)

Sets a JuMP solvehook for integer SDDP to stage problem `sp` that will use non-convex cuts,
with maximum Lipschitz constant `lip`.

If tents=true, then we use reverse 1-norm cuts, with Lipschitz constant lip.

If tents=false, then we use strengthened augmented Benders cuts.
The value of the Lagrangian augmentation parameter rho is given by a*niter + b,
where (a,b) = `rho_line`.
If doSB=true, then also add a linear strengthened Benders cut;
if doSB=false, then also add a linear Benders cut.
Generally, one should use doSB=true.

You should specify an LP/MIP solver if you are using SA-Benders cuts,
and your solver can't solve both MIPs and LPs.


For future use, we have two extra parameters.
The solver will store at most `maxcuts` ALD cuts (0 = infinite),
and will remove ALD cuts according to `drop`.
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
