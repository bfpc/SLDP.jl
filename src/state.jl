# This file includes modified source code from https://github.com/lkapelevich/SDDiP.jl
# commit a59e55ba227533a399a6b535270f6ce1ba519862
# which itself comes from https://github.com/odow/SDDP.jl
# commit 5b3ba3910f6347765708dd6e7058e2fcf7d13ae5

#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

bound_error(x) =  error("You must provide both upper and lower bounds on $(x).")

# This is called by setASDDiPsolver!
function ald_statevariable!(m::JuMP.Model, xin::JuMP.Variable, xout::JuMP.Variable,
                            cc)
    xin0 = JuMP.getvalue(xin)
    lb = m.colLower[xout.col]
    ub = m.colUpper[xout.col]
    if isinf(lb) || isinf(ub)
      bound_error(xout)
    end

    av = @variable(m, lowerbound=0, upperbound=Inf)
    push!(aldstates(m),
          ALDState(xin, xout, av, cc,
                   JuMP.@constraint(m, xin - xin0 <= av),
                   JuMP.@constraint(m, -xin + xin0 <= av),
                   lb,ub)
         )
end

