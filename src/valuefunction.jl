function add_aldcut!(m::SDDP.SDDPModel, sp::JuMP.Model, v, l, rho)
  _theta = SDDP.valueoracle(sp).theta
  t = SDDP.ext(sp).stage
  xhat = SDDP.getstage(m, t).state
  niter = length(m.log)

  push!(aldcuts(sp), ALDCut(copy(xhat), v, l, rho))

  aff_expr = v
  for (st,xi,li) in zip(aldstates(sp), xhat, l)
    pos_part = @variable(sp, lowerbound=0, basename="theta_pos_"*string(niter))
    neg_part = @variable(sp, lowerbound=0, basename="theta_neg_"*string(niter))
    bin_choose = @variable(sp, category=:Bin, basename="theta_bin_"*string(niter))

    @constraint(sp, pos_part <= bin_choose * (st.ub - xi))
    @constraint(sp, neg_part <= (1-bin_choose) * (xi - st.lb))
    @constraint(sp, pos_part - neg_part == st.xout - xi)
    aff_expr += li*(st.xout - xi) - rho*sum(pos_part+neg_part)
  end
  @constraint(sp, _theta >= aff_expr)
end

function modify_ald_valuefunction!(m::SDDP.SDDPModel{V}, settings::SDDP.Settings, sp::JuMP.Model) where V<:SDDP.DefaultValueFunction
    ex = SDDP.ext(sp)
    vf = SDDP.valueoracle(sp)

    t = ex.stage
    vs = zeros(0)
    ls = []
    rhos = zeros(0)
    for sp_next in SDDP.subproblems(m,t+1)
      append!(vs, sp_next.ext[:ALD].vstore)
      append!(ls, sp_next.ext[:ALD].lstore)
      append!(rhos, sp_next.ext[:ALD].rhostore)
    end

    I = 1:length(vs)
    # TODO: improve this to reduce the extra memory usage
    current_transition = copy(m.storage.probability.data[I])
    Pij = SDDP.getstage(m, t+1).transitionprobabilities
    for i in I
        current_transition[i] *= Pij[ex.markovstate, m.storage.markov[i]]
    end
    modified_probability = zeros(current_transition)
    @timeit SDDP.TIMER "risk measure" begin
        SDDP.modifyprobability!(ex.riskmeasure,
            modified_probability,
            current_transition,
            vs,
            m,
            sp
        )
    end

    vbar = dot(modified_probability, vs)
    lbar = sum(mpi*li for (mpi,li) in zip(modified_probability, ls))
    rhobar = dot(modified_probability, rhos)

    if rhobar > 0
      add_aldcut!(m, sp, vbar, lbar, rhobar)
    end
end


