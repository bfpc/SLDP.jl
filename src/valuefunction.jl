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


