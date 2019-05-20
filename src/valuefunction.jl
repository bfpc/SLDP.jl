# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function add_aldcut!(sp::JuMP.Model, cut::ALDCut)
  push!(aldcuts(sp), cut)
  idx = length(aldcuts(sp))

  aff_expr = cut.v
  for (st,xi,li) in zip(aldstates(sp), cut.xhat, cut.l)
    pos_part = @variable(sp, lowerbound=0, basename="theta_pos_"*string(idx))
    neg_part = @variable(sp, lowerbound=0, basename="theta_neg_"*string(idx))
    bin_choose = @variable(sp, category=:Bin, basename="theta_bin_"*string(idx))

    @constraint(sp, pos_part <= bin_choose * (st.ub - xi))
    @constraint(sp, neg_part <= (1-bin_choose) * (xi - st.lb))
    @constraint(sp, pos_part - neg_part == st.xout - xi)
    aff_expr += li*(st.xout - xi) - cut.rho*sum(pos_part+neg_part)
  end

  _theta = SDDP.valueoracle(sp).theta
  @constraint(sp, _theta >= aff_expr)
end


function add_aldcut!(m::SDDP.SDDPModel, sp::JuMP.Model, v, l, rho)
  t = SDDP.ext(sp).stage
  xhat = SDDP.getstage(m, t).state

  add_aldcut!(sp, ALDCut(copy(xhat), v, l, rho))
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
        SDDP.modify_probability(ex.riskmeasure,
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

"""
    load_aldcuts!(m::SDDPModel, filename::String)

Load ALD cuts from the file created using save_aldcuts(m, filename).

### Example

    m = SDDPModel() do ... end
    status = solve(m; cut_output_file="path/to/m.cuts")`
    save_aldcuts!(m; filename="path/to/m_ald.cuts")
    m2 = SDDPModel() do ... end
    loadcuts!(m2, "path/to/m.cuts")
    load_aldcuts!(m; filename="path/to/m_ald.cuts")

"""
function load_aldcuts!(m::SDDP.SDDPModel{SDDP.DefaultValueFunction{C}}, filename::String) where C
    open(filename, "r") do file
        while true
            line      = readline(file)
            line == nothing || line == "" && break
            items     = split(line, ",")
            stage     = parse(Int, items[1])
            ms        = parse(Int, items[2])
            dim       = div(length(items) - 4, 2)
            xhat      = [parse(Float64, i) for i in items[3:2+dim]]
            intercept = parse(Float64, items[3+dim])
            coefficients = [parse(Float64, i) for i in items[4+dim:end-1]]
            rho       = parse(Float64, items[end])
            cut = ALDCut(xhat, intercept, coefficients, rho)
            sp = SDDP.getsubproblem(m, stage, ms)
            add_aldcut!(sp, cut)
        end
    end
end

function load_aldcuts!(m::SDDP.SDDPModel{SDDP.DefaultValueFunction{C}}, cutlist::Vector{Vector{Vector{ALDCut}}}; lastn=nothing::Union{Void,Int64}) where C
    for (stage,stagecuts) in enumerate(cutlist)
        for (ms,spcuts) in enumerate(stagecuts)
            sp = SDDP.getsubproblem(m, stage, ms)
            selectedcuts = ALDCut[]
            if lastn == nothing || length(spcuts) < lastn
                selectedcuts = spcuts
            else
                selectedcuts = spcuts[end-lastn:end]
            end
            for cut in selectedcuts
                add_aldcut!(sp, cut)
            end
        end
    end
end

aldcuts(stage::SDDP.Stage) = aldcuts.(stage.subproblems)
aldcuts(m::SDDP.SDDPModel) = aldcuts.(m.stages)

function write_aldcut!(io::IO, stage::Int, markovstate::Int, cut::ALDCut)
    write(io, string(stage), ",", string(markovstate))
    for xi in cut.xhat
        write(io, ",", string(xi))
    end
    write(io, ",", string(cut.v))
    for li in cut.l
        write(io, ",", string(li))
    end
    write(io, ",", string(cut.rho), "\n")
end

function save_aldcuts!(m::SDDP.SDDPModel{SDDP.DefaultValueFunction{C}}, filename::String) where C
    open(filename, "w") do file
        for (t, stage) in enumerate(SDDP.stages(m))
            for (j, sp) in enumerate(SDDP.subproblems(stage))
                for cut in aldcuts(sp)
                    write_aldcut!(file, t, j, cut)
                end
            end
        end
    end
end
