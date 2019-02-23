# Analytic solution exploiting symmetry
cost(x::Float64) = abs(x)
next_step(x::Float64) = (x >=0 ? x-1 : x+1)

function Q2_bar(x2::T, nmax::Int, noise::AbstractArray{T,1}, beta=0.5::Float64, t=1::Int) where T
    if t >= nmax
        return 0.0
    end

    nsamples = length(noise)
    ans = 0.0
    for xi in noise
        ans += solve_bin_sym(x2+xi, nmax, noise, beta, t+1)
    end
    return ans/nsamples
end

function solve_bin_sym(x::T, nmax::Int, noise::AbstractArray{T,1}, beta=0.5::Float64, t=1::Int) where T
    if t > nmax
        return 0.0
    end

    x2 = next_step(x)
    return cost(x2) + beta*Q2_bar(x2, nmax, noise, beta, t)
end

function graph_fcfs(m, t, ts, QTrue; filename=nothing)
    qe = QTrue
    qt = ASDDiP.Qtilde(m,t,1,ts)
    qf = ASDDiP.Qfrak(m,t,1,ts)

    fig, (ax1,ax2) = PyPlot.subplots(ncols=2, figsize=(12,4))
    PyPlot.suptitle("Stage $t")

    ax1[:plot](ts, qe, label=L"$\overline{Q}$: exact")
    ax1[:plot](ts, qt, label=L"$\tilde{Q}$: mean of next stage approximations")
    ax1[:plot](ts, qf, label=L"$\mathfrak{Q}$: current stage approximation")
    ax1[:legend]()
    ax1[:set_title]("Future cost functions")
    ax1[:grid]()

    ax2[:plot](ts, qt - qf, label=L"$\tilde{Q} - \mathfrak{Q}$")
    ax2[:plot](ts, qe - qf, label=L"$\overline{Q} - \mathfrak{Q}$")
    ax2[:plot](ts, qe - qt, label=L"$\overline{Q} - \tilde{Q}$")
    ax2[:set_title]("Differences")
    ax2[:legend]()
    ax2[:grid]()

    if filename != nothing
        PyPlot.savefig(filename)
    end
end

function compare_qfrak(models, t, ts, QTrue; filename=nothing)
    qe = QTrue
    qf = [ASDDiP.Qfrak(m,t,1,ts) for m in models]

    fig, (ax1,ax2) = PyPlot.subplots(ncols=2, figsize=(12,4))
    PyPlot.suptitle("Stage $t")

    for (i,qfi) in enumerate(qf)
        ax1[:plot](ts, qfi, label="$i")
    end
    ax1[:plot](ts, qe, label=L"\overline{Q}")
    ax1[:legend]()
    ax1[:set_title](L"$\mathfrak{Q}$ comparison")
    ax1[:grid]()

    for (i, qfi) in enumerate(qf)
        ax2[:plot](ts, qe - qfi, label="$i")
    end
    ax2[:set_title](L"Differences $\overline{Q} - \mathfrak{Q}$")
    ax2[:legend]()
    ax2[:grid]()

    if filename != nothing
        PyPlot.savefig(filename)
    end
end
