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