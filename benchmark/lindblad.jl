using OpenSpinDynamics
using BenchmarkTools
using SparseArrays

function setup_lindblad(N)
    coupling = NearestNeighborCoupling(0.0, N)

    model = SpinModel(
        N;
        Jxy=1.0,
        Jz=1.0,
        coupling=coupling,
    )

    ψ0 = neel_state(N)
    ρ0 = sparse(ψ0 * ψ0')

    ops = spin_operators(N)

    γ = 0.1
    lindblad_ops = [
        sqrt(γ) * ops.minus[i]
        for i in 1:N
    ]

    times = collect(range(0.0, 1.0; length=21))
    observables = [ops.z[1]]

    return model, ρ0, times, observables, lindblad_ops
end


function benchmark_lindblad(N; seconds=3)
    model, ρ0, times, observables, lindblad_ops =
        setup_lindblad(N)

    # Warm-up
    evolve(
        model,
        ρ0,
        times,
        observables;
        method=:expm,
        lindblad_ops=lindblad_ops,
    )

    evolve(
        model,
        ρ0,
        times,
        observables;
        method=:ode,
        lindblad_ops=lindblad_ops,
    )

    expm_bench = @benchmarkable evolve(
        $model,
        $ρ0,
        $times,
        $observables;
        method=:expm,
        lindblad_ops=$lindblad_ops,
    )

    ode_bench = @benchmarkable evolve(
        $model,
        $ρ0,
        $times,
        $observables;
        method=:ode,
        lindblad_ops=$lindblad_ops,
    )

    expm_trial = run(expm_bench; seconds=seconds)
    ode_trial = run(ode_bench; seconds=seconds)

    return expm_trial, ode_trial
end


function print_result(N, method, trial)
    estimate = median(trial)

    println(
        N, ",",
        2^N, ",",
        4^N, ",",
        method, ",",
        round(estimate.time / 1e6; digits=3), ",",
        round(estimate.memory / 1024^2; digits=3), ",",
        estimate.allocs,
    )
end


println(
    "N,hilbert_dim,liouville_dim,method,time_ms,memory_MiB,allocations"
)

# Start conservatively: Liouville space grows as 4^N.
for N in (2, 3, 4, 5)
    expm_trial, ode_trial = benchmark_lindblad(N)

    print_result(N, "expm", expm_trial)
    print_result(N, "ode", ode_trial)
end
