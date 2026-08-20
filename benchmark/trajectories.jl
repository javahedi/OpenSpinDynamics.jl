using OpenSpinDynamics
using BenchmarkTools

function setup_trajectory(N)
    coupling = NearestNeighborCoupling(0.0, N)

    model = SpinModel(
        N;
        Jxy=1.0,
        Jz=1.0,
        coupling=coupling,
    )

    ψ0 = neel_state(N)
    ops = spin_operators(N)

    γ = 0.1

    lindblad_ops = [
        sqrt(γ) * ops.minus[i]
        for i in 1:N
    ]

    times = collect(
        range(0.0, 1.0; length=21)
    )

    observables = [ops.z[1]]

    return model, ψ0, times, observables, lindblad_ops
end


function benchmark_trajectory(
    N;
    num_samples=100,
    seconds=3,
)
    model, ψ0, times, observables, lindblad_ops =
        setup_trajectory(N)

    # Warm-up
    evolve(
        model,
        ψ0,
        times,
        observables;
        method=:trajectories,
        lindblad_ops=lindblad_ops,
        num_samples=num_samples,
    )

    bench = @benchmarkable evolve(
        $model,
        $ψ0,
        $times,
        $observables;
        method=:trajectories,
        lindblad_ops=$lindblad_ops,
        num_samples=$num_samples,
    )

    return run(
        bench;
        seconds=seconds,
    )
end


function print_result(
    N,
    num_samples,
    trial,
)
    estimate = median(trial)

    println(
        N, ",",
        2^N, ",",
        num_samples, ",",
        round(
            estimate.time / 1e6;
            digits=3,
        ), ",",
        round(
            estimate.memory / 1024^2;
            digits=3,
        ), ",",
        estimate.allocs,
    )
end


println(
    "N,hilbert_dim,num_samples,time_ms,memory_MiB,allocations"
)

for N in (2, 3, 4, 5, 6, 8, 10)
    trial = benchmark_trajectory(
        N;
        num_samples=100,
    )

    print_result(
        N,
        100,
        trial,
    )
end



println(
    "N,hilbert_dim,num_samples,time_ms,memory_MiB,allocations"
)

N = 8

for num_samples in (1, 10, 100, 500)
    trial = benchmark_trajectory(
        N;
        num_samples=num_samples,
    )

    print_result(
        N,
        num_samples,
        trial,
    )
end