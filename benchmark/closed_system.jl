using OpenSpinDynamics
using BenchmarkTools

function setup_closed(N)
    coupling = NearestNeighborCoupling(0.0, N)

    model = SpinModel(
        N;
        Jxy=1.0,
        Jz=1.0,
        coupling=coupling,
    )

    ψ0 = neel_state(N)
    ops = spin_operators(N)

    times = collect(
        range(0.0, 1.0; length=21)
    )

    observables = [ops.z[1]]

    return model, ψ0, times, observables
end


function benchmark_closed(N; seconds=3)
    model, ψ0, times, observables = setup_closed(N)

    # Warm-up so compilation is excluded.
    evolve(
        model,
        ψ0,
        times,
        observables;
        method=:krylov,
    )

    evolve(
        model,
        ψ0,
        times,
        observables;
        method=:arnoldi,
    )

    krylov_bench = @benchmarkable evolve(
        $model,
        $ψ0,
        $times,
        $observables;
        method=:krylov,
    )

    arnoldi_bench = @benchmarkable evolve(
        $model,
        $ψ0,
        $times,
        $observables;
        method=:arnoldi,
    )

    krylov = run(
        krylov_bench;
        seconds=seconds,
    )

    arnoldi = run(
        arnoldi_bench;
        seconds=seconds,
    )

    return krylov, arnoldi
end


function print_result(N, method, trial)
    estimate = median(trial)

    println(
        N, ",",
        method, ",",
        round(estimate.time / 1e6; digits=3), ",",
        round(estimate.memory / 1024^2; digits=3), ",",
        estimate.allocs,
    )
end


println("N,method,time_ms,memory_MiB,allocations")

for N in (4, 6, 8, 10)
    krylov, arnoldi = benchmark_closed(N)

    print_result(
        N,
        "krylov",
        krylov,
    )

    print_result(
        N,
        "arnoldi",
        arnoldi,
    )
end