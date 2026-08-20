using SparseArrays
using Random

"""
    evolve(model, ψ0, time_points, observables; kwargs...)

Evolve a pure state under a `SpinModel`.

Supported methods:

- `:krylov`
- `:arnoldi`
- `:trajectories`
"""
function evolve(
    model::SpinModel,
    ψ0::SparseVector{Float64, Int64},
    time_points::Vector{Float64},
    observables::Vector{SparseMatrixCSC{Float64, Int64}};
    method::Symbol=:krylov,
    lindblad_ops=nothing,
    num_samples::Int=1,
    rng::AbstractRNG=Random.default_rng(),
)
    if method in (:krylov, :arnoldi)
        lindblad_ops === nothing ||
            throw(ArgumentError(
                "lindblad_ops is only used with method=:trajectories",
            ))

        num_samples == 1 ||
            throw(ArgumentError(
                "num_samples is only used with method=:trajectories",
            ))

        solver = KrylovArnoldiSolver.KrylovArnoldiSystem(
            model.hamiltonian,
        )

        return evolve(
            solver,
            ψ0,
            time_points,
            observables;
            method=method,
        )

    elseif method === :trajectories
        lindblad_ops === nothing &&
            throw(ArgumentError(
                "lindblad_ops is required for method=:trajectories",
            ))

        solver = TrajectorySystem(
            model,
            lindblad_ops,
        )

        return evolve(
            solver,
            ψ0,
            time_points,
            observables;
            num_samples=num_samples,
            rng=rng,
        )

    else
        throw(ArgumentError(
            "method must be :krylov, :arnoldi, or :trajectories",
        ))
    end
end


"""
    evolve(model, ρ0, time_points, observables; lindblad_ops, method=:expm)

Evolve a density matrix under Lindblad dynamics.

Supported methods:

- `:expm`
- `:ode`
"""
function evolve(
    model::SpinModel,
    ρ0::SparseMatrixCSC{Float64, Int64},
    time_points::Vector{Float64},
    observables::Vector{SparseMatrixCSC{Float64, Int64}};
    method::Symbol=:expm,
    lindblad_ops=nothing,
)
    lindblad_ops === nothing &&
        throw(ArgumentError(
            "lindblad_ops is required for Lindblad evolution",
        ))

    method in (:expm, :ode) ||
        throw(ArgumentError(
            "method must be :expm or :ode for density-matrix evolution",
        ))

    solver = LindbladSystem(
        model,
        lindblad_ops,
    )

    return evolve(
        solver,
        ρ0,
        time_points,
        observables;
        method=method,
    )
end
