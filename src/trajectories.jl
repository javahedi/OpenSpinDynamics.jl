module TrajectorySolver

using LinearAlgebra
using SparseArrays
using Statistics

import ..evolve
using ..SpinModels: SpinModel

export TrajectorySystem

struct TrajectorySystem
    hamiltonian::SparseMatrixCSC{ComplexF64, Int64}
    lindblad_ops::Vector{SparseMatrixCSC{ComplexF64, Int64}}
    effective_hamiltonian::SparseMatrixCSC{ComplexF64, Int64}
end

function TrajectorySystem(
    model::SpinModel,
    lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}},
)
    return TrajectorySystem(
        model.hamiltonian,
        lindblad_ops,
    )
end

function TrajectorySystem(
    hamiltonian::SparseMatrixCSC{Float64, Int64},
    lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}},
)
    dim = size(hamiltonian, 1)

    size(hamiltonian, 2) == dim ||
        throw(ArgumentError("Hamiltonian must be square"))

    all(size(L) == (dim, dim) for L in lindblad_ops) ||
        throw(ArgumentError(
            "Lindblad operators must have the same dimensions as the Hamiltonian",
        ))

    H = complex.(hamiltonian)
    Ls = [complex.(L) for L in lindblad_ops]

    Heff = copy(H)

    for L in Ls
        Heff -= (1im / 2) * (L' * L)
    end

    return TrajectorySystem(H, Ls, Heff)
end

# function _stochastic_evolution(
#     solver::TrajectorySystem,
#     ψ0::Vector{ComplexF64},
#     time_points::Vector{Float64},
#     observables::Vector{SparseMatrixCSC{ComplexF64, Int64}},
#     num_samples::Int,
# )
#     m = length(time_points)

#     LdagL = [L' * L for L in solver.lindblad_ops]

#     values = zeros(
#         Float64,
#         m,
#         length(observables),
#         num_samples,
#     )

#     for sample in 1:num_samples
#         ψ = copy(ψ0)

#         for (j, observable) in enumerate(observables)
#             values[1, j, sample] =
#                 real(dot(ψ, observable * ψ))
#         end

#         for i in 2:m
#             dt = time_points[i] - time_points[i - 1]

#             dps = [
#                 real(dt * dot(ψ, Ldag * ψ))
#                 for Ldag in LdagL
#             ]

#             dP = sum(dps)

#             if rand() > dP
#                 ψ = (I - 1im * solver.effective_hamiltonian * dt) * ψ
#             else
#                 probabilities = cumsum(dps) / dP
#                 k = searchsortedfirst(probabilities, rand())
#                 ψ = solver.lindblad_ops[k] * ψ
#             end

#             ψ ./= norm(ψ)

#             for (j, observable) in enumerate(observables)
#                 values[i, j, sample] =
#                     real(dot(ψ, observable * ψ))
#             end
#         end
#     end

#     mean_values = dropdims(mean(values, dims=3), dims=3)

#     std_values = if num_samples == 1
#         zeros(Float64, m, length(observables))
#     else
#         dropdims(std(values, dims=3), dims=3)
#     end

#     return mean_values, std_values
# end

function _stochastic_evolution(
    solver::TrajectorySystem,
    ψ0::Vector{ComplexF64},
    time_points::Vector{Float64},
    observables::Vector{SparseMatrixCSC{ComplexF64, Int64}},
    num_samples::Int,
)
    m = length(time_points)
    nobs = length(observables)
    njumps = length(solver.lindblad_ops)

    # Construct L†L once per evolve call.
    LdagL = [
        L' * L
        for L in solver.lindblad_ops
    ]

    values = zeros(
        Float64,
        m,
        nobs,
        num_samples,
    )

    # Reusable working storage.
    ψ = copy(ψ0)
    tmp = similar(ψ0)
    dps = zeros(Float64, njumps)

    for sample in 1:num_samples
        copyto!(ψ, ψ0)

        # Initial observables.
        for (j, observable) in enumerate(observables)
            mul!(tmp, observable, ψ)
            values[1, j, sample] =
                real(dot(ψ, tmp))
        end

        for i in 2:m
            dt = time_points[i] - time_points[i - 1]

            # Jump probabilities without allocating temporary vectors.
            dP = 0.0

            for k in eachindex(LdagL)
                mul!(tmp, LdagL[k], ψ)

                dp = real(
                    dt * dot(ψ, tmp)
                )

                dps[k] = dp
                dP += dp
            end

            if rand() > dP
                # ψ ← (I - i Heff dt) ψ
                #
                # Compute Heff * ψ into preallocated storage.
                mul!(
                    tmp,
                    solver.effective_hamiltonian,
                    ψ,
                )

                @. ψ =
                    ψ -
                    1im * dt * tmp
            else
                # Select a jump without cumsum(dps) / dP.
                threshold = rand() * dP
                cumulative = 0.0
                selected = lastindex(dps)

                for k in eachindex(dps)
                    cumulative += dps[k]

                    if threshold <= cumulative
                        selected = k
                        break
                    end
                end

                # Apply jump into reusable storage.
                mul!(
                    tmp,
                    solver.lindblad_ops[selected],
                    ψ,
                )

                copyto!(ψ, tmp)
            end

            ψ ./= norm(ψ)

            # Observables without observable * ψ allocations.
            for (j, observable) in enumerate(observables)
                mul!(tmp, observable, ψ)

                values[i, j, sample] =
                    real(dot(ψ, tmp))
            end
        end
    end

    mean_values =
        dropdims(
            mean(values, dims=3),
            dims=3,
        )

    std_values = if num_samples == 1
        zeros(
            Float64,
            m,
            nobs,
        )
    else
        dropdims(
            std(values, dims=3),
            dims=3,
        )
    end

    return mean_values, std_values
end


function evolve(
    solver::TrajectorySystem,
    ψ0::SparseVector{Float64, Int64},
    time_points::Vector{Float64},
    observables::Vector{SparseMatrixCSC{Float64, Int64}};
    num_samples::Int=1,
)
    isempty(time_points) &&
        throw(ArgumentError("time_points must not be empty"))

    issorted(time_points) ||
        throw(ArgumentError("time_points must be sorted"))

    num_samples > 0 ||
        throw(ArgumentError("num_samples must be positive"))

    length(ψ0) == size(solver.hamiltonian, 1) ||
        throw(ArgumentError(
            "Initial state must have compatible dimensions with the Hamiltonian",
        ))

    all(size(obs) == size(solver.hamiltonian) for obs in observables) ||
        throw(ArgumentError(
            "Observables must have the same dimensions as the Hamiltonian",
        ))

    ψ = Vector{ComplexF64}(ψ0)
    obs = [complex.(observable) for observable in observables]

    return _stochastic_evolution(
        solver,
        ψ,
        time_points,
        obs,
        num_samples,
    )
end


end
