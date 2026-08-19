module Solvers

using OpenSpinDynamics
using SparseArrays

export solver_function

function solver_function(
    model,
    Cop,
    ψ0,
    observables,
    TIMEPOINTS,
    solver_type::Symbol,
    method::Symbol,
    n_sample::Int,
)
    solver_type in (:KrylovArnoldi, :Lindblad, :swf) ||
        throw(ArgumentError(
            "solver_type must be :KrylovArnoldi, :Lindblad, or :swf",
        ))

    if solver_type === :KrylovArnoldi
        return _solve_with_KrylovArnoldi(
            model,
            ψ0,
            observables,
            TIMEPOINTS;
            method=method,
        )

    elseif solver_type === :Lindblad
        ρ0 = sparse(ψ0 * ψ0')

        return _solve_with_Lindblad(
            model,
            Cop,
            ρ0,
            observables,
            TIMEPOINTS;
            method=method,
        )

    else
        return _solve_with_StochasticWavefunction(
            model,
            Cop,
            ψ0,
            observables,
            TIMEPOINTS;
            num_samples=n_sample,
        )
    end
end

function _solve_with_KrylovArnoldi(
    model,
    ψ0,
    observables,
    TIMEPOINTS;
    method::Symbol=:krylov,
)
    solver = KrylovArnoldiSystem(model.hamiltonian)

    return evolve(
        solver,
        ψ0,
        TIMEPOINTS,
        observables;
        method=method,
    )
end

function _solve_with_Lindblad(
    model,
    Cop,
    ρ0,
    observables,
    TIMEPOINTS;
    method::Symbol=:expm,
)
    solver = LindbladSystem(model.hamiltonian, Cop)

    return evolve(
        solver,
        ρ0,
        TIMEPOINTS,
        observables;
        method=method,
    )
end

function _solve_with_StochasticWavefunction(
    model,
    Cop,
    ψ0,
    observables,
    TIMEPOINTS;
    num_samples::Int=100,
)
    solver = StochasticWavefunctionSystem(model.hamiltonian, Cop)

    mean_values, _ = evolve(
        solver,
        ψ0,
        TIMEPOINTS,
        observables;
        num_samples=num_samples,
    )

    return mean_values
end

end