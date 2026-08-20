module LindbladSolver

    using SparseArrays
    using LinearAlgebra
    using OrdinaryDiffEqTsit5: Tsit5
    using SciMLBase: ODEProblem, solve
    using Expokit
    import ..evolve
    using ..SpinModels: SpinModel

    export LindbladSystem

    struct LindbladSystem
        superop::SparseMatrixCSC{ComplexF64, Int64}
        hamiltonian::SparseMatrixCSC{Float64, Int64}
        lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}}
    end


    function LindbladSystem(
        model::SpinModel,
        lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}},
    )
        return LindbladSystem(
            model.hamiltonian,
            lindblad_ops,
        )
    end

    # Constructor for Lindblad system
    function LindbladSystem(
        hamiltonian::SparseMatrixCSC{Float64, Int64},
        lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}},
    )
        dim = size(hamiltonian, 1)

        size(hamiltonian, 2) == dim ||
            throw(ArgumentError("Hamiltonian must be square"))

        all(size(op) == (dim, dim) for op in lindblad_ops) ||
            throw(ArgumentError(
                "Lindblad operators must have the same dimensions as the Hamiltonian",
            ))

        identity = sparse(I, dim, dim)

        superH = -1im * (
            kron(identity, hamiltonian) -
            kron(transpose(hamiltonian), identity)
        )

        superL = spzeros(ComplexF64, dim^2, dim^2)

        for op in lindblad_ops
            opdag_op = op' * op

            superL += kron(conj(op), op)

            superL -= 0.5 * kron(identity, opdag_op)

            superL -= 0.5 * kron(
                transpose(op) * conj(op),
                identity,
            )
        end

        superop = superH + superL

        return LindbladSystem(
            superop,
            hamiltonian,
            lindblad_ops,
        )
    end


    function _validate_dimensions(
        solver::LindbladSystem,
        ρ0,
        observables,
    )
        dim = size(solver.hamiltonian, 1)

        size(ρ0) == (dim, dim) ||
            throw(DimensionMismatch(
                "Initial density matrix must have size ($dim, $dim)",
            ))

        all(size(obs) == (dim, dim) for obs in observables) ||
            throw(DimensionMismatch(
                "Observables must have the same dimensions as the Hamiltonian",
            ))

        return nothing
    end

    function _rho_dot!(dρ, state, superop, time)
        dρ .= superop * state
    end
   
   function _solver_ode(
                solver::LindbladSystem,
                ρ0::SparseMatrixCSC{ComplexF64, Int},
                time_points::Vector{Float64},
                observables::Vector{SparseMatrixCSC{Float64, Int}},
                solver_algorithm=Tsit5(),
            )

        # Validate dimensions
        #@assert size(ρ0[:], 1) == size(solver.superop, 1) "Initial state must have compatible dimensions with the superoperator."
        #@assert all(size(obs) == size(ρ0) for obs in observables) "Observables must have the same dimensions as the initial state."
        _validate_dimensions(solver, ρ0, observables)
        
        time_span = (time_points[1],time_points[end])

        #ρ0_vector =  Vector(ρ0[:]) # Vector{Float64}
        # ρ0[:] --> SparseVector{Float64, Int64}
        problem = ODEProblem((dρ, ρ, params, time) -> _rho_dot!(dρ, ρ, solver.superop, time), ρ0[:], time_span)
        solution = solve(problem, solver_algorithm, saveat=time_points, verbose=false)

        outputs = Matrix{Float64}(undef, length(time_points), length(observables))
        for (j, observable) in enumerate(observables)
            for i in 1:length(time_points)
                ρt = reshape(solution.u[i], size(ρ0) )
                outputs[i, j] = real(tr(ρt * observable))
            end
        end
    
        return outputs
    end

    function _propagate_ρt(solver::LindbladSystem, 
                        ρ0::SparseMatrixCSC{ComplexF64, Int}, 
                        time::Float64)
        #ρ0_vector =  Vector(ρ0[:]) # Vector{Float64}
        # ρ0[:] --> SparseVector{Float64, Int64}
        ρt = expmv(time, solver.superop, Vector(ρ0[:]))
        return reshape(ρt, size(ρ0))
    end

   
   function _solver_expm(
        solver::LindbladSystem,
        ρ0::SparseMatrixCSC{ComplexF64, Int},
        time_points::Vector{Float64},
        observables::Vector{SparseMatrixCSC{Float64, Int}},
    )
        #@assert size(ρ0[:], 1) == size(solver.superop, 1)
        #@assert all(size(obs) == size(ρ0) for obs in observables)
        _validate_dimensions(solver, ρ0, observables)

        shifted_times = time_points .- time_points[1]

        outputs = Matrix{Float64}(
            undef,
            length(time_points),
            length(observables),
        )

        for (i, time) in enumerate(shifted_times)
            ρt = _propagate_ρt(solver, ρ0, time)

            for (j, observable) in enumerate(observables)
                outputs[i, j] = real(tr(ρt * observable))
            end
        end

        return outputs
    end
    
    
    function evolve(
        solver::LindbladSystem,
        ρ0::SparseMatrixCSC{Float64, Int},
        time_points::Vector{Float64},
        observables::Vector{SparseMatrixCSC{Float64, Int}};
        method::Symbol=:expm,
    )
        isempty(time_points) &&
            throw(ArgumentError("time_points must not be empty"))

        issorted(time_points) ||
            throw(ArgumentError("time_points must be sorted in ascending order"))

        method in (:expm, :ode) ||
            throw(ArgumentError("method must be :expm or :ode"))

        if method === :expm
            return _solver_expm(
                solver,
                complex.(ρ0),
                time_points,
                observables,
            )
        else
            return _solver_ode(
                solver,
                complex.(ρ0),
                time_points,
                observables,
            )
        end
    end

end  # End of LiouvillianSolver module