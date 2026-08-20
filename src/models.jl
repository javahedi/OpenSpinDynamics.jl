
module SpinModels

using SparseArrays
using ..Operators: spin_operators
export SpinModel, update_model!

mutable struct SpinModel
    N::Int
    hamiltonian::SparseMatrixCSC{Float64, Int}
    Jxy::Float64
    Jz::Float64
    hx::Vector{Float64}
    hz::Vector{Float64}
    Jmn::AbstractMatrix{Float64}
end


"""
    SpinModel(N; Jxy=1.0, Jz=1.0, hx=zeros(N), hz=zeros(N), Jmn=zeros(N, N))

Construct a quantum spin model.

# Arguments
- `N`: Number of spins.

# Keywords
- `Jxy`: XY coupling strength.
- `Jz`: Ising coupling strength.
- `hx`: Local x-fields.
- `hz`: Local z-fields.
- `Jmn`: Coupling matrix.

# Returns
A `SpinModel` containing the Hamiltonian and model parameters.
"""
function SpinModel(
    N::Int;
    Jxy::Real=1.0,
    Jz::Real=1.0,
    hx::AbstractVector{<:Real}=zeros(N),
    hz::AbstractVector{<:Real}=zeros(N),
    Jmn::AbstractMatrix{<:Real}=zeros(N, N),
)
    N > 0 ||
        throw(ArgumentError("N must be positive"))

    length(hx) == N ||
        throw(DimensionMismatch("hx must have length $N"))

    length(hz) == N ||
        throw(DimensionMismatch("hz must have length $N"))

    size(Jmn) == (N, N) ||
        throw(DimensionMismatch("Jmn must have size ($N, $N)"))

    hx_vec = Float64.(hx)
    hz_vec = Float64.(hz)
    J_matrix = Matrix{Float64}(Jmn)

    hamiltonian = spzeros(Float64, 2^N, 2^N)
    ops = spin_operators(N)

    for i in 1:N
        hamiltonian += hx_vec[i] * ops.x[i]
        hamiltonian += hz_vec[i] * ops.z[i]
    end

    for i in 1:N-1
        for j in i+1:N
            hamiltonian +=
                2.0 * Float64(Jxy) *
                (
                    ops.plus[i] * ops.minus[j] +
                    ops.minus[i] * ops.plus[j]
                ) *
                J_matrix[i, j]

            hamiltonian +=
                Float64(Jz) *
                (ops.z[i] * ops.z[j]) *
                J_matrix[i, j]
        end
    end

    return SpinModel(
        N,
        hamiltonian,
        Float64(Jxy),
        Float64(Jz),
        hx_vec,
        hz_vec,
        J_matrix,
    )
end




"""
    update_model!(model::SpinModel, Jxy::Float64, Jz::Float64; kwargs...)

Updates the Hamiltonian of an existing `SpinModel` object with new parameters.

Arguments:
  - `model`: The `SpinModel` instance to update.
  - `Jxy`: XY coupling constant.
  - `Jz`: Z coupling constant.

Keyword arguments:
  - `hx`: New magnetic field in the x-direction.
  - `hz`: New magnetic field in the z-direction.
  - `Jmn`: New coupling matrix.

Returns:
  Nothing. Modifies the `SpinModel` in place.
"""
function update_model!(
    model::SpinModel;
    Jxy::Real=model.Jxy,
    Jz::Real=model.Jz,
    hx::AbstractVector{<:Real}=model.hx,
    hz::AbstractVector{<:Real}=model.hz,
    Jmn::AbstractMatrix{<:Real}=model.Jmn,
)
    N = model.N

    length(hx) == N ||
        throw(DimensionMismatch("hx must have length $N"))

    length(hz) == N ||
        throw(DimensionMismatch("hz must have length $N"))

    size(Jmn) == (N, N) ||
        throw(DimensionMismatch("Jmn must have size ($N, $N)"))

    hx_vec = Float64.(hx)
    hz_vec = Float64.(hz)
    J_matrix = Matrix{Float64}(Jmn)

    hamiltonian = spzeros(Float64, 2^N, 2^N)
    ops = spin_operators(N)

    for i in 1:N
        hamiltonian += hx_vec[i] * ops.x[i]
        hamiltonian += hz_vec[i] * ops.z[i]
    end

    for i in 1:N-1
        for j in i+1:N
            hamiltonian +=
                2.0 * Float64(Jxy) *
                (
                    ops.plus[i] * ops.minus[j] +
                    ops.minus[i] * ops.plus[j]
                ) *
                J_matrix[i, j]

            hamiltonian +=
                Float64(Jz) *
                (ops.z[i] * ops.z[j]) *
                J_matrix[i, j]
        end
    end

    model.hamiltonian = hamiltonian
    model.Jxy = Float64(Jxy)
    model.Jz = Float64(Jz)
    model.hx = hx_vec
    model.hz = hz_vec
    model.Jmn = J_matrix

    return model
end

end  # End of SpinModel module
