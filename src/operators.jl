module Operators

using SparseArrays
using LinearAlgebra

export SpinOperators, spin_operators


"""
    SpinOperators

Collection of single-site spin operators for an `N`-spin Hilbert space.

The fields contain sparse operators acting on the full Hilbert space:

- `x[i]`: Pauli ``σ_x`` operator acting on site `i`.
- `z[i]`: Pauli ``σ_z`` operator acting on site `i`.
- `plus[i]`: raising operator ``σ_+`` acting on site `i`.
- `minus[i]`: lowering operator ``σ_-`` acting on site `i`.
"""
struct SpinOperators
    x::Vector{SparseMatrixCSC{Float64, Int}}
    z::Vector{SparseMatrixCSC{Float64, Int}}
    plus::Vector{SparseMatrixCSC{Float64, Int}}
    minus::Vector{SparseMatrixCSC{Float64, Int}}
end


"""
    spin_operators(N::Int)

Construct single-site spin operators for an `N`-spin system.

Returns a [`SpinOperators`](@ref) object containing the
``σ_x``, ``σ_z``, ``σ_+``, and ``σ_-`` operators for every site,
represented as sparse matrices on the full ``2^N``-dimensional Hilbert space.

# Arguments

- `N`: Number of spins. Must be positive.

# Example

```julia
ops = spin_operators(4)

ops.z[1]
ops.plus[2]
```
"""
function spin_operators(N::Int)
    N > 0 || throw(ArgumentError("N must be positive"))

    σx = sparse([0.0 1.0; 1.0 0.0])
    σz = sparse([1.0 0.0; 0.0 -1.0])
    σplus = sparse([0.0 1.0; 0.0 0.0])
    σminus = sparse([0.0 0.0; 1.0 0.0])

    x = Vector{SparseMatrixCSC{Float64, Int}}()
    z = Vector{SparseMatrixCSC{Float64, Int}}()
    plus = Vector{SparseMatrixCSC{Float64, Int}}()
    minus = Vector{SparseMatrixCSC{Float64, Int}}()

    for i in 0:N-1
        leftdim = 2^i
        rightdim = 2^(N - i - 1)

        I_left = sparse(I, leftdim, leftdim)
        I_right = sparse(I, rightdim, rightdim)

        push!(x, kron(kron(I_left, σx), I_right))
        push!(z, kron(kron(I_left, σz), I_right))
        push!(plus, kron(kron(I_left, σplus), I_right))
        push!(minus, kron(kron(I_left, σminus), I_right))
    end

    return SpinOperators(
        x,
        z,
        plus,
        minus,
    )
end

end # module Operators
