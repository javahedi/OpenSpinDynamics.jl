module Operators

using SparseArrays
using LinearAlgebra

export SpinOperators, spin_operators

struct SpinOperators
    x::Vector{SparseMatrixCSC{Float64, Int}}
    z::Vector{SparseMatrixCSC{Float64, Int}}
    plus::Vector{SparseMatrixCSC{Float64, Int}}
    minus::Vector{SparseMatrixCSC{Float64, Int}}
end

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

    return SpinOperators(x, z, plus, minus)
end

end # module Operators