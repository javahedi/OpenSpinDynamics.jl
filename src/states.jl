module QuantumState

using SparseArrays

export neel_state, polarized_state

function _basis_states(direction::Symbol)
    direction in (:x, :z) ||
        throw(ArgumentError("direction must be :x or :z"))

    up_z = sparse([1.0, 0.0])
    dn_z = sparse([0.0, 1.0])

    if direction === :z
        return up_z, dn_z
    end

    factor = inv(sqrt(2.0))

    up_x = factor .* (up_z + dn_z)
    dn_x = factor .* (up_z - dn_z)

    return up_x, dn_x
end


"""
    neel_state(N; direction=:z)

Construct a normalized Néel product state for `N` spins.

The first spin points along the positive `direction`, followed by
alternating spins. `direction` may be `:z` or `:x`.
"""
function neel_state(
    N::Int;
    direction::Symbol=:z,
)
    N > 0 ||
        throw(ArgumentError("N must be positive"))

    up, dn = _basis_states(direction)

    ψ = up

    for i in 2:N
        ψ = kron(ψ, iseven(i) ? dn : up)
    end

    return sparse(ψ)
end


"""
    polarized_state(N; direction=:z)

Construct a normalized fully polarized product state for `N` spins.

`direction` may be `:z` or `:x`.
"""
function polarized_state(
    N::Int;
    direction::Symbol=:z,
)
    N > 0 ||
        throw(ArgumentError("N must be positive"))

    up, _ = _basis_states(direction)

    ψ = up

    for _ in 2:N
        ψ = kron(ψ, up)
    end

    return sparse(ψ)
end

end

