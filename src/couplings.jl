module Coupling

using Random

export LongRangeCouplingDisorder,
       LongRangeCouplingClean,
       NearestNeighborCoupling

abstract type AbstractCoupling end


"""
    LongRangeCouplingDisorder

Random long-range power-law coupling between `N` spins distributed over
a lattice of size `L`.

Couplings decay with distance as

```math
J_{ij} = \\frac{1}{|r_i-r_j|^\\alpha}.
```

The generated coupling matrix is symmetric with zero diagonal.

# Fields

- `matrix`: Coupling matrix.
- `α`: Power-law exponent.
- `L`: Lattice size.
- `N`: Number of spins.
- `reordered`: Whether sampled sites are reordered.
"""
struct LongRangeCouplingDisorder <: AbstractCoupling
    matrix::Matrix{Float64}
    α::Float64
    L::Int
    N::Int
    reordered::Bool
end


"""
    LongRangeCouplingClean

Clean long-range power-law coupling between `N` spins.

The interaction between sites `i` and `j` is

```math
J_{ij} = \\frac{1}{|i-j|^\\alpha},
```

with zero diagonal.
"""
struct LongRangeCouplingClean <: AbstractCoupling
    matrix::Matrix{Float64}
    α::Float64
    N::Int
end


"""
    NearestNeighborCoupling

Nearest-neighbor coupling with alternating bond strengths.

For neighboring sites, the coupling is

```math
J_{n,n+1} = 1 + (-1)^n \\delta,
```

while all non-nearest-neighbor couplings vanish.
"""
struct NearestNeighborCoupling <: AbstractCoupling
    matrix::Matrix{Float64}
    δ::Float64
    N::Int
end


"""
    LongRangeCouplingDisorder(
        rng::AbstractRNG,
        α::Float64,
        L::Int,
        N::Int;
        reordered::Bool=true,
    )

Construct a reproducible random long-range power-law coupling using
the supplied random-number generator `rng`.

`N` spin positions are sampled from a lattice containing `L` sites,
and the interaction between two selected spins decays as ``1/r^α``.

# Arguments

- `rng`: Random-number generator.
- `α`: Power-law exponent.
- `L`: Number of available lattice sites.
- `N`: Number of spins.

# Keywords

- `reordered=true`: Reorder sampled positions using the internal
  pairing procedure.
"""
function LongRangeCouplingDisorder(
    rng::AbstractRNG,
    α::Float64,
    L::Int,
    N::Int;
    reordered::Bool=true,
)
    L > 0 || throw(ArgumentError("L must be positive."))
    N > 0 || throw(ArgumentError("N must be positive."))
    N <= L || throw(ArgumentError("N must not exceed L."))

    r_index = sort(randperm(rng, L)[1:N])

    if reordered
        _, r_index = reorder_indices(r_index)
    end

    Jmn = zeros(Float64, N, N)

    for i in 1:N-1
        for j in i+1:N
            dij = abs(r_index[i] - r_index[j])
            Jmn[i, j] = 1.0 / dij^α
        end
    end

    Jmn .+= Jmn'

    return LongRangeCouplingDisorder(Jmn, α, L, N, reordered)
end


"""
    LongRangeCouplingDisorder(α, L, N; reordered=true)

Construct a random long-range power-law coupling.

`N` spin positions are sampled from a lattice containing `L` sites.
The interaction between two selected spins decays as ``1/r^α``.

# Arguments

- `α`: Power-law exponent.
- `L`: Number of available lattice sites.
- `N`: Number of spins.

# Keywords

- `reordered=true`: Reorder sampled positions using the internal
  pairing procedure.

# Example

```julia
coupling = LongRangeCouplingDisorder(2.0, 100, 10)

coupling.matrix
```
"""
function LongRangeCouplingDisorder(
    α::Float64,
    L::Int,
    N::Int;
    reordered::Bool=true,
)
    return LongRangeCouplingDisorder(
        Random.default_rng(),
        α,
        L,
        N;
        reordered=reordered,
    )
end


"""
    LongRangeCouplingClean(α, N)

Construct a clean long-range power-law coupling for `N` spins.

The interaction strength decays with distance as ``1/r^α``.

# Arguments

- `α`: Power-law exponent.
- `N`: Number of spins.

# Example

```julia
coupling = LongRangeCouplingClean(2.0, 6)

coupling.matrix
```
"""
function LongRangeCouplingClean(α::Float64, N::Int)
    N > 0 || throw(ArgumentError("N must be positive."))

    Jmn = zeros(Float64, N, N)

    for i in 1:N-1
        for j in i+1:N
            dij = abs(i - j)
            Jmn[i, j] = 1.0 / dij^α
        end
    end

    Jmn .+= Jmn'

    return LongRangeCouplingClean(Jmn, α, N)
end


"""
    NearestNeighborCoupling(δ, N)

Construct an alternating nearest-neighbor coupling for `N` spins.

The bond strengths alternate between ``1-δ`` and ``1+δ``.

# Arguments

- `δ`: Alternation strength.
- `N`: Number of spins.

# Example

```julia
coupling = NearestNeighborCoupling(0.2, 6)

coupling.matrix
```
"""
function NearestNeighborCoupling(δ::Float64, N::Int)
    N > 0 || throw(ArgumentError("N must be positive."))

    Jmn = zeros(Float64, N, N)

    for n in 1:N-1
        Jmn[n, n + 1] = 1.0 + (-1.0)^n * δ
    end

    Jmn .+= Jmn'

    return NearestNeighborCoupling(Jmn, δ, N)
end


function reorder_indices(
    r_index,
    i_index=collect(1:length(r_index)),
)
    r_remaining = collect(r_index)
    i_remaining = collect(i_index)

    reordered_r_index = Int[]
    reordered_i_index = Int[]

    while length(r_remaining) >= 2
        min_distance = Inf
        min_indices = nothing

        for i in 1:length(r_remaining)-1
            for j in i+1:length(r_remaining)
                distance = abs(r_remaining[i] - r_remaining[j])

                if distance < min_distance
                    min_distance = distance
                    min_indices = (i, j)
                end
            end
        end

        i, j = min_indices

        push!(reordered_r_index, r_remaining[i], r_remaining[j])
        push!(reordered_i_index, i_remaining[i], i_remaining[j])

        keep = setdiff(eachindex(r_remaining), (i, j))

        r_remaining = r_remaining[keep]
        i_remaining = i_remaining[keep]
    end

    # Preserve an unpaired site for odd N.
    if length(r_remaining) == 1
        push!(reordered_r_index, r_remaining[1])
        push!(reordered_i_index, i_remaining[1])
    end

    return reordered_i_index, reordered_r_index
end

end # module Coupling
