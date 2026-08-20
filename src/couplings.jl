module Coupling

    using Random

    export LongRangeCouplingDisorder,
        LongRangeCouplingClean,
        NearestNeighborCoupling
        

    abstract type AbstractCoupling end

    # Define concrete types
    struct LongRangeCouplingDisorder <: AbstractCoupling
        matrix::Matrix{Float64}  # Coupling matrix
        α::Float64               # Power-law exponent
        L::Int                   # Lattice size
        N::Int                   # Number of spins
        reordered::Bool          # Whether indices are reordered
    end

    struct LongRangeCouplingClean <: AbstractCoupling
        matrix::Matrix{Float64}  # Coupling matrix
        α::Float64               # Power-law exponent
        N::Int                   # Number of spins
    end

    struct NearestNeighborCoupling <: AbstractCoupling
        matrix::Matrix{Float64}  # Coupling matrix
        δ::Float64               # Alternating strength
        N::Int                   # Number of spins
    end

    # # Constructor for long-range random coupling
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






    # Constructor for long-range coupling
    function LongRangeCouplingClean(α::Float64, N::Int)
       
        
        # Initialize the coupling matrix Jmn
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


    # Constructor for nearest-neighbor coupling
    function NearestNeighborCoupling(δ::Float64, N::Int)
        

        # Initialize the coupling matrix Jmn
        Jmn = zeros(Float64, N, N)
        for n in 1:N-1
            Jmn[n, n+1] = (1.0 + (-1.0)^n * δ)
        end

        Jmn .+= Jmn'
        return NearestNeighborCoupling(Jmn, δ, N)
    end



    function reorder_indices(r_index, i_index=collect(1:length(r_index)))
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

    # Common interface functions
    get_matrix(coupling::AbstractCoupling) = coupling.matrix
    get_N(coupling::AbstractCoupling) = coupling.N

end  # End of Coupling module