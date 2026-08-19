
using LinearAlgebra
using SparseArrays

function arnoldi(H::SparseMatrixCSC{Float64, Int64}, 
                ψ, t::Float64; 
                max_dim::Int=200, tol::Float64=1e-12)

    m, n = size(H)
    if m != n
        throw(ArgumentError("H must be a square matrix"))
    end
    if length(ψ) != m
        throw(ArgumentError("The length of ψ must match the dimensions of H"))
    end

    # Initialize variables
    h = zeros(Complex{Float64}, max_dim + 1, max_dim)
    Q = zeros(Complex{Float64}, m, max_dim + 1)
    normb = norm(ψ)
    q = ψ / normb
    Q[:, 1] = q
    expiHt__ = zeros(Complex{Float64}, m)  # Initialize as zero vector
    residue = Inf

    for n in 1:max_dim
        # Arnoldi iteration
        v = -1im * t * H * q
        for j in 1:n
            h[j, n] = dot(Q[:, j], v)
            v -= h[j, n] * Q[:, j]
        end
        h[n + 1, n] = norm(v)

        # Check for convergence
        if abs(h[n + 1, n]) < tol
            Hn = h[1:n, 1:n]

            expiHt =
                Q[:, 1:n] *
                exp(Hn)[:, 1] *
                normb

            return expiHt, Dict(
                "n_iterations" => n,
                "residue" => 0.0,
            )
        end
        normb = norm(ψ)

        normb > 0 ||
            throw(ArgumentError("ψ must have nonzero norm"))

        # Update basis vector
        q = v / h[n + 1, n]
        Q[:, n + 1] = q

        # Compute matrix exponential approximation
        expiHt = Q[:, 1:n+1] * exp(h[1:n+1, 1:n+1])[:, 1] * normb

        # Compute residue
        residue = norm(expiHt - expiHt__)
        if residue < tol
            return expiHt, Dict("n_iterations" => n, "residue" => residue)
        else
            expiHt__ = expiHt
        end
    end

    return expiHt__, Dict("n_iterations" => max_dim, "residue" => residue)
end
