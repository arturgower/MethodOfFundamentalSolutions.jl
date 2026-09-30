# Cholesky of a symmetric positive definite matrix, retried with a small diagonal jitter when
# round-off makes the factorization fail. LinearAlgebra offers the non-throwing
# `cholesky(A; check = false)` used here, but no jittered retry.
function _safe_cholesky(A::AbstractMatrix)
    S = Symmetric(Matrix(A))
    C = cholesky(S; check = false)
    issuccess(C) && return C
    jitter = 1e-12 * max(1.0, maximum(abs, diag(S)))
    return cholesky(Symmetric(S + jitter * I))
end

# A left square root W of a precision matrix P, so that P = WᵀW. The precision is given
# either as a vector of reciprocal variances (independent noise, the fast path) or as a full
# precision matrix.
_whitener(precision::AbstractVector) = Diagonal(sqrt.(precision))
_whitener(precision::AbstractMatrix) = Matrix(transpose(_safe_cholesky(precision).L))

# Average distance between neighbouring boundary points, estimated from a few sampled points:
# the natural length scale of a point cloud, used both to place sources and to keep them clear
# of the boundary.
function _boundary_spacing(points)
    len = length(points)
    sampled_rng = LinRange(1, len, min(6, len)) .|> round .|> Int
    return mean(map(points[sampled_rng]) do p
        dists = [norm(p - q) for q in points]
        idx = sortperm(dists)[2:min(3, len)]
        mean(dists[idx])
    end)
end

function interior_points_along_coordinate(points; offset_percent = 0.1,
    coordinate_number = 2, segments = 10)

    if offset_percent < 0.0 || offset_percent > 1.0
        throw(ArgumentError("offset_percent must be between 0.0 and 1.0"))
    end
    
    cs = [xy[coordinate_number] for xy in points]

    minc = minimum(cs)
    maxc = maximum(cs)

    edges = collect(range(minc, stop=maxc, length=segments+1))

    # collect indices for each y-bin
    bin_indices = [Int[] for _ in 1:segments]
    
    for (idx, y) in enumerate(cs)
        b = clamp(searchsortedlast(edges, y), 1, segments)
        push!(bin_indices[b], idx)
    end

    bin_indices = filter(v -> !isempty(v), bin_indices)

    interior_points = map(bin_indices) do inds
        p = mean(points[inds])
        
        inner_low = minc + offset_percent * (maxc - minc)
        inner_high = maxc - offset_percent * (maxc - minc)
        pc = p[coordinate_number]
        t = clamp((pc - minc) / (maxc - minc), 0.0, 1.0)

        p = @set p[coordinate_number] = inner_low + t * (inner_high - inner_low)
        p
    end

    return interior_points
end
