"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

# --- Representations of possibly-uncertain point sets ---
# The boundary points and the normals of a `BoundaryShape` (and the fields of a
# `BoundaryData`) can each be given in one of three forms:
#   1. a plain vector of (static) vectors                      — deterministic values;
#   2. a single `MvNormal` over the flattened values           — a joint Gaussian;
#   3. a vector of independent `MvNormal`s, one per point      — independent Gaussians.
# The helpers below convert between these forms: `struct_points` returns the mean as a
# `Vector{SVector{Dim}}`, `flat_points` returns the mean flattened, and the standard `cov`
# (from Statistics, re-exported by Distributions) is extended to return the covariance of
# the flattened values; plain vectors are deterministic, so their covariance is `0.0 * I`.

import Statistics: cov

cov(x::AbstractVector{<:AbstractVector}) = 0.0 * I
cov(x::AbstractVector{<:AbstractMvNormal}) = cat(cov.(x)...; dims = (1, 2))

struct_points(points::AbstractVector, Dim) = points
function struct_points(points::AbstractMvNormal, Dim)
    m = mean(points)
    return Vector(reinterpret(SVector{Dim, eltype(m)}, m))
end
struct_points(points::AbstractVector{<:AbstractMvNormal}, Dim) =
    [SVector{Dim, Float64}(mean(d)) for d in points]

flat_points(points::AbstractVector) = vcat(points...)
flat_points(points::AbstractMvNormal) = mean(points)
flat_points(points::AbstractVector{<:AbstractMvNormal}) = vcat(mean.(points)...)

# the fields are flattened and restructured exactly like the points
const struct_fields = struct_points
const flat_fields = flat_points

"""
    BoundaryShape{Dim,P,N} <: Shape{Dim}

The geometry of a domain, defined by a set of points on its boundary with no particular
order. The boundary points and the normals may each be deterministic (a vector of vectors)
or uncertain: a single `MvNormal` over the flattened values, or a vector of `MvNormal`s,
one per point. Use [`mean_points`](@ref) and [`mean_normals`](@ref) to obtain their means
as `Vector{SVector{Dim,Float64}}` regardless of the representation, and `cov` for the
covariance of the flattened values (`0.0 * I` when deterministic).

# Fields
- `boundary_points::P`: the points on the boundary, in any of the forms above.
- `normals::N`: the outward normals, aligned by index with the boundary points, in any of
  the forms above.
- `interior_points::Vector{SVector{Dim,Float64}}`: points in the interior of the domain,
  used to determine what is inside the domain, see `in(x::AbstractVector, ::BoundaryShape)`.

# Keyword constructor
    BoundaryShape(; boundary_points, normals = nothing, interior_points = nothing, Dim = 2)

Omitted `interior_points` default to the mean of the boundary points; omitted `normals` are
estimated with [`compute_outward_normals`](@ref). Deterministic normals are normalized to
unit length. `Dim` is deduced from the boundary points except when they are a single
`MvNormal` over the flattened values, where it cannot be and must be given.
"""
struct BoundaryShape{Dim, P, N} <: Shape{Dim}
    boundary_points::P
    normals::N
    interior_points::Vector{SVector{Dim,Float64}}
end

function BoundaryShape(;
        boundary_points = [zeros(Float64, 2)],
        normals = nothing,
        interior_points = nothing,
        Dim = 2
    )

    # the spatial dimension can be deduced except for a single flattened MvNormal
    actual_Dim = (boundary_points isa AbstractMvNormal) ? Dim : length(first(boundary_points))

    pts = struct_points(boundary_points, actual_Dim)

    interior = isnothing(interior_points) ? [mean(pts)] : interior_points
    interior = [SVector{actual_Dim, Float64}(p) for p in interior]

    actual_normals = isnothing(normals) ? compute_outward_normals(pts, interior) : normals

    # deterministic normals are normalized; uncertain ones are stored as given
    if actual_normals isa AbstractVector{<:AbstractVector{<:Real}}
        actual_normals = [SVector{actual_Dim, Float64}(n / norm(n)) for n in actual_normals]
    end

    return BoundaryShape(boundary_points, actual_normals, interior)
end

"""
    mean_points(shape::BoundaryShape)

The mean of the boundary points as a `Vector{SVector{Dim,Float64}}`, whichever of the three
representations (deterministic, joint `MvNormal`, per-point `MvNormal`s) they are stored in.
"""
mean_points(shape::BoundaryShape{Dim}) where Dim = struct_points(shape.boundary_points, Dim)

"""
    mean_normals(shape::BoundaryShape)

The mean of the normals as a `Vector{SVector{Dim,Float64}}`; see [`mean_points`](@ref).
"""
mean_normals(shape::BoundaryShape{Dim}) where Dim = struct_points(shape.normals, Dim)

"""
    BoundaryData(field_type::F; fields = , boundary_shape = , boundary_points = , normals = , interior_points = )

The physical data on a boundary: a [`BoundaryShape`](@ref) together with the values
`fields` of the physical field on the boundary points.

# Type parameters
- `F <: FieldType`: the physical field, e.g. displacement or traction.
- `Dim`: the spatial dimension.

# Fields
- `fieldtype::F`
- `boundary_shape::BS` where `BS <: BoundaryShape{Dim}`: the geometry of the boundary.
- `fields::FS`: `fields[i]` is the value of the physical field at boundary point `i`. Like
  the boundary points, the fields may be a plain vector of vectors, a single `MvNormal`
  over the flattened values (whose covariance is the sensor noise), or a vector of
  `MvNormal`s, one per point.

The keyword constructor either takes a ready-made `boundary_shape`, or builds one from
`boundary_points`, `normals` and `interior_points` (see [`BoundaryShape`](@ref)).
`BoundaryData` is itself a `Shape`: `in`, `bounding_box` and `points_in_shape` delegate to
its `boundary_shape`.
"""
struct BoundaryData{F <: FieldType, Dim, BS <: BoundaryShape{Dim}, FS} <: Shape{Dim}
    fieldtype::F
    boundary_shape::BS
    fields::FS
end

function BoundaryData(field_type::F;
        boundary_shape = nothing,
        boundary_points = [zeros(Float64, 2)],
        fields = nothing,
        normals = nothing,
        interior_points = nothing,
        Dim = 2
    ) where {F <: FieldType}

    shape = isnothing(boundary_shape) ?
        BoundaryShape(;
            boundary_points = boundary_points,
            normals = normals,
            interior_points = interior_points,
            Dim = Dim
        ) :
        boundary_shape

    resolved_fields = isnothing(fields) ? [zero(first(mean_points(shape)))] : fields

    return BoundaryData(field_type, shape, resolved_fields)
end

mean_points(bd::BoundaryData) = mean_points(bd.boundary_shape)
mean_normals(bd::BoundaryData) = mean_normals(bd.boundary_shape)

import MultipleScattering: name
name(shape::BoundaryShape) = "BoundaryShape"
name(shape::BoundaryData) = "BoundaryData"

bounding_box(shape::BoundaryShape) = Box(mean_points(shape))
bounding_box(bd::BoundaryData) = bounding_box(bd.boundary_shape)

import Base.in
function in(x::AbstractVector, shape::BoundaryShape)::Bool

    # find nearest interior point q to x. Then find the point p on the boundary which is closest to crossing the line through x and q. The point x is in the interior if norm(x - q) < norm(p - q)

    pts = mean_points(shape)

    dists = [sum((x - p) .^2) for p in shape.interior_points];
    q = shape.interior_points[argmin(dists)]

    vec = (x - q) ./ norm(x - q)

    # q + t .* vec == p
    dist_from_line = map(pts) do p
        t = dot(vec, p - q)
        t = (t < 0) ? 0.0 : t
        p_on_line = q + t .* vec
        norm(p_on_line - p)
    end
    p = pts[argmin(dist_from_line)]

    return norm(x - q) < norm(p - q)
end

in(x::AbstractVector, bd::BoundaryData)::Bool = in(x, bd.boundary_shape)

import Base.issubset
function issubset(shape::Union{BoundaryShape, BoundaryData}, box::Box)
    return issubset(bounding_box(shape), box)
end

"""
    issubset(box::Box, shape::BoundaryShape)

Returns true if the corners of the box are contained within the shape, false otherwise.
"""
function issubset(box::Box, shape::Union{BoundaryShape, BoundaryData})
    return all(c ∈ shape for c in corners(box))
end

function compute_outward_normals(boundary_points, interior_points;
        number_of_neighbours = min(2 * length(boundary_points[1]), max(1, length(boundary_points)-1))
    )

    # number of neighbors to use (at least Dim+1, at most n-1)
    k = number_of_neighbours

    pts = boundary_points
    n = length(pts)
    if n == 0
        return Vector{typeof(pts[1])}()
    end
    Dim = length(pts[1])
    T = eltype(pts[1])


    normals = Vector{typeof(pts[1])}(undef, n)
    # neighbour_distances = Vector{T}(undef, n)

    # for orientation choose the closest interior point to each boundary point
    for i in 1:n
        p = pts[i]

        # find k nearest neighbours (including the point itself)
        dists = [sum((p - q) .^2) for q in pts]
        idx = sortperm(dists)[1:min(k+1, n)]   # +1 because p itself is at distance 0
        neighbors = pts[idx]

        # neighbour_distances[i] = mean(dists[idx])

        # center data and compute covariance
        m = length(neighbors)
        μ = mean(neighbors)
        X = hcat([μ - q for q in neighbors]...) # Dim x k+1 matrix
        C = (X * transpose(X)) / max(1, m-1)   # covariance-like matrix (Dim x Dim)

        # eigen-decomposition: smallest eigenvalue eigenvector is normal to a (Dim-1)-manifold
        E = eigen(C)
        jmin = argmin(E.values)
        nvec = E.vectors[:, jmin]
        nvec = nvec / norm(nvec)

        # orient normal to point outward (away from interior)
        # pick closest interior point
        intdists = [sum((p - q) .^2) for q in interior_points]
        intp = interior_points[argmin(intdists)]
        interior_vector = intp .- p    # points from boundary point into interior
        # if dot(nvec, interior_vector) > 0 then nvec points inward -> flip
        if dot(nvec, interior_vector) > 0
            nvec = -nvec
        end

        # convert to same element type as pts
        normals[i] = convert(SVector{Dim,T}, nvec)
    end

    # neighbour_distance = mean(neighbour_distances) - std(neighbour_distances) / 2
    # return normals, neighbour_distance
    return normals
end
