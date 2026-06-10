"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

"""
    BoundaryData(field_type::F; boundary_points = , fields = , interior_points = , outward_normals = ) where {F <: FieldType}

A [`Shape`](@ref) defined by a set of points on the boundary with no particular order in the data. 

Type parameters
- F <: FieldType: the physical field witehr displacement or traction.
- Dim: Integer-sized compile-time dimension of the spatial dimension

Fields
- `interior_points::Vector{SVector{Dim,Float64}}`
    A vector of spatial points in the interior of the domain used to determine what is inside the domain. See
@doc in(x::AbstractVector, cloud::PointCloud).
- `points::Vector{SVector{Dim,Float64}}`
    A vector of spatial points on the boundary
- `fields::Vector{SVector{Dim,Float64}}`
    `fields[i]` is the value of the physical field at `points[i]`

Notes
- It is expected that `length(points) == length(fields)` and that entries are aligned by index.
"""
struct BoundaryData{F <: FieldType, Dim, FieldDim} <: Shape{Dim}
    fieldtype::F
    boundary_points::Vector{SVector{Dim,Float64}}
    fields::Vector{SVector{FieldDim,Union{Float64, ComplexF64}}}
    outward_normals::Vector{SVector{Dim,Float64}}
    interior_points::Vector{SVector{Dim,Float64}}
    fields_covariance::Union{AbstractMatrix{Float64}, UniformScaling{Float64}}
    boundary_points_covariance::Union{AbstractMatrix{Float64}, UniformScaling{Float64}}
end

function BoundaryData(field_type::F; 
        boundary_points::Vector = [zeros(Float64,2)], 
        fields::Vector{FV} = [p .* 0.0 for p in boundary_points],
        interior_points::Vector = [mean(boundary_points)],
        outward_normals::Vector = outward_normals(boundary_points,interior_points),
        fields_covariance::Union{AbstractMatrix{Float64}, UniformScaling{Float64}} = I,
        boundary_points_covariance::Union{AbstractMatrix{Float64}, UniformScaling{Float64}} = I,
    ) where {F <: FieldType, FV <: AbstractVector}

    Dim = length(boundary_points[1])
    FieldDim = length(fields[1])

    # All outward normals should have unit length
    outward_normals = [n / norm(n) for n in outward_normals]

    BoundaryData{F,Dim,FieldDim}(field_type, boundary_points, fields, outward_normals, interior_points, fields_covariance, boundary_points_covariance)
end

import MultipleScattering: name
name(shape::BoundaryData) = "BoundaryData"

bounding_box(cloud::BoundaryData) = Box(cloud.boundary_points)

# function BoundaryData(medium::Ph, bd::BoundaryData, psol::P) where {Ph <: PhysicalMedium,P <: ParticularSolution}
#     fs = field(medium,bd,psol)
   
#     bd_particular = @set bd.fields = bd.fields + fs
#     return bd_particular
# end

import Base.in
function in(x::AbstractVector, cloud::BoundaryData)::Bool

    # find nearest interior point q to x. Then  find the point p on the boundary which is closest to crossing the line through x and q. The point x is in the interior if norm(x - q) < norm(p - q)

    dists = [sum((x - p) .^2) for p in cloud.interior_points];
    q = cloud.interior_points[argmin(dists)]

    vec = (x - q) ./ norm(x - q)

    # q + t .* vec == p
    dist_from_line = map(cloud.boundary_points) do p
        t = dot(vec, p - q)
        t = (t < 0) ? 0.0 : t
        p_on_line = q + t .* vec
        norm(p_on_line - p)
    end
    p = cloud.boundary_points[argmin(dist_from_line)]

    inside = norm(x - q) < norm(p - q) ? true : false
    return inside
end

import Base.issubset
function issubset(cloud::BoundaryData, box::Box)
    cloud_box = bounding_box(cloud)
    return issubset(cloud_box,box)
end

"""
    issubset(box::Box, poly::PointCloud)

Returns true if the corners of the box are contained within polygon, false otherwise.
"""
function issubset(box::Box, cloud::BoundaryData)
    return all(c ∈ cloud for c in corners(box))
end

function outward_normals(boundary_points, interior_points; 
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
