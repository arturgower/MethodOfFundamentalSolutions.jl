"""
    AbstractSolver

Abstract type for different solution methods for the Method of Fundamental Solutions.
"""
abstract type AbstractSolver end

"""
    ParticularSolution

A type used to specify the type of particular solution to add to the boundary data.
"""
abstract type ParticularSolution end

struct NoParticularSolution <: ParticularSolution end

"""
    ParticularSolutions

A vector of [`ParticularSolution`](@ref)s whose fields are summed. Accepted anywhere a
single `ParticularSolution` is (e.g. as the `particular_solution` of a [`Simulation`](@ref)).
"""
const ParticularSolutions = AbstractVector{<:ParticularSolution}

# a single ParticularSolution or a vector of them (whose fields are summed)
const AnyParticularSolution = Union{ParticularSolution, ParticularSolutions}

"""
    TikhonovSolver{T<:Real} <: AbstractSolver

Tikhonov regularization solver for MFS.

# Parameters
- `λ::T`: Regularization parameter see equation below (default:-1 means will be auto-computed from tolerance)
- `tolerance::T`: Tolerance for auto-computing λ if not specified

The solution x minimises
```math
x = \\argmin_x \\|A x - b\\|^2 + \\lambda \\|x\\|^2
```
"""
struct TikhonovSolver{T<:Real} <: AbstractSolver
    λ::T
    tolerance::T
    function TikhonovSolver(; λ::Real = - one(Float64), tolerance::Real = 1e-11)
        T = typeof(λ)
        new{T}(λ, tolerance)
    end
end

"""
    Simulation{S,Dim,P,PS,BD}

An MFS problem: a `medium`, the `boundary_data` to be matched, a `solver`, and the MFS
`source_positions` (the columns of the system matrix). The angular frequency ω is taken
from the `medium`, not stored separately.

`boundary_data` may be a single [`BoundaryData`](@ref), or a `Tuple` of `BoundaryData` that
share the same `source_positions` but impose different boundary conditions on (possibly)
different points — e.g. traction on part of the boundary and displacement on another. Each
`BoundaryData` in the tuple contributes its own block of rows to the system matrix, stacked
in the order given.

`particular_solution` may be a single [`ParticularSolution`](@ref) or a vector of them
([`ParticularSolutions`](@ref)), in which case their fields are summed.
"""
struct Simulation{S <: AbstractSolver, Dim, P<:PhysicalMedium{Dim}, PS <: AnyParticularSolution, BD}
    solver::S
    medium::P
    boundary_data::BD
    particular_solution::PS
    source_positions::Vector{SVector{Dim,Float64}}
end

# treat a lone BoundaryData as a one-element tuple so the assembly code is written once
_as_tuple(bd::BoundaryData) = (bd,)
_as_tuple(bds::Tuple) = bds

function Simulation(medium::P, bd::BoundaryData{<:FieldType,Dim};
        solver::AbstractSolver = TikhonovSolver(),
        particular_solution::AnyParticularSolution = NoParticularSolution(),
        source_positions = source_positions(bd; relative_source_distance = 1.2),
    ) where {Dim, P <: PhysicalMedium{Dim}}

    sp = [SVector{Dim, Float64}(p) for p in source_positions]
    return Simulation(solver, medium, bd, particular_solution, sp)
end

function Simulation(medium::P, bds::Tuple{Vararg{BoundaryData{<:FieldType,Dim}}};
        solver::AbstractSolver = TikhonovSolver(),
        particular_solution::AnyParticularSolution = NoParticularSolution(),
        source_positions = source_positions(bds; relative_source_distance = 1.2),
    ) where {Dim, P <: PhysicalMedium{Dim}}

    sp = [SVector{Dim, Float64}(p) for p in source_positions]
    return Simulation(solver, medium, bds, particular_solution, sp)
end

system_matrix(sim::Simulation) = system_matrix(sim.source_positions, sim.medium, sim.boundary_data)

# a tuple of BoundaryData stacks its blocks vertically; all blocks share the same columns
function system_matrix(source_pos::AbstractVector{<:SVector}, medium::PhysicalMedium, bds::Tuple)
    return reduce(vcat, map(bd -> system_matrix(source_pos, medium, bd), bds))
end

function system_matrix(
    source_pos::AbstractVector{<:SVector{Dim}}, 
    medium::P, 
    bd::BoundaryData
) where {Dim, P<:PhysicalMedium{Dim}}

    points = mean_points(bd)
    normals = mean_normals(bd)

    # The comprehension block automatically handles the return type
    Ms = [
        begin
            # 1. This distance vector naturally promotes to Dual during optimization
            r_vec = points[i] - x 
            
            # 2. Extract whatever type r_vec ended up being (Float64 or Dual)
            NumType = eltype(r_vec) 
            
            # 3. Cast the boundary normal to perfectly match the Dual/Float64 type
            normal_vec = SVector{Dim, NumType}(normals[i])
            
            # 4. Call greens. Now r_vec and normal_vec share the exact same type!
            greens(bd.fieldtype, medium, r_vec, normal_vec)
        end
        for i in eachindex(points), x in source_pos
    ]

    Ms = (typeof(Ms[1]) <: AbstractMatrix) ? Matrix(mortar(Ms)) : Ms
    
    return Ms
end


function system_matrix_gradient(
    source_positions::AbstractVector{<:SVector{Dim}}, 
    medium::P, 
    bd::BoundaryData
) where {Dim, P<:PhysicalMedium{Dim}}

    points = mean_points(bd)
    normals = mean_normals(bd)

    n_sensors = length(points)
    n_sources = length(source_positions)

    T_points = eltype(eltype(points))
    T_sources = eltype(eltype(source_positions))
    NumType = promote_type(T_points, T_sources)

    # 1. Run test evaluation
    r_vec_test = points[1] - source_positions[1]
    normal_test = SVector{Dim, NumType}(normals[1])
    G_sample = greens_gradient(bd.fieldtype, medium, r_vec_test, normal_test)
    
    out_shape = size(G_sample)

    # 2. Dynamically determine Degrees of Freedom (DOF) based on kernel shape
    d_m_out = length(out_shape) == 1 ? 1 : out_shape[1]
    d_m_in  = length(out_shape) == 1 ? 1 : out_shape[2]
    
    N = n_sensors * d_m_out
    K = n_sources * d_m_in

    # 3. Preallocate the flattened (N, K, Dim) array. The element type follows the kernel as
    #    well as the geometry: complex physics (e.g. acoustics) has a complex gradient even
    #    though the points and sources are real.
    grad_M = zeros(promote_type(NumType, eltype(G_sample)), N, K, Dim)
    
    # --- MATRIX ASSEMBLY ---
    for j in 1:n_sources
        x = source_positions[j]
        for i in 1:n_sensors
            r_vec = points[i] - x
            normal_vec = SVector{Dim, NumType}(normals[i])
            
            G_grad = greens_gradient(bd.fieldtype, medium, r_vec, normal_vec)
            
            # 4. Map local tensor to global flat indices
            if length(out_shape) == 1
                # Scalar physics (Laplace)
                for k in 1:Dim
                    grad_M[i, j, k] = G_grad[k]
                end
            else
                # Matrix/Tensor physics (Elasticity)
                for i_ch in 1:d_m_out, j_ch in 1:d_m_in, k in 1:Dim
                    row = d_m_out * (i - 1) + i_ch
                    col = d_m_in  * (j - 1) + j_ch
                    grad_M[row, col, k] = G_grad[i_ch, j_ch, k]
                end
            end
        end
    end

    return grad_M
end

system_matrix_gradient(sim::Simulation) = system_matrix_gradient(sim.source_positions, sim.medium, sim.boundary_data)

# a tuple of BoundaryData stacks its gradient blocks along the row dimension
function system_matrix_gradient(source_pos::AbstractVector{<:SVector}, medium::PhysicalMedium, bds::Tuple)
    blocks = map(bd -> system_matrix_gradient(source_pos, medium, bd), bds)
    return reduce((a, b) -> cat(a, b; dims = 1), blocks)
end

function solve(medium::P, bd::BoundaryData; kwargs... ) where P <: PhysicalMedium
    sim = Simulation(medium, bd; kwargs...)
    return solve(sim)
end

"""
    boundary_forcing(sim::Simulation)

The right-hand side that `sim`'s fundamental solution must match on the boundary: the
(flattened) boundary `fields` minus the contribution of the `particular_solution`. For a
tuple `boundary_data` the blocks are stacked in the same order as [`system_matrix`](@ref).
"""
function boundary_forcing(sim::Simulation)
    bds = _as_tuple(sim.boundary_data)
    fields = reduce(vcat, map(bd -> flat_fields(bd.fields), bds))
    particular = reduce(vcat,
        map(bd -> vcat(field(sim.medium, bd, sim.particular_solution)...), bds))
    return fields - particular
end

# Regularised least-squares solve of `M * coes = forcing`, shared by the single-domain and
# transmission solvers. Returns the coefficients, the relative boundary error and cond(M).
function tikhonov_solve(M::AbstractMatrix, forcing::AbstractVector, solver::TikhonovSolver)
    condM = cond(M)
    sqrtλ = solver.λ < zero(eltype(solver.λ)) ?
        condM * sqrt(solver.tolerance) :
        sqrt(solver.λ)

    bigM = [M; sqrtλ * I]
    coes = bigM \ [forcing; zeros(size(M)[2])]

    relative_error = norm(M * coes - forcing) / norm(forcing)

    return coes, relative_error, condM
end

# Implement Tikhonov solver
function solve(sim::Simulation{TikhonovSolver{T}}) where T

    M = system_matrix(sim)
    forcing = boundary_forcing(sim)

    coes, relative_error, condM = tikhonov_solve(M, forcing, sim.solver)

    @info "Solved the system with condition number $(condM), and with a relative error of the boundary data of $(relative_error), using the tolerance $(sim.solver.tolerance)"

    return FundamentalSolution(sim.medium;
        positions = sim.source_positions,
        coefficients = coes,
        particular_solution = sim.particular_solution,
        relative_boundary_error = relative_error
    )
end


"""
    source_positions(cloud::BoundaryData; α=1.0)

Return source positions for MFS from some `BoundaryData`.

- α: scale factor for the distance of the source from the boundary d = α * h, where h is the average distance between consecutive points on the boundary.
"""
function source_positions(cloud::Union{BoundaryShape, BoundaryData}; relative_source_distance = 1.0)

    points = mean_points(cloud)
    normals = mean_normals(cloud)

    source_distance = _boundary_spacing(points) * relative_source_distance

    positions = map(eachindex(points)) do i
        points[i] + normals[i] .* source_distance
    end

    # Occasionally the normal direction is wrong. In which case, do not add a source inside the body!
    return filter(p -> p ∉ cloud, positions)
end

"""
    source_positions(bds::Tuple; relative_source_distance = 1.0)

Source positions shared by a tuple of [`BoundaryData`](@ref). The boundary conditions in
`bds` all use the same MFS sources (columns), so the sources are placed relative to the
union of their boundaries. See the single-argument [`source_positions`](@ref).
"""
function source_positions(bds::Tuple{Vararg{Union{BoundaryShape, BoundaryData}}}; relative_source_distance = 1.0)

    _interior(bd::BoundaryData) = bd.boundary_shape.interior_points
    _interior(shape::BoundaryShape) = shape.interior_points

    points = reduce(vcat, map(mean_points, bds))
    normals = reduce(vcat, map(mean_normals, bds))
    interiors = reduce(vcat, map(_interior, bds))

    combined = BoundaryShape(
        boundary_points = points,
        normals = normals,
        interior_points = interiors
    )

    return source_positions(combined; relative_source_distance = relative_source_distance)
end

"""
    grid_source_positions(bd::BoundaryData; n = 15, scale = 2.0, clearance = 1.0)

Candidate MFS source positions "everywhere": a regular `n × n` grid covering the bounding
box of the boundary enlarged by `scale`, keeping only points outside the domain and further
than `clearance` times the average boundary spacing from the boundary. Intended as an
overcomplete set of candidates for a [`VariationalBayesianSolver`](@ref), which selects the
useful sources one at a time. See also [`source_positions`](@ref), which instead places one
source behind every boundary point.
"""
function grid_source_positions(bd::BoundaryData{F, 2}; n::Int = 15, scale::Real = 2.0, clearance::Real = 1.0) where F
    pts = mean_points(bd)

    xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
    centre = SVector((minimum(xs) + maximum(xs)) / 2, (minimum(ys) + maximum(ys)) / 2)
    halfwidth = SVector(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) ./ 2

    spacing = _boundary_spacing(pts)

    grid = [
        centre + SVector(2u - 1, 2v - 1) .* (scale .* halfwidth)
    for u in LinRange(0, 1, n), v in LinRange(0, 1, n)]

    return filter(vec(grid)) do p
        p ∉ bd && minimum(norm(p - q) for q in pts) > clearance * spacing
    end
end
