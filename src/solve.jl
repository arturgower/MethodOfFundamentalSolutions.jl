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

struct Simulation{S <: AbstractSolver, Dim, P<:PhysicalMedium{Dim}, PS <:ParticularSolution, BD <: BoundaryData}
    solver::S
    medium::P
    boundary_data::BD
    particular_solution::PS
    source_positions::Vector{SVector{Dim,Float64}}
    ω::Float64
end

function Simulation(medium::P, bd::BD; 
        solver::S = TikhonovSolver(),
        particular_solution::PS = NoParticularSolution(),
        source_positions = source_positions(bd; relative_source_distance = 1.2),
        ω::Float64 = 2pi * 1.0 
    ) where {
        S <: AbstractSolver, Dim, 
        P <: PhysicalMedium{Dim}, PS <: ParticularSolution, 
        BD <: BoundaryData{<:FieldType,Dim}
    }

    return Simulation{S,Dim,P,PS,BD}(solver, medium, bd, particular_solution, source_positions, ω)
end

system_matrix(sim::Simulation) = system_matrix(sim.source_positions, sim.medium, sim.boundary_data; ω=sim.ω)

function system_matrix(source_positions::Vector{SVector{Dim,Float64}}, medium::P, bd::BoundaryData; ω::Float64=2pi) where {Dim,P<:PhysicalMedium{Dim}}

    points = bd.boundary_points

    Ms = [
        greens(bd.fieldtype, medium, points[i] - x, bd.outward_normals[i]; ω=ω)    
    for i in eachindex(points), x in source_positions]

    return mortar(Ms)
end

function solve(medium::P, bd::BoundaryData; kwargs... ) where P <: PhysicalMedium
    sim = Simulation(medium, bd; kwargs...)
    return solve(sim)
end

# Implement Tikhonov solver
function solve(sim::Simulation{TikhonovSolver{T}}) where T

    M = system_matrix(sim)

    forcing = vcat(sim.boundary_data.fields...)
    forcing_particular = field(sim.medium, sim.boundary_data, sim.particular_solution)
    forcing = forcing - vcat(forcing_particular...)
    

    # Tikinov solution
    condM = cond(M)
    sqrtλ = if sim.solver.λ < zero(eltype(sim.solver.λ)) 
        condM * sqrt(sim.solver.tolerance)
    else sqrt(sim.solver.λ)
    end

    bigM = [M; sqrtλ * I];
    coes = bigM \ [forcing; zeros(size(M)[2])]

    relative_error = norm(M * coes - forcing) / norm(forcing)

    println("Solved the system with condition number:$(condM), and with a relative error of boundary data: $(norm(M * coes - forcing) / norm(forcing)) with a tolerance of $(sim.solver.tolerance)")

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
function source_positions(cloud::BoundaryData; relative_source_distance = 1.0) 

    points = cloud.boundary_points 
    len = points |> length
    
    # Note this could be calculated at the same time as the outward normals. But that would make the code quite ugly!
    # Sample just a few number of points to approximate the distance between neighbours
    sampled_rng = LinRange(1,len, min(6,len)) .|> round .|> Int
    neighbors_dists = map(points[sampled_rng]) do p 
        dists = [norm(p - q) for q in points]
        idx = sortperm(dists)[2:min(3, len)]
        mean(dists[idx])
    end
    source_distance = mean(neighbors_dists) * relative_source_distance
    
    positions = map(cloud.boundary_points |> eachindex) do i
        cloud.boundary_points[i] + cloud.outward_normals[i] .* source_distance 
    end

    # Occasionally the normal direction is wrong. In which case, do not add a source inside the body!
    return filter(p -> p ∉ cloud, positions)
end
