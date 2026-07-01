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

"""
    BayesianSolver{T<:Real} <: AbstractSolver

Bayesian solver for MFS.

# Parameters
- `prior::P`: The prior distribution for the solution
The solution is the posterior distribution over the coefficients given the boundary data and the prior. 
"""
struct BayesianSolver{P<:ContinuousMultivariateDistribution} <: AbstractSolver
    prior::P
    optimise_source_positions_flag::Bool 
    use_greens_gradient_analytical_flag::Bool  
end

function BayesianSolver(
    prior::ContinuousMultivariateDistribution; 
    optimise_source_positions_flag::Bool = false, 
    use_greens_gradient_analytical_flag::Bool = true
)
    return BayesianSolver{typeof(prior)}(prior, optimise_source_positions_flag, use_greens_gradient_analytical_flag)
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

system_matrix(sim::Simulation) = system_matrix(sim.source_positions, sim.medium, sim.boundary_data)

function system_matrix(
    source_positions::AbstractVector{<:SVector{Dim}}, 
    medium::P, 
    bd::BoundaryData
) where {Dim, P<:PhysicalMedium{Dim}}

    points = bd.boundary_points isa AbstractMvNormal ? 
          struct_points(bd.boundary_points, Dim) : 
          bd.boundary_points
    normals = bd.outward_normals

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
        for i in eachindex(points), x in source_positions
    ]

    return Matrix(mortar(Ms))
end


function system_matrix_gradient(
    source_positions::AbstractVector{<:SVector{Dim}}, 
    medium::P, 
    bd::BoundaryData
) where {Dim, P<:PhysicalMedium{Dim}}

    points = bd.boundary_points isa AbstractMvNormal ? 
          struct_points(bd.boundary_points, Dim) : 
          bd.boundary_points
          
    normals = bd.outward_normals
    
    n_sensors = length(points)
    n_sources = length(source_positions)

    T_points = eltype(eltype(points))
    T_sources = eltype(eltype(source_positions))
    NumType = promote_type(T_points, T_sources)

    # --- THE MAGIC: DYNAMIC ADAPTATION ---
    # 1. Run a single test evaluation to see what the physics kernel outputs
    r_vec_test = points[1] - source_positions[1]
    normal_test = SVector{Dim, NumType}(normals[1])
    G_sample = greens_gradient(bd.fieldtype, medium, r_vec_test, normal_test)
    
    # 2. Extract its exact data type and shape 
    # Laplace will report shape `(2,)`. Elasticity might report `(2, 2)`.
    OutType = eltype(G_sample)
    out_shape = size(G_sample)

    # 3. Preallocate an N-dimensional array using the exact shape reported!
    # For Laplace, this builds an Array{OutType, 3} of size (N, M, 2).
    # For Elasticity, it builds an Array{OutType, 4} of size (N, M, 2, 2).
    grad_M = Array{OutType}(undef, n_sensors, n_sources, out_shape...)

    # --- MATRIX ASSEMBLY ---
    for j in 1:n_sources
        x = source_positions[j]
        for i in 1:n_sensors
            r_vec = points[i] - x
            normal_vec = SVector{Dim, NumType}(normals[i])
            
            # Evaluate the physics
            G_grad = greens_gradient(bd.fieldtype, medium, r_vec, normal_vec)
            
            # CartesianIndices automatically iterates over whatever shape G_grad is!
            for K in CartesianIndices(G_grad)
                grad_M[i, j, K] = G_grad[K]
            end
        end
    end

    return grad_M
end

system_matrix_gradient(sim::Simulation) = system_matrix_gradient(sim.source_positions, sim.medium, sim.boundary_data)

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

function solve(
    sim::Simulation{<:BayesianSolver{<:AbstractMvNormal}, Dim}
    ) where {Dim}
    
    # 1. Determine Source Positions (chi)
    if sim.solver.optimise_source_positions_flag
        println("Optimizing source positions...")
        best_source_positions = optimise_source_positions(sim)
    else
        best_source_positions = vcat(sim.source_positions...)
    end

    # 2. Compute Posterior Coefficients
    μ_post, Σ_post = compute_coefficient_posterior(
            sim, best_source_positions
        )
    
    new_source_positions = [
    SVector{Dim, Float64}(best_source_positions[i : i + Dim - 1]) 
    for i in 1:Dim:length(best_source_positions)
    ]
    # 3. Return the solution
    return FundamentalSolution(
        sim.medium; 
        positions = new_source_positions,
        coefficients = μ_post, 
        coefficients_covariance = Σ_post,
        particular_solution = sim.particular_solution
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
