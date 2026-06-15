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

struct BayesianSolver{P<:Prior} <: AbstractSolver
    prior::P
end

function BayesianSolver(prior::P) where {P<:Prior}
    BayesianSolver{P}(prior)
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

function system_matrix(source_positions::Vector{SVector{Dim,Float64}}, medium::P, bd::BoundaryData) where {Dim,P<:PhysicalMedium{Dim}}

    points = bd.boundary_points

    Ms = [
        greens(bd.fieldtype, medium, points[i] - x, bd.outward_normals[i])    
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

function solve(sim::Simulation{<:BayesianSolver{<:GaussianPrior}},
    system_matrix_function::Function; 
    gradient_system_matrix_function::Union{Function, Nothing} = nothing,
    optimize_source_positions::Bool = false) 
    
    bd = sim.boundary_data
    prior = sim.solver.prior
    
    μ_a = vcat(prior.mean...)
    Σ_a = prior.covariance_matrix
    

    
    g = vcat(sim.boundary_data.fields...)
    g_particular = field(sim.medium, sim.boundary_data, sim.particular_solution)
    g = g - vcat(g_particular...)
    Σ_sensor = bd.fields_covariance
    
    x0_sensors = vcat(bd.boundary_points...)
    Σ_x=bd.boundary_points_covariance

    init_source_positions = vcat(sim.source_positions...)

    if optimize_source_positions
        if isnothing(gradient_system_matrix_function)
            # Calculate WITHOUT the analytical gradient function
       
            best_source_positions = optimize_hyperparameters(
                                    g, 
                                    x0_sensors, 
                                    Σ_a, 
                                    Σ_sensor, 
                                    Σ_x, 
                                    init_source_positions, 
                                    system_matrix_function
                            )
        # Add your source position optimization logic here
        
        else
            best_source_positions = optimize_hyperparameters(
                                     g, 
                                     x0_sensors, 
                                     Σ_a, 
                                     Σ_sensor, 
                                     Σ_x, 
                                     init_source_positions, 
                                     system_matrix_function, 
                                     gradient_system_matrix_function
                                    )
        # Add your source position optimization logic here
            println("Solving with gradient...")
            # You can call it here: gradient_system_matrix_function(...)
        end
    else
        best_source_positions = init_source_positions 
    end

    μ_post, Σ_post = compute_coefficient_posterior(g, x0_sensors, best_source_positions, Σ_a, Σ_sensor, Σ_x, laplace_M, laplace_grad_M)

    return FundamentalSolution(sim.medium; 
        positions = best_source_positions,
        coefficients_mean = μ_post, 
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
