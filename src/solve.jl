"""
    AbstractSolver

Abstract type for different solution methods for the Method of Fundamental Solutions.
"""
abstract type AbstractSolver end

"""
    SolverOptions

Options shared by the solvers that can optimise the source positions
([`BayesianSolver`](@ref) and [`VariationalBayesianSolver`](@ref)), stored in their
`options` field. Each solver's keyword constructor accepts these as keywords directly.

# Fields
- `optimise_source_positions_flag::Bool`: optimise the source positions χ.
- `use_greens_gradient_analytical_flag::Bool`: use the analytic `greens_gradient` where
  available; otherwise finite differences are used.
- `update_geometry_flag::Bool`: update the boundary factor and re-center the boundary
  (used by [`VariationalBayesianSolver`](@ref)).
- `learn_prior_flag::Bool`: learn the prior (used by [`VariationalBayesianSolver`](@ref)).
- `max_iters::Int`: maximum number of iterations of the solver's outer loop.
- `source_position_iters::Int`: inner iterations per source-position optimisation step.
"""
struct SolverOptions
    optimise_source_positions_flag::Bool
    use_greens_gradient_analytical_flag::Bool
    update_geometry_flag::Bool
    learn_prior_flag::Bool
    max_iters::Int
    source_position_iters::Int
end

function SolverOptions(;
        optimise_source_positions_flag::Bool = false,
        use_greens_gradient_analytical_flag::Bool = true,
        update_geometry_flag::Bool = false,
        learn_prior_flag::Bool = true,
        max_iters::Int = 50,
        source_position_iters::Int = 5
    )
    return SolverOptions(
        optimise_source_positions_flag, use_greens_gradient_analytical_flag,
        update_geometry_flag, learn_prior_flag, max_iters, source_position_iters
    )
end

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
    BayesianSolver{P} <: AbstractSolver

Bayesian solver for MFS.

# Parameters
- `prior::P`: The prior distribution for the solution
- `options::SolverOptions`: the options shared with the other solvers, see
  [`SolverOptions`](@ref); the keyword constructor accepts them as keywords directly.
The solution is the posterior distribution over the coefficients given the boundary data and the prior.
"""
struct BayesianSolver{P<:ContinuousMultivariateDistribution} <: AbstractSolver
    prior::P
    options::SolverOptions
    gradient_tol::Float64
    objective_function_tol::Float64
end

function BayesianSolver(
    prior::ContinuousMultivariateDistribution;
    optimise_source_positions_flag::Bool = false,
    use_greens_gradient_analytical_flag::Bool = true,
    gradient_tol::Float64 = 1e-3,
    objective_function_tol::Float64 = 1e-4,
    max_iters::Int = 50
)
    options = SolverOptions(;
        optimise_source_positions_flag = optimise_source_positions_flag,
        use_greens_gradient_analytical_flag = use_greens_gradient_analytical_flag,
        max_iters = max_iters
    )
    return BayesianSolver{typeof(prior)}(prior, options, gradient_tol, objective_function_tol)
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

    # 3. Preallocate the flattened (N, K, Dim) array using the statically known NumType
    grad_M = zeros(NumType, N, K, Dim)
    
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

function solve(medium::P, bd::BoundaryData; kwargs... ) where P <: PhysicalMedium
    sim = Simulation(medium, bd; kwargs...)
    return solve(sim)
end

# Implement Tikhonov solver
function solve(sim::Simulation{TikhonovSolver{T}}) where T

    M = system_matrix(sim)

    forcing = flat_fields(sim.boundary_data.fields)
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

    @info "Solved the system with condition number $(condM), and with a relative error of the boundary data of $(relative_error), using the tolerance $(sim.solver.tolerance)"

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
    if sim.solver.options.optimise_source_positions_flag
        @info "Optimizing source positions..."
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
function source_positions(cloud::Union{BoundaryShape, BoundaryData}; relative_source_distance = 1.0)

    points = mean_points(cloud)
    normals = mean_normals(cloud)
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

    positions = map(eachindex(points)) do i
        points[i] + normals[i] .* source_distance
    end

    # Occasionally the normal direction is wrong. In which case, do not add a source inside the body!
    return filter(p -> p ∉ cloud, positions)
end
