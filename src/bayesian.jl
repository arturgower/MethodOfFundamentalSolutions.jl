# Bayesian Inference Routines for MFS Inverse Problems.
# Defining prior type 
abstract type ProbabilityDistribution end

struct GaussianDistribution{Dim, FieldDim} <: ProbabilityDistribution
    mean::Vector{SVector{Dim,Float64}}
    covariance::Union{AbstractMatrix{Float64}, UniformScaling{Float64}}
end

function GaussianDistribution(
    mean::Vector{<:AbstractVector{<:Real}},
    covariance::Union{AbstractMatrix{Float64}, UniformScaling{Float64}}
)
    Dim = length(first(mean))

    mean_svec = [SVector{Dim}(v) for v in mean]

    GaussianDistribution{Dim, length(mean)}(mean_svec, covariance)
end

function Distributions.MvNormal(dist::GaussianDistribution)
    # 1. Flatten the mean: MvNormal requires a 1D Vector, not a Vector of SVectors
    flat_mean = vcat(dist.mean...)
    
    # 2. Wrap the covariance in a PDMat to satisfy Distributions.jl
    # This handles both dense matrices and UniformScaling (diagonal/isotropic)
    pd_cov = isa(dist.covariance, UniformScaling) ? 
             PDMat(dist.covariance.λ * I(length(flat_mean))) : 
             PDMat(dist.covariance)
             
    return MvNormal(flat_mean, pd_cov)
end

function GaussianDistribution(mv_dist::MvNormal, Dim::Int)
    # 1. Extract the flat mean vector from the MvNormal distribution
    flat_mean = mean(mv_dist)
    
    # 2. Unflatten the mean back into a Vector of SVectors
    # We use `reinterpret` to chunk the flat array by `Dim` (e.g., pairs for 2D),
    # and `collect` to ensure it returns a standard Vector instead of a ReinterpretArray.
    T = eltype(flat_mean)
    structured_mean = collect(reinterpret(SVector{Dim, T}, flat_mean))
    
    # 3. Extract the raw covariance matrix
    # The `cov()` function safely unwraps the PDMat back into a standard Matrix
    cov_mat = cov(mv_dist)
    
    # 4. Construct your custom struct
    return GaussianDistribution(structured_mean, cov_mat)
end

#Helper function to compute the matrix (Cx).
#Shared across optimization and posterior routines.
#Works seamlessly for both scalar and matrix-valued Green's functions.

function compute_Cx(
    xb_flat::AbstractVector, 
    chi::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_x_input::AbstractMatrix, 
    M_func::F
) where {F}

    M_nom = M_func(xb_flat, chi)
    N, K = size(M_nom) # N: total measurement rows, K: total weight columns
    
    n_sensors = div(length(xb_flat), 2)
    d_m = div(N, n_sensors) # Automatically detects degrees of freedom per sensor (1 or 2)
    
    # Differentiate the fully assembled matrix M with respect to sensor positions
    jac_M = ForwardDiff.jacobian(x -> M_func(x, chi), xb_flat)
    
    Cx = zeros(eltype(chi), N, N)
    
    # Check if a full global covariance matrix was passed or just a 2x2 block
    is_full_sigma = size(Sigma_x_input, 1) == length(xb_flat)
    
    for r in 1:N, s in 1:N
        # Identify which physical sensors own row 'r' and row 's'
        i = div(r - 1, d_m) + 1
        j = div(s - 1, d_m) + 1
        
        # If sensors jitter independently, cross-sensor terms (i != j) are zero.
        if i == j 
            # Dynamically pull the correct 2x2 block if the full matrix was provided
            Sigma_x_block = if is_full_sigma
                idx_range = (2i-1):(2i)
                Sigma_x_input[idx_range, idx_range]
            else
                Sigma_x_input # Already a 2x2 block
            end

            J_r = zeros(eltype(chi), K, 2)
            J_s = zeros(eltype(chi), K, 2)
            
            for c in 1:K
                idx_r = r + (c - 1) * N
                idx_s = s + (c - 1) * N
                
                # Extract derivatives relative to the 2D coordinates of the owning sensor
                J_r[c, 1] = jac_M[idx_r, 2i - 1]; J_r[c, 2] = jac_M[idx_r, 2i]
                J_s[c, 1] = jac_M[idx_s, 2j - 1]; J_s[c, 2] = jac_M[idx_s, 2j]
            end
            
            # Account for cross-channel coupling (e.g., how x-jitter affects y-measurement)
            Cx[r, s] = tr(Sigma_a * J_r * Sigma_x_block * J_s')
        end
    end
    return Cx
end


function compute_Cx_analytical(
    xb_flat::AbstractVector, 
    chi::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_x_input::AbstractMatrix, 
    grad_M_func::F
) where {F}

    # grad_M_func returns a 3D tensor of shape (N_measurements, K_weights, 2)
    grad_M = grad_M_func(xb_flat, chi)
    N, K, _ = size(grad_M)
    
    n_sensors = div(length(xb_flat), 2)
    d_m = div(N, n_sensors) # Automatically detects measurement channels per sensor (1 or 2)
    
    Cx = zeros(eltype(chi), N, N)
    
    # Check if a full global covariance matrix was passed or just a 2x2 block
    is_full_sigma = size(Sigma_x_input, 1) == length(xb_flat)
    
    for r in 1:N, s in 1:N
        # Identify which physical sensors own row 'r' and row 's'
        i = div(r - 1, d_m) + 1
        j = div(s - 1, d_m) + 1
        
        # If sensors jitter independently, cross-sensor terms are zero
        if i == j 
            # Dynamically pull the correct 2x2 block if the full matrix was provided
            Sigma_x_block = if is_full_sigma
                idx_range = (2i-1):(2i)
                Sigma_x_input[idx_range, idx_range]
            else
                Sigma_x_input
            end

            # Extract the (K x 2) spatial derivative matrices directly
            J_r = grad_M[r, :, :] 
            J_s = grad_M[s, :, :]
            
            # Project 2D spatial error into measurement covariance
            Cx[r, s] = tr(Sigma_a * J_r * Sigma_x_block * J_s')
        end
    end
    return Cx
end


#Objective function: Compute the negative log-marginal likelihood.

function log_marginal_likelihood(
    chi::AbstractVector, 
    g::AbstractVector, 
    xb_flat::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_sensor::AbstractMatrix, 
    Sigma_x_block::AbstractMatrix, 
    M_func::F
) where {F}

    M_nom = M_func(xb_flat, chi)
    Cx = compute_Cx(xb_flat, chi, Sigma_a, Sigma_x_block, M_func)
    
    Sigma_gg = Sigma_sensor + Cx + M_nom * Sigma_a * M_nom'
    return 0.5 * (g' * (Sigma_gg \ g) + logdet(Sigma_gg) + length(g) * log(2π))
end

function log_marginal_likelihood(
    chi::AbstractVector, 
    g::AbstractVector, 
    xb_flat::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_sensor::AbstractMatrix, 
    Sigma_x_block::AbstractMatrix, 
    M_func::F, 
    grad_M_func::G
) where {F, G}

    M_nom = M_func(xb_flat, chi)
    Cx = compute_Cx_analytical(xb_flat, chi, Sigma_a, Sigma_x_block, grad_M_func)
    
    Sigma_gg = Sigma_sensor + Cx + M_nom * Sigma_a * M_nom'
    return 0.5 * (g' * (Sigma_gg \ g) + logdet(Sigma_gg) + length(g) * log(2π))
end
#Step 1: Hyperparameter Optimization.
#Finds the best parameter vector chi (e.g., basis properties).

function optimise_hyperparameters(
    g::AbstractVector, 
    xb_flat::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_sensor::AbstractMatrix, 
    Sigma_x_block::AbstractMatrix, 
    init_chi::AbstractVector, 
    M_func::F
) where {F}

    obj(chi) = log_marginal_likelihood(chi, g, xb_flat, Sigma_a, Sigma_sensor, Sigma_x_block, M_func)
    res = optimize(obj, init_chi, LBFGS(), autodiff=AutoForwardDiff())
    return Optim.minimizer(res)
end

function optimise_hyperparameters(
    g::AbstractVector, 
    xb_flat::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_sensor::AbstractMatrix, 
    Sigma_x_block::AbstractMatrix, 
    init_chi::AbstractVector, 
    M_func::F, 
    grad_M_func::G
) where {F, G}

    obj(chi) = log_marginal_likelihood(chi, g, xb_flat, Sigma_a, Sigma_sensor, Sigma_x_block, M_func, grad_M_func)
    # Optimization can still use forward-mode AD on the hyperparameter scalar 'chi' itself!
    res = optimize(obj, Vector(init_chi), LBFGS(), autodiff=AutoForwardDiff())
    return Optim.minimizer(res)
end

# Step 2: Calculate Posterior Coefficient Distribution.
# Computes the mean and covariance for the weights 'a'.

function compute_coefficient_posterior(
    g::AbstractVector, 
    xb_flat::AbstractVector, 
    chi::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_sensor::AbstractMatrix, 
    Sigma_x_block::AbstractMatrix, 
    M_func::F
) where {F}

    M_nom = M_func(xb_flat, chi)
    Cx = compute_Cx(xb_flat, chi, Sigma_a, Sigma_x_block, M_func)
    
    # Total effective noise matrix
    R_eff = Sigma_sensor + Cx
    R_inv = inv(R_eff)
    
    # Linear Gaussian update rules for coefficients
    Sigma_post = inv(M_nom' * R_inv * M_nom + inv(Sigma_a))
    mu_post = Sigma_post * (M_nom' * R_inv * g)
    
    return mu_post, Sigma_post
end

function compute_coefficient_posterior(
    g::AbstractVector, 
    xb_flat::AbstractVector, 
    chi::AbstractVector, 
    Sigma_a::AbstractMatrix, 
    Sigma_sensor::AbstractMatrix, 
    Sigma_x_block::AbstractMatrix, 
    M_func::F, 
    grad_M_func::G
) where {F, G}

    M_nom = M_func(xb_flat, chi)
    Cx = compute_Cx_analytical(xb_flat, chi, Sigma_a, Sigma_x_block, grad_M_func)
    
    R_eff = Sigma_sensor + Cx
    R_inv = inv(R_eff)
    
    Sigma_post = inv(M_nom' * R_inv * M_nom + inv(Sigma_a))
    mu_post = Sigma_post * (M_nom' * R_inv * g)
    
    return mu_post, Sigma_post
end

# Step 3: Reconstruct Full Continuous Field.
# Projects coefficient statistics back across the entire spatial grid domain.

function reconstruct_full_field(
    x_grid_flat::AbstractVector, 
    chi::AbstractVector, 
    mu_post::AbstractVector, 
    Sigma_post::AbstractMatrix, 
    phi_func::F; 
    return_full_cov::Bool=false
) where {F}

    # Map domain points through your custom function shape template
    Phi = phi_func(x_grid_flat, chi)
    
    # Mean field: μ_u = Φ * μ_post
    mu_u = Phi * mu_post
    
    if return_full_cov
        # Full covariance: Σ_u = Φ * Σ_post * Φᵀ
        Sigma_u = Phi * Sigma_post * Phi'
        return mu_u, Sigma_u
    else
        # Row-by-row variance computation for performance/memory savings
        variance_u = [dot(Phi[i, :], Sigma_post * Phi[i, :]) for i in 1:size(Phi, 1)]
        return mu_u, variance_u
    end
end

function _extract_bayesian_components(sim)
    bd = sim.boundary_data
    prior = sim.solver.prior
    
    # Extract Covariances
    Σ_a = prior.covariance
    Σ_sensor = bd.fields_covariance
    Σ_x = bd.boundary_points_covariance
    
    # Extract Vectors
    xb_flat = vcat(bd.boundary_points...)
    source_positions = vcat(sim.source_positions...)
    
    # Calculate effective measurement vector 'g'
    g = vcat(bd.fields...)
    g_particular = field(sim.medium, bd, sim.particular_solution)
    g = g - vcat(g_particular...)
    
    return g, xb_flat, source_positions, Σ_a, Σ_sensor, Σ_x
end

function _build_physics_closure(::Nothing, sim, xb_flat, chi_template)
  # 1. Handle the case where no gradient function is provided
    return nothing
end

function _build_physics_closure(func, sim, xb_flat, chi_template)
   # 2. FLAT SIGNATURE TEST (e.g., laplace_M)
    # Does the function accept (medium, flat_sensors, flat_sources)?
    if applicable(func, sim.medium, xb_flat, chi_template)
        return (xb, chi) -> func(sim.medium, xb, chi)
    end
    
   # 3. STRUCTURED SIGNATURE FALLBACK (e.g., your system_matrix)
    # If it isn't flat, we automatically adapt it to (structured_sources, medium, boundary_data)
    Dim = typeof(sim.medium).parameters[1]
    return (xb, chi) -> begin
        structured_sources = reinterpret(SVector{Dim, eltype(chi)}, chi)
        return func(structured_sources, sim.medium, sim.boundary_data)
    end
end


# Overload for Hyperparameter (Source Position) Optimization
function optimise_source_positions(
    sim, 
    system_matrix_function;
    gradient_system_matrix_function = nothing
) 
    g, xb_flat, init_source_positions, Σ_a, Σ_sensor, Σ_x = _extract_bayesian_components(sim)
    
    # Automatically build the correct closures based on the function passed!
    physics_matrix = _build_physics_closure(
        system_matrix_function, sim, xb_flat, init_source_positions
    )
    physics_gradient = _build_physics_closure(
        gradient_system_matrix_function, sim, xb_flat, init_source_positions
    )
    
    if isnothing(physics_gradient)
        return optimise_hyperparameters(
            g, xb_flat, Σ_a, Σ_sensor, Σ_x, init_source_positions, physics_matrix
        )
    else
        return optimise_hyperparameters(
            g, xb_flat, Σ_a, Σ_sensor, Σ_x, init_source_positions, physics_matrix, physics_gradient
        )
    end
end

function compute_coefficient_posterior(
    sim, 
    chi::AbstractVector, 
    system_matrix_function;
    gradient_system_matrix_function = nothing
) 
    g, xb_flat, _, Σ_a, Σ_sensor, Σ_x = _extract_bayesian_components(sim)
    
    # Automatically build the correct closures
    physics_matrix = _build_physics_closure(
        system_matrix_function, sim, xb_flat, chi
    )
    physics_gradient = _build_physics_closure(
        gradient_system_matrix_function, sim, xb_flat, chi
    )
    
    if isnothing(physics_gradient)
        return compute_coefficient_posterior(
            g, xb_flat, chi, Σ_a, Σ_sensor, Σ_x, physics_matrix
        )
    else
        return compute_coefficient_posterior(
            g, xb_flat, chi, Σ_a, Σ_sensor, Σ_x, physics_matrix, physics_gradient
        )
    end
end

