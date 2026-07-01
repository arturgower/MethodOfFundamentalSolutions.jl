# Bayesian Inference Routines for MFS Inverse Problems.

#Helper function to compute the matrix (Cx).
#Shared across optimization and posterior routines.
#Works seamlessly for both scalar and matrix-valued Green's functions.

function extract_bayesian_components(sim)
    bd = sim.boundary_data
    prior = sim.solver.prior
    
    # Extract Covariances using the helper functions
    Σ_a = cov(prior)
    Σ_sensor = cov_fields(bd.fields)
    Σ_x = cov_points(bd.boundary_points)
    
    # Extract Flat Vectors using the helper functions
    xb_flat = flat_points(bd.boundary_points)
    source_positions = vcat(sim.source_positions...)
    
    # Calculate effective measurement vector 'g'
    g = flat_fields(bd.fields)
    g_particular = field(sim.medium, bd, sim.particular_solution)
    g = g - vcat(g_particular...)
    
    return g, xb_flat, source_positions, Σ_a, Σ_sensor, Σ_x
end


function geometric_covariance(
    sim, 
    structured_chi::Vector{SVector{2, T}};
    h::Real = 1e-5
) where T
    # Clean extraction of baseline variables
    _, xb_flat, _, Sigma_a, _, Sigma_x_input = extract_bayesian_components(sim)
    
    n_xb = length(xb_flat)
    n_sensors = div(n_xb, 2)
    
    if sim.solver.use_greens_gradient_analytical_flag
        
        # ---------------------------------
        #    Analytical Gradient Branch
        # ---------------------------------
        #grad_M = grad_M_func(xb_flat, chi)
        grad_M = system_matrix_gradient(structured_chi, sim.medium, sim.boundary_data)
        N, K, _ = size(grad_M)
        d_m = div(N, n_sensors)
        
        Cx = zeros(T, N, N)
        is_full_sigma = size(Sigma_x_input, 1) == n_xb
        
        for r in 1:N, s in 1:N
            i = div(r - 1, d_m) + 1
            j = div(s - 1, d_m) + 1
            
            if i == j 
                Sigma_x_block = is_full_sigma ? Sigma_x_input[(2i-1):(2i), (2i-1):(2i)] : Sigma_x_input
                J_r = grad_M[r, :, :] 
                J_s = grad_M[s, :, :]
                
                Cx[r, s] = tr(Sigma_a * J_r * Sigma_x_block * J_s')
            end
        end
        return Cx
        
    else
        
        # ---------------------------------
        #    Finite Difference Branch
        # ---------------------------------
        #M_nom = M_func(xb_flat, chi)
        M_nom = system_matrix(structured_chi, sim.medium, sim.boundary_data)
        N, K = size(M_nom)
        d_m = div(N, n_sensors)
        
        jac_M = zeros(T, N * K, n_xb)
        # 1. Allocate the structured boundary points ONCE outside the loop.
        # This gives system_matrix the exact heap-allocated Vector it demands.
        xb_structured = Vector(reinterpret(SVector{2, eltype(xb_flat)}, xb_flat))

        for v in 1:n_xb
            # 2. Map the flat index (v) to the SVector index (s_idx) and the coordinate (x=1, y=2)
            s_idx = div(v - 1, 2) + 1
            coord = mod1(v, 2)
    
            # Save the original SVector
            orig_svec = xb_structured[s_idx]
            
            # --- Forward Step ---
            # Create a new SVector on the stack (zero heap allocations) and overwrite the slot
            if coord == 1
                xb_structured[s_idx] = SVector(orig_svec[1] + h, orig_svec[2])
            else
                xb_structured[s_idx] = SVector(orig_svec[1], orig_svec[2] + h)
            end
            
            bd_fw = BoundaryData(sim.boundary_data.fieldtype; boundary_points = xb_structured, fields = sim.boundary_data.fields,outward_normals = sim.boundary_data.outward_normals, interior_points = sim.boundary_data.interior_points)
            M_fw = system_matrix(structured_chi, sim.medium, bd_fw)
            
            # --- Backward Step ---
            if coord == 1
                xb_structured[s_idx] = SVector(orig_svec[1] - h, orig_svec[2])
            else
                xb_structured[s_idx] = SVector(orig_svec[1], orig_svec[2] - h)
            end
            
            bd_bw = BoundaryData(sim.boundary_data.fieldtype; boundary_points = xb_structured, fields = sim.boundary_data.fields,outward_normals = sim.boundary_data.outward_normals, interior_points = sim.boundary_data.interior_points)
            M_bw = system_matrix(structured_chi, sim.medium, bd_bw)
            
            # --- Reset ---
            # Put the original SVector back before moving to the next coordinate
            xb_structured[s_idx] = orig_svec
            
            # Apply math safely (using the corrected explicit broadcast)
            @views jac_M[:, v] .= (vec(M_fw) .- vec(M_bw)) ./ (2 * h)
        end
        
        Cx = zeros(T, N, N)
        is_full_sigma = (Sigma_x_input isa AbstractMatrix) && size(Sigma_x_input, 1) == n_xb
        
        for r in 1:N, s in 1:N
            i = div(r - 1, d_m) + 1
            j = div(s - 1, d_m) + 1
            
            if i == j 
                Sigma_x_block = is_full_sigma ? Sigma_x_input[(2i-1):(2i), (2i-1):(2i)] : Sigma_x_input
                J_r = zeros(T, K, 2)
                J_s = zeros(T, K, 2)
                
                for c in 1:K
                    idx_r = r + (c - 1) * N
                    idx_s = s + (c - 1) * N
                    
                    J_r[c, 1] = jac_M[idx_r, 2i - 1]; J_r[c, 2] = jac_M[idx_r, 2i]
                    J_s[c, 1] = jac_M[idx_s, 2j - 1]; J_s[c, 2] = jac_M[idx_s, 2j]
                end
                
                Cx[r, s] = tr(Sigma_a * J_r * Sigma_x_block * J_s')
            end
        end
        return Cx
    end
end
#Objective function: Compute the negative log-marginal likelihood.

function log_marginal_likelihood(chi::AbstractVector, sim)
    g, _, _, Sigma_a, Sigma_sensor, _ = extract_bayesian_components(sim)

    structured_chi=Vector(reinterpret(SVector{2, eltype(chi)}, chi))

    M_nom = system_matrix(structured_chi, sim.medium, sim.boundary_data)
    Cx = geometric_covariance(sim, structured_chi)
    
    Sigma_gg = Sigma_sensor + Cx + M_nom * Sigma_a * M_nom'
    return 0.5 * (g' * (Sigma_gg \ g) + logdet(Sigma_gg) + length(g) * log(2π))
end

function log_marginal_likelihood(chi::AbstractVector, diagSigma_a:: AbstractVector, sim)
    g, _, _, _, Sigma_sensor, _ = extract_bayesian_components(sim)

    Sigma_a=Diagonal(diagSigma_a)

    structured_chi=Vector(reinterpret(SVector{2, eltype(chi)}, chi))

    M_nom = system_matrix(structured_chi, sim.medium, sim.boundary_data)
    Cx = geometric_covariance(sim, structured_chi)
    
    Sigma_gg = Sigma_sensor + Cx + M_nom * Sigma_a * M_nom'
    return 0.5 * (g' * (Sigma_gg \ g) + logdet(Sigma_gg) + length(g) * log(2π))
end

function optimise_source_positions(sim, init_chi::AbstractVector)
    obj(chi) = log_marginal_likelihood(chi, sim)
    res = optimize(obj, Vector(init_chi), LBFGS())
    return Optim.minimizer(res)
end

function optimise_source_positions(sim)
    optimise_source_positions(sim, vcat(sim.source_positions...))
end

function compute_coefficient_posterior(sim, chi::AbstractVector)
    g, _, _, Sigma_a, Sigma_sensor, _ = extract_bayesian_components(sim)

    structured_chi=Vector(reinterpret(SVector{2, eltype(chi)}, Vector(chi)))

    M_nom = system_matrix(structured_chi, sim.medium, sim.boundary_data)
    
    Cx = geometric_covariance(sim, structured_chi)
    
    R_eff = Sigma_sensor + Cx
    R_inv = inv(R_eff)
    
    Sigma_post = inv(M_nom' * R_inv * M_nom + inv(Sigma_a))
    mu_post = Sigma_post * (M_nom' * R_inv * g)
    
    return mu_post, Sigma_post
end



function construct_prior(
    sim
)
    # 1. Extract all current components from the simulation
    g, xb_flat, init_chi, Σ_a_init, Σ_sensor, Σ_x = extract_bayesian_components(sim)
    

    # 2. Extract initial diagonal variances from the MvNormal prior
    # We add a tiny nugget (1e-12) to prevent log(0) if the initial variance is strictly zero
    init_variances = diag(Σ_a_init)
    init_omega = log.(init_variances .+ 1e-12) 
    
    # 3. Concatenate into a single master parameter vector
    theta_init = [init_chi; init_omega]
    n_chi = length(init_chi)
    
    # 4. Define the joint objective function
    function joint_obj(theta)
        # Unpack the master vector
        chi_current = theta[1:n_chi]
        omega_current = theta[n_chi+1:end]
        
        # Reconstruct the strictly positive diagonal covariance matrix
        Sigma_a_current = Diagonal(exp.(omega_current))
        
        # Call the appropriate typed log_marginal_likelihood
        
        return log_marginal_likelihood(chi_current, omega_current, sim)
        
    end
    
    # 6. Run the joint optimization
    # Note: If you want to use ForwardDiff, change to LBFGS(), autodiff=AutoForwardDiff()
    res = optimize(joint_obj, Vector(theta_init), LBFGS())
    theta_opt = Optim.minimizer(res)
    
    # 7. Unpack the optimized results
    chi_opt = theta_opt[1:n_chi]

    chi_opt_matrix=reshape(chi_opt, 2, :)

    chi_opt=[chi_opt_matrix[:,i] for i in 1:size(chi_opt_matrix,2)]

    omega_opt = theta_opt[n_chi+1:end]
    
    # 8. Reconstruct the optimized MvNormal prior
    opt_variances = exp.(omega_opt)
    opt_prior = MvNormal(zeros(length(opt_variances)), Diagonal(opt_variances))
    
    # 9. Create the updated BayesianSolver
    opt_solver = BayesianSolver(
        opt_prior; 
        optimise_source_positions_flag = sim.solver.optimise_source_positions_flag,
        use_greens_gradient_analytical_flag = sim.solver.use_greens_gradient_analytical_flag
    )
    
    opt_sim = Simulation(
    sim.medium, 
    sim.boundary_data; 
    solver =opt_solver,
    source_positions = chi_opt,
    particular_solution = NoParticularSolution(),
    ω = 2pi * 1.0
    )


    return opt_sim
end




