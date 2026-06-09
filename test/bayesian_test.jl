@testset "calculate_Cx_laplace_convergence_and_consistency_test" begin
    # 1. Setup physical geometry and distributions
    θs = range(0, 2π, length=101)[1:end-1]
    x0_sensors = [cos.(θs)'; sin.(θs)'] # 100 sensors in a circle
    x0_sensors_flat = x0_sensors[:]
    n_sensors = length(θs)

    σ_x = 0.01
    Σ_x_block = σ_x^2 * I(2) 
    Σ_x = kron(I(n_sensors), Σ_x_block) 

    noise_x_distribution = MvNormal(zeros(length(x0_sensors_flat)), Σ_x)
    noise_x = rand(noise_x_distribution)
    x0_sensors_noisy_flat = x0_sensors_flat + noise_x

    n_sources = 10
    σ_a = 1.0
    Σ_a = σ_a^2 * I(n_sources)

    init_source_positions_true = [2.0*cos.(range(0, 2π, length=n_sources+1)[1:end-1])'; 2.0*sin.(range(0, 2π, length=n_sources+1)[1:end-1])'] 
    init_source_positions_flat = init_source_positions_true[:]
    
    σ_source = 0.00
    Σ_source_block = σ_source^2 * I(2)
    Σ_source = kron(I(n_sources), Σ_source_block) 

    source_distribution = MvNormal(init_source_positions_flat, Σ_source)
    init_source_positions = rand(source_distribution)

    # Wrap the base M evaluation to fit (xb_flat, chi) -> M
    M_func = (xb, chi) -> laplace_M(xb, chi)
    # Wrap your analytical gradient function to fit (xb_flat, chi) -> grad_M tensor
    grad_M_func = (xb, chi) -> laplace_grad_M(xb, chi)

    # 2. Compute Exact Matrices via your two unchanged functions
    Cx_forwarddiff = compute_Cx(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, M_func)
    Cx_analytical  = compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, grad_M_func)

    # Check 1: The two exact methods must match down to near machine precision
    @test isapprox(Cx_forwarddiff, Cx_analytical, rtol=1e-12, atol=1e-12)

    # 3. Define a vector of shrinking finite-difference step sizes
    hs = [1e-2, 1e-3, 1e-4, 1e-5]
    relative_errors = Float64[]
    norm_analytical = norm(Cx_analytical)

    # Determine dimensions based on the nominal evaluation
    M_nom = M_func(x0_sensors_noisy_flat, init_source_positions)
    N, K = size(M_nom)

    # 4. Generate finite-difference tensor evaluations to test convergence
    for h in hs
        # Allocate a 3D tensor to hold the numerical central differences of M
        # tracking spatial changes per sensor: [N_measurements, K_weights, 2 (x and y)]
        grad_M_num = zeros(N, K, 2)
        
        for s_idx in 1:n_sensors
            # Indicies for x and y components of the current sensor in the flattened array
            idx_x = 2s_idx - 1
            idx_y = 2s_idx
            
            # --- Perturb X Coordinate ---
            x_plus_x = copy(x0_sensors_noisy_flat);  x_plus_x[idx_x] += h
            x_minus_x = copy(x0_sensors_noisy_flat); x_minus_x[idx_x] -= h
            M_plus_x = M_func(x_plus_x, init_source_positions)
            M_minus_x = M_func(x_minus_x, init_source_positions)
            
            # --- Perturb Y Coordinate ---
            x_plus_y = copy(x0_sensors_noisy_flat);  x_plus_y[idx_y] += h
            x_minus_y = copy(x0_sensors_noisy_flat); x_minus_y[idx_y] -= h
            M_plus_y = M_func(x_plus_y, init_source_positions)
            M_minus_y = M_func(x_minus_y, init_source_positions)
            
            # Central difference maps directly into the spatial slices of the tensor
            # for the specific measurement rows influenced by this sensor
            for r in 1:N
                grad_M_num[r, :, 1] .+= (M_plus_x[r, :] - M_minus_x[r, :]) / (2h)
                grad_M_num[r, :, 2] .+= (M_plus_y[r, :] - M_minus_y[r, :]) / (2h)
            end
        end

        # Feed the numerical tensor into your unchanged analytical layout function
        Cx_num = compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, (xb, chi) -> grad_M_num)

        # Track the tracking error compared to your exact analytical evaluation
        err = norm(Cx_num - Cx_analytical) / norm_analytical
        push!(relative_errors, err)
    end

    # 5. Calculate empirical convergence rates: Δlog(Error) / Δlog(h)
    rates = [log(relative_errors[i+1] / relative_errors[i]) / log(hs[i+1] / hs[i]) for i in 1:(length(hs)-1)]

    # 6. Assertions
    # A: Error must strictly drop as h gets smaller
    @test all(diff(relative_errors) .< 0)

    # B: Rate must be approximately 2.0 (quadratic convergence of central differences)
    expected_order = 2.0 
    for (j, rate) in enumerate(rates)
        @info "Testing tensor step h = $(hs[j]) -> $(hs[j+1]) | Observed Convergence Rate: $(round(rate, digits=4))"
        @test isapprox(rate, expected_order, atol=0.15)
    end
end

