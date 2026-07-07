# ==============================================================================
# TEST SET 1: Consistency (Analytical vs Finite Difference)
# ==============================================================================
@testset "geometric_covariance_consistency_test_laplace" begin
    # ==============================================================================
    # 1. Geometry Setup & Spatial Uncertainty
    # ==============================================================================
    n_sensors = 100
    θs = range(0, 2π, length=n_sensors+1)[1:end-1]
    x0_sensors_true = [[cos(θ), sin(θ)] for θ in θs]

    σ_x = 0.01
    Σ_x_block = σ_x^2 * I(2)
    Σ_x = kron(I(n_sensors), Σ_x_block)

    noise_x_distribution = MvNormal(zeros(2 * n_sensors), Σ_x)
    noise_x_flat = rand(noise_x_distribution)
    noise_x_structured = [noise_x_flat[i:i+1] for i in 1:2:length(noise_x_flat)]
    x0_sensors_noisy = x0_sensors_true .+ noise_x_structured

    flat_noisy_sensors = vcat(x0_sensors_noisy...)
    sensor_distribution = MvNormal(flat_noisy_sensors, Σ_x)

    # ==============================================================================
    # 2. Boundary Data & Measurement Uncertainty
    # ==============================================================================
    g_true = [p[1] for p in x0_sensors_true] 

    σ_sensor = 0.05
    Σ_sensor = σ_sensor^2 * I(n_sensors)

    measurement_generator = MvNormal(g_true, Σ_sensor)
    g_noisy = rand(measurement_generator)
    field_distribution = MvNormal(g_noisy, Σ_sensor)

    bd = BoundaryData(
        DirichletType();
        boundary_points = sensor_distribution,  
        fields = field_distribution              
    )

    # ==============================================================================
    # 3. Solver & Simulation Setup
    # ==============================================================================
    n_sources = 10
    σ_a = 1.0
    Σ_a = σ_a^2 * I(n_sources)
    prior = MvNormal(zeros(n_sources), Σ_a)

    θ_sources = range(0, 2π, length=n_sources+1)[1:end-1]
    init_source_positions = [[2.0*cos(θ), 2.0*sin(θ)] for θ in θ_sources] 

    medium = LaplaceMedium{2, Float64}()

    # --- Simulation 1: Analytical ---
    solver_analytical = BayesianSolver(
        prior;
        optimise_source_positions_flag = true,
        use_greens_gradient_analytical_flag = true
    )

    sim = Simulation(
        medium, bd; 
        solver = solver_analytical,
        source_positions = init_source_positions,
        particular_solution = NoParticularSolution(),
        ω = 2pi * 1.0
    )

    # --- Simulation 2: Finite Difference ---
    solver_fd = BayesianSolver(
        prior;
        optimise_source_positions_flag = true,
        use_greens_gradient_analytical_flag = false
    )

    sim_fd = Simulation(
        medium, bd; 
        solver = solver_fd,
        source_positions = init_source_positions,
        particular_solution = NoParticularSolution(),
        ω = 2pi * 1.0
    )
    
    # Define source positions as SVectors
    structured_chi = [SVector{2, Float64}(pos[1], pos[2]) for pos in init_source_positions]
    
    # Compute Exact Analytical Matrix
    # (Assuming `sim` is set up with use_greens_gradient_analytical_flag = true)
    Cx_analytical = geometric_covariance(sim, structured_chi)
    
    # Compute Numerical Matrix
    # (Assuming `sim_fd` is set up with use_greens_gradient_analytical_flag = false)
    h_test = 1e-5
    Cx_fd = geometric_covariance(sim_fd, structured_chi; h = h_test)

    # Check 1: The two methods must match closely (Tolerance tuned for FD at h=1e-5)
    @test isapprox(Cx_analytical, Cx_fd, rtol=1e-10)
end

@testset "geometric_covariance_consistency_test_elastostatic" begin
    # 1. Medium and Nominal Geometry Setup
    medium = Elastostatic(2; ρ = 1.0, cp = 3.0, cs = 2.0)
    
    n = 12; L = 0.5; H = 0.5
    x = [
        LinRange(-L, L, n+2)[2:end-1]; zeros(n) .+ L; 
        LinRange(L, -L, n+2)[2:end-1]; zeros(n) .- L
    ]
    y = [
        zeros(n); LinRange(0, H, n+2)[2:end-1];
        zeros(n) .+ H; LinRange(H, 0, n+2)[2:end-1]
    ]
    points_true = [[x[i], y[i]] for i in eachindex(x)]
    n_sensors = length(points_true)

    # 2. Spatial Uncertainty (Geometric Noise for Sensors)
    σ_x = 0.01
    Σ_x_block = σ_x^2 * I(2)
    Σ_x = kron(I(n_sensors), Σ_x_block)

    noise_x_distribution = MvNormal(zeros(2 * n_sensors), Σ_x)
    noise_x_flat = rand(noise_x_distribution)
    noise_x_structured = [noise_x_flat[i:i+1] for i in 1:2:length(noise_x_flat)]
    points_noisy = points_true .+ noise_x_structured

    # Sensor positions are typically interleaved: [x1, y1, x2, y2, ...]
    flat_noisy_sensors = vcat(points_noisy...)
    sensor_distribution = MvNormal(flat_noisy_sensors, Σ_x)

    # 3. Boundary Data & Measurement Uncertainty
    weight = 2L * H * medium.ρ * 9.81
    traction_magnitude = weight / (2L)
    
    # Block Format for Vector Field (All X, then All Y)
    traction_x = zeros(n_sensors)
    traction_y = zeros(n_sensors)
    traction_y[1:n] .= traction_magnitude # Apply to the bottom wall only
    flat_traction_true = vcat(traction_x, traction_y)

    σ_traction = max(0.01 * maximum(abs.(flat_traction_true)), 1e-5) 
    Σ_traction = (σ_traction^2) * I(2 * n_sensors)
    
    measurement_generator = MvNormal(flat_traction_true, Σ_traction)
    g_noisy = rand(measurement_generator)
    field_distribution = MvNormal(g_noisy, Σ_traction)

    bd = BoundaryData(TractionType(); 
        boundary_points = sensor_distribution, 
        fields = field_distribution
    )

    # 4. Source Positions & Priors
    ns_per_side = 12
    source_positions = Vector{Float64}[]
    for x_val in LinRange(-2.0, 2.0, ns_per_side+1)[1:end-1]; push!(source_positions, [x_val, -1.75]); end
    for y_val in LinRange(-1.75, 2.25, ns_per_side+1)[1:end-1]; push!(source_positions, [2.0, y_val]); end
    for x_val in LinRange(2.0, -2.0, ns_per_side+1)[1:end-1]; push!(source_positions, [x_val, 2.25]); end
    for y_val in LinRange(2.25, -1.75, ns_per_side+1)[1:end-1]; push!(source_positions, [-2.0, y_val]); end
    n_sources = length(source_positions)

    σ_prior = 1.0
    Σ_prior = (σ_prior^2) * I(2 * n_sources) # 2 DoF per source for elastostatics
    prior_distribution = MvNormal(zeros(2 * n_sources), Σ_prior)

    # 5. Solver & Simulation Setup
    solver_analytical = BayesianSolver(
        prior_distribution;
        optimise_source_positions_flag = true, 
        use_greens_gradient_analytical_flag = true
    )

    solver_fd = BayesianSolver(
        prior_distribution;
        optimise_source_positions_flag = true, 
        use_greens_gradient_analytical_flag = false
    )

    sim_analytical = Simulation(medium, bd; 
        solver = solver_analytical,
        source_positions = source_positions,
        particular_solution = ParticularGravity(height = H)
    )

    sim_fd = Simulation(medium, bd; 
        solver = solver_fd,
        source_positions = source_positions,
        particular_solution = ParticularGravity(height = H)
    )

    # 6. Geometric Covariance Evaluation
    structured_chi = [SVector{2, Float64}(pos[1], pos[2]) for pos in source_positions]
    
    Cx_analytical = geometric_covariance(sim_analytical, structured_chi)
    
    h_test = 1e-5
    Cx_fd = geometric_covariance(sim_fd, structured_chi; h = h_test)

    # The matrices must match closely at h=1e-5
    @test isapprox(Cx_analytical, Cx_fd, rtol=1e-10)
end
# ==============================================================================
# TEST SET 2: Finite-Difference Convergence Rate
# ==============================================================================
@testset "geometric_covariance_convergence_test" begin
    # ==============================================================================
    # 1. Geometry Setup & Spatial Uncertainty
    # ==============================================================================
    n_sensors = 100
    θs = range(0, 2π, length=n_sensors+1)[1:end-1]
    x0_sensors_true = [[cos(θ), sin(θ)] for θ in θs]

    σ_x = 0.01
    Σ_x_block = σ_x^2 * I(2)
    Σ_x = kron(I(n_sensors), Σ_x_block)

    noise_x_distribution = MvNormal(zeros(2 * n_sensors), Σ_x)
    noise_x_flat = rand(noise_x_distribution)
    noise_x_structured = [noise_x_flat[i:i+1] for i in 1:2:length(noise_x_flat)]
    x0_sensors_noisy = x0_sensors_true .+ noise_x_structured

    flat_noisy_sensors = vcat(x0_sensors_noisy...)
    sensor_distribution = MvNormal(flat_noisy_sensors, Σ_x)

    # ==============================================================================
    # 2. Boundary Data & Measurement Uncertainty
    # ==============================================================================
    g_true = [p[1] for p in x0_sensors_true] 

    σ_sensor = 0.05
    Σ_sensor = σ_sensor^2 * I(n_sensors)

    measurement_generator = MvNormal(g_true, Σ_sensor)
    g_noisy = rand(measurement_generator)
    field_distribution = MvNormal(g_noisy, Σ_sensor)

    bd = BoundaryData(
    DirichletType();
    boundary_points = sensor_distribution,  
    fields = field_distribution              
    )

    # ==============================================================================
    # 3. Solver & Simulation Setup
    # ==============================================================================
    n_sources = 10
    σ_a = 1.0
    Σ_a = σ_a^2 * I(n_sources)
    prior = MvNormal(zeros(n_sources), Σ_a)

    θ_sources = range(0, 2π, length=n_sources+1)[1:end-1]
    init_source_positions = [[2.0*cos(θ), 2.0*sin(θ)] for θ in θ_sources] 

    medium = LaplaceMedium{2, Float64}()

    # --- Simulation 1: Analytical ---
    solver_analytical = BayesianSolver(
    prior;
    optimise_source_positions_flag = true,
    use_greens_gradient_analytical_flag = true
    )

    sim = Simulation(
    medium, bd; 
    solver = solver_analytical,
    source_positions = init_source_positions,
    particular_solution = NoParticularSolution(),
    ω = 2pi * 1.0
    )

    # --- Simulation 2: Finite Difference ---
    solver_fd = BayesianSolver(
    prior;
    optimise_source_positions_flag = true,
    use_greens_gradient_analytical_flag = false
    )

    sim_fd = Simulation(
    medium, bd; 
    solver = solver_fd,
    source_positions = init_source_positions,
    particular_solution = NoParticularSolution(),
    ω = 2pi * 1.0
    )
    structured_chi = [SVector{2, Float64}(pos[1], pos[2]) for pos in init_source_positions]
    
    # Baseline analytical computation
    Cx_analytical = geometric_covariance(sim, structured_chi)
    norm_analytical = norm(Cx_analytical)

    # Sequence of step sizes
    hs = [1e-2, 1e-3, 1e-4, 1e-5]
    relative_errors = Float64[]

    for h in hs
        # Generate the numerical covariance matrix for the current step size
        Cx_num = geometric_covariance(sim_fd, structured_chi; h = h)

        # Track the relative error compared to the exact analytical evaluation
        err = norm(Cx_num - Cx_analytical) / norm_analytical
        push!(relative_errors, err)
    end

    # Calculate empirical convergence rates: Δlog(Error) / Δlog(h)
    rates = [log(relative_errors[i+1] / relative_errors[i]) / log(hs[i+1] / hs[i]) for i in 1:(length(hs)-1)]

    # Assertions
    @test all(diff(relative_errors) .< 0) # Error must drop strictly as h decreases

    expected_order = 2.0 # Central finite difference is O(h^2)
    
    for (j, rate) in enumerate(rates)
        @info "Testing tensor step h = $(hs[j]) -> $(hs[j+1]) | Observed Convergence Rate: $(round(rate, digits=4))"
        
        # Test that the convergence rate holds near 2.0
        @test isapprox(rate, expected_order, atol=0.15)
    end
end

@testset "geometric_covariance_convergence_test_elastostatic" begin
    # 1. Geometry Setup (Reusing smaller variables for speed)
    medium = Elastostatic(2; ρ = 1.0, cp = 3.0, cs = 2.0)
    
    n = 10; L = 0.5; H = 0.5
    x = [LinRange(-L, L, n+2)[2:end-1]; zeros(n) .+ L; LinRange(L, -L, n+2)[2:end-1]; zeros(n) .- L]
    y = [zeros(n); LinRange(0, H, n+2)[2:end-1]; zeros(n) .+ H; LinRange(H, 0, n+2)[2:end-1]]
    points_true = [[x[i], y[i]] for i in eachindex(x)]
    n_sensors = length(points_true)

    # 2. Geometric Uncertainty
    σ_x = 0.01
    Σ_x_block = σ_x^2 * I(2)
    Σ_x = kron(I(n_sensors), Σ_x_block)
    noise_x_distribution = MvNormal(zeros(2 * n_sensors), Σ_x)
    points_noisy = points_true .+ [rand(noise_x_distribution)[i:i+1] for i in 1:2:(2*n_sensors)]
    sensor_distribution = MvNormal(vcat(points_noisy...), Σ_x)

    # 3. Traction Uncertainty (Block Format)
    traction_x = zeros(n_sensors); traction_y = zeros(n_sensors)
    traction_y[1:n] .= (2L * H * medium.ρ * 9.81) / (2L) 
    flat_traction_true = vcat(traction_x, traction_y)

    Σ_traction = (max(0.01 * maximum(abs.(flat_traction_true)), 1e-5)^2) * I(2 * n_sensors)
    field_distribution = MvNormal(rand(MvNormal(flat_traction_true, Σ_traction)), Σ_traction)

    bd = BoundaryData(TractionType(); boundary_points = sensor_distribution, fields = field_distribution)

    # 4. Sources and Prior
    ns_per_side = 10
    source_positions = Vector{Float64}[]
    for x_val in LinRange(-2.0, 2.0, ns_per_side+1)[1:end-1]; push!(source_positions, [x_val, -1.75]); end
    for y_val in LinRange(-1.75, 2.25, ns_per_side+1)[1:end-1]; push!(source_positions, [2.0, y_val]); end
    for x_val in LinRange(2.0, -2.0, ns_per_side+1)[1:end-1]; push!(source_positions, [x_val, 2.25]); end
    for y_val in LinRange(2.25, -1.75, ns_per_side+1)[1:end-1]; push!(source_positions, [-2.0, y_val]); end

    prior = MvNormal(zeros(2 * length(source_positions)), (1.0^2) * I(2 * length(source_positions)))

    # 5. Simulations Setup
    sim_analytical = Simulation(medium, bd; 
        solver = BayesianSolver(prior; optimise_source_positions_flag=true, use_greens_gradient_analytical_flag=true),
        source_positions = source_positions, particular_solution = ParticularGravity(height = H)
    )

    sim_fd = Simulation(medium, bd; 
        solver = BayesianSolver(prior; optimise_source_positions_flag=true, use_greens_gradient_analytical_flag=false),
        source_positions = source_positions, particular_solution = ParticularGravity(height = H)
    )

    # 6. Evaluation Loop
    structured_chi = [SVector{2, Float64}(pos[1], pos[2]) for pos in source_positions]
    
    Cx_analytical = geometric_covariance(sim_analytical, structured_chi)
    norm_analytical = norm(Cx_analytical)

    hs = [1e-2, 1e-3, 1e-4, 1e-5]
    relative_errors = Float64[]

    for h in hs
        Cx_num = geometric_covariance(sim_fd, structured_chi; h = h)
        err = norm(Cx_num - Cx_analytical) / norm_analytical
        push!(relative_errors, err)
    end

    rates = [log(relative_errors[i+1] / relative_errors[i]) / log(hs[i+1] / hs[i]) for i in 1:(length(hs)-1)]

    # Error must drop strictly as h decreases
    @test all(diff(relative_errors) .< 0) 

    expected_order = 2.0 # Central finite difference is O(h^2)
    
    for (j, rate) in enumerate(rates)
        @info "Testing Elastostatic tensor step h = $(hs[j]) -> $(hs[j+1]) | Observed Convergence Rate: $(round(rate, digits=4))"
        @test isapprox(rate, expected_order, atol=0.15)
    end
end

# ==============================================================================
# TEST SET 3: Field Prediction Validation (Analytical vs Bayesian)
# ==============================================================================



@testset "Full Grid Bayesian Validation Square Gravity without position optimisation" begin
    # ==============================================================================
    # 1. Medium and Physical Domain Setup (Interleaved Layout)
    # ==============================================================================
    medium = Elastostatic(2; ρ = 1.0, cp = 3.0, cs = 2.0)

    n = 20; L = 0.5; H = 0.5
    x = [
        LinRange(-L, L, n+2)[2:end-1]; 
        zeros(n) .+ L; 
        LinRange(L, -L, n+2)[2:end-1]; 
        zeros(n) .- L
    ]
    y = [
        zeros(n);
        LinRange(0, H, n+2)[2:end-1];
        zeros(n) .+ H;
        LinRange(H, 0, n+2)[2:end-1]
    ]
    points = [[x[i], y[i]] for i in eachindex(x)]

    # Create normals (optional for traction type here, but kept for completeness)
    normals = [
        [ [0.0, -1.0] for i = 1:n]; 
        [ [1.0,  0.0] for i = 1:n]; 
        [ [0.0,  1.0] for i = 1:n]; 
        [ [-1.0, 0.0] for i = 1:n]
    ]

    weight = 2L * H * medium.ρ * 9.81
    traction_deterministic = [
        [ [0.0, weight / (2L)] for i = 1:n]; 
        [ [0.0, 0.0]           for i = 1:n]; 
        [ [0.0, 0.0]           for i = 1:n]; 
        [ [0.0, 0.0]           for i = 1:n]
    ]

    # Flatten using Interleaved Layout: [x1, y1, x2, y2, ...]
    flat_points = vcat(points...)
    flat_traction = vcat(vcat(traction_deterministic...)...)

    # ==============================================================================
    # 2. Bayesian Bounding Box: Large Source Square
    # ==============================================================================
    ns_per_side = 20
    source_positions = Vector{Float64}[]

    # Bottom Side
    for x_val in LinRange(-2.0, 2.0, ns_per_side+1)[1:end-1]
        push!(source_positions, [x_val, -1.75])
    end
    # Right Side
    for y_val in LinRange(-1.75, 2.25, ns_per_side+1)[1:end-1]
        push!(source_positions, [2.0, y_val])
    end
    # Top Side
    for x_val in LinRange(2.0, -2.0, ns_per_side+1)[1:end-1]
        push!(source_positions, [x_val, 2.25])
    end
    # Left Side
    for y_val in LinRange(2.25, -1.75, ns_per_side+1)[1:end-1]
        push!(source_positions, [-2.0, y_val])
    end

    n_sources = length(source_positions)

    # ==============================================================================
    # 3. Statistical Distributions & Solver Configuration
    # ==============================================================================
    # Define measurement uncertainty for the traction data
    σ_traction = 0.01 * maximum(abs.(flat_traction))
    Σ_traction = (σ_traction^2) * I(length(flat_traction))
    field_distribution = MvNormal(flat_traction, Σ_traction)

    # Define coordinate uncertainty for the boundary points
    σ_points = 0.01
    Σ_points = (σ_points^2) * I(length(flat_points))
    points_distribution = MvNormal(flat_points, Σ_points)

    # Define the prior distribution for the source charges
    σ_prior = 100.0
    Σ_prior = (σ_prior^2) * I(2 * n_sources)
    prior_distribution = MvNormal(zeros(2 * n_sources), Σ_prior)

    # Package into Bayesian Configurations
    bd = BoundaryData(TractionType(); 
        boundary_points = points_distribution, 
        fields = field_distribution
    )

    solver_bayesian = BayesianSolver(
        prior_distribution;
        optimise_source_positions_flag = false, 
        use_greens_gradient_analytical_flag = true
    )

    # ==============================================================================
    # 4. Simulation and Execution
    # ==============================================================================
    sim = Simulation(medium, bd; 
        particular_solution = ParticularGravity(height = H),
        solver = solver_bayesian,
        source_positions = source_positions
    )

    fsol = solve(sim)

    # ==============================================================================
    # 5. Field Prediction and Plotting
    # ==============================================================================
    # Set up a deterministic BoundaryData just to fetch the internal grid safely
    bd_plot = BoundaryData(TractionType(); 
        boundary_points = points, 
        fields = field_distribution 
    )

    grid, idx = points_in_shape(bd_plot)
    interior_points = grid[idx]
    normal_up = [0.0, 1.0]

    # Evaluate mean field and standard deviation
    fs = [field(TractionType(), fsol, p, normal_up) for p in interior_points]
    stds = [field_std(TractionType(), fsol, p, normal_up) for p in interior_points]
        
    # Extract the σ_yy component (index 2) and pad with NaN for clean plotting
    field_mat = [[NaN] for x in grid]
    field_mat[idx] = [[fs[i][2]] for i in eachindex(fs)]

    std_mat = [[NaN] for x in grid]
    std_mat[idx] = [[stds[i][2]] for i in eachindex(stds)]

    field_predict = FieldResult(grid, [field_mat[i] for i in eachindex(field_mat)])
    std_predict = FieldResult(grid, [std_mat[i] for i in eachindex(std_mat)])

    total_tested = length(interior_points)
    passed_tests = 0
    # ==============================================================================
    # 6. Bayesian Validation Test
    # ==============================================================================
    println("\n=== Running Full Grid Statistical Test ===")
    
    for (i, p) in enumerate(interior_points)
        y_coord = p[2]
        
        # Analytical vertical stress profile for a block under gravity
        analytical_σ_yy = medium.ρ * 9.81 * (y_coord - H)
        
        # Predicted stress and standard deviation from earlier calculations
        predicted_σ_yy = fs[i][2]
        predicted_std_σ_yy = stds[i][2]
        
        # 3-Sigma validation check
        abs_error = abs(predicted_σ_yy - analytical_σ_yy)
        three_sigma_bound = 3 * predicted_std_σ_yy
        
        @test abs_error <= three_sigma_bound
        
        if abs_error <= three_sigma_bound
            passed_tests += 1
        end
    end
    
    println("Test completed for the entire interior grid!")
    println("Successfully validated: $passed_tests / $total_tested interior points.")

    # ==============================================================================
    # 7. Boundary Residual Check (Self-Consistency)
    # ==============================================================================

    total_boundary_points = length(points)
    total_components = 2 * total_boundary_points # Checking both X and Y directions
    passed_boundary_tests = 0
    max_traction_error = 0.0
    
    println("\n=== Running Boundary Residual Check ===")
    
    for i in 1:total_boundary_points
        p = points[i]
        n_vec = normals[i]
        applied_traction = traction_deterministic[i]
        
        # Predict the traction vector [Tx, Ty] and its standard deviation at the boundary
        predicted_traction = field(TractionType(), fsol, p, n_vec)
        predicted_std = field_std(TractionType(), fsol, p, n_vec)
        
        # Check both X and Y components of the traction vector
        for dim in 1:2
            abs_error = abs(predicted_traction[dim] - applied_traction[dim])
            three_sigma_bound = 3 * predicted_std[dim]
            
            # Track the maximum error just to see how well it fits overall
            max_traction_error = max(max_traction_error, abs_error)
            
            @test abs_error <= three_sigma_bound
            
            if abs_error <= three_sigma_bound
                passed_boundary_tests += 1
            end
        end
    end
    
    println("Maximum boundary traction residual error: $max_traction_error")
    println("Successfully validated boundary components: $passed_boundary_tests / $total_components")


end
