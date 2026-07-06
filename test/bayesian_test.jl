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
    @test isapprox(Cx_analytical, Cx_fd, rtol=1e-4)
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
    @test isapprox(Cx_analytical, Cx_fd, rtol=1e-4)
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