
using MethodOfFundamentalSolutions
using Plots
using Distributions
using StaticArrays
using Accessors
using LinearAlgebra
using ForwardDiff

# include("../../src/bayesian.jl")
# include("../../src/physics/laplace.jl")
# include("../../plot/plot.jl")
# include("../../src/utils.jl")

# ==============================================================================
# 1. Geometry Setup & Noise Generation
# ==============================================================================
n_sensors = 100
θs = range(0, 2π, length=n_sensors+1)[1:end-1]

# Create structured SVectors for the true sensor positions (unit circle)
x0_sensors_true = [SVector{2, Float64}(cos(θ), sin(θ)) for θ in θs]

# Noise properties for sensor locations
σ_x = 0.01
Σ_x_block = σ_x^2 * I(2)
Σ_x = kron(I(n_sensors), Σ_x_block)

# Add noise to sensor positions
noise_x_distribution = MvNormal(zeros(2 * n_sensors), Σ_x)
noise_x_flat = rand(noise_x_distribution)
# Zero-cost reshape flat noise back to structured SVector array
noise_x_structured = reinterpret(SVector{2, Float64}, noise_x_flat) 

x0_sensors_noisy = x0_sensors_true .+ noise_x_structured

# ==============================================================================
# 2. Boundary Data & Measurement Noise
# ==============================================================================
# Analytical boundary condition: u(x,y) = x
g_true = [p[1] for p in x0_sensors_true] 

σ_sensor = 0.05
Σ_sensor = σ_sensor^2 * I(n_sensors)

measurement_distribution = MvNormal(g_true, Σ_sensor)
g_noisy = rand(measurement_distribution)
# Explicitly ensure the types are correct
# Ensure your points are specifically a Vector of SVector{2, Float64}
# This matches the FV <: AbstractVector constraint
pts = SVector{2, Float64}.(x0_sensors_noisy)
fds = [SVector{1, Float64}(g) for g in g_noisy] 

bd = BoundaryData(
    DirichletType(),
    boundary_points = pts,
    fields = fds,
    fields_covariance = Σ_sensor,
    boundary_points_covariance = Σ_x_block
)
# ==============================================================================
# 3. Source Positions (Hyperparameters) & Prior
# ==============================================================================
n_sources = 10
σ_a = 1.0
Σ_a = σ_a^2 * I(n_sources)

# Place initial guess for sources in a larger circle (r = 2.0)
θ_sources = range(0, 2π, length=n_sources+1)[1:end-1]
init_source_positions = [SVector{2, Float64}(2.0*cos(θ), 2.0*sin(θ)) for θ in θ_sources]

# Create the Gaussian prior for the coefficients
prior_mean = [SVector{1, Float64}(0.0) for _ in 1:n_sources]
prior = GaussianPrior(prior_mean, Σ_a)

# ==============================================================================
# 4. Simulation Assembly & Solving
# ==============================================================================
# Define the physical medium (Assuming standard acoustic/laplace dummy medium)
medium = MethodOfFundamentalSolutions.Acoustic(2; ρ=1.0, c=1.0) 

# Package into Simulation
sim = Simulation(
    medium, 
    bd; 
    solver = BayesianSolver(prior),
    source_positions = init_source_positions,
    particular_solution = NoParticularSolution(),
    ω = 2pi * 1.0
)

# Call the clean solver API we built!
println("Starting optimization and posterior computation...")
sol = solve(
    sim, 
    laplace_M; 
    gradient_system_matrix_function = laplace_grad_M, 
    optimise_source_positions_flag = true
)

# ==============================================================================
# 5. Field Prediction & Verification
# ==============================================================================
# Get points inside the unit disk for field prediction
using MultipleScattering
grid, idx = points_in_shape(Circle([0.0, 0.0], 1.0))
points = grid[idx]
points_flat = vcat(points...)
best_source_positions_flat = vcat(sol.positions...)

# Reconstruct the full field (using your existing utility function)
mean_field, variance_field = reconstruct_full_field(
    points_flat, 
    best_source_positions_flat, 
    sol.coefficients, 
    sol.coefficients_covariance, 
    (x, y) -> laplace_M(sim.medium, x, y) # <--- The fix
)

# Compute relative error against the true analytical solution u(x,y) = x
true_field = [p[1] for p in points]
rel_error = norm(mean_field - true_field) / norm(true_field)
println("Relative error of the mean field: ", round(rel_error * 100, digits=3), "%")

# ==============================================================================
# 6. Plotting
# ==============================================================================
# Extract x and y coordinates from the optimized source SVectors
x_sources = [p[1] for p in sol.positions]
y_sources = [p[2] for p in sol.positions]

# Plot optimized sources and domain boundary
scatter(x_sources, y_sources, 
        label="Optimized Sources", 
        color=:red, 
        markersize=6, 
        marker=:circle, 
        aspect_ratio=:equal, 
        grid=true) 

plot!(Circle([0.0, 0.0], 1.0), label="Domain Boundary")
# savefig("images/laplace_optimized_sources.pdf")

x0_sensors_flat = vcat(x0_sensors_noisy...)

# Call your existing plot functions
plot_reconstructed_fields_mean(
    points_flat, mean_field, x0_sensors_flat, 
    source_positions_flat=best_source_positions_flat
)

plot_reconstructed_fields_variance(
    points_flat, variance_field, x0_sensors_flat, 
    source_positions_flat=best_source_positions_flat
)

plot_reconstructed_fields(
    points_flat, mean_field, variance_field, x0_sensors_flat, 
    source_positions_flat=best_source_positions_flat
)