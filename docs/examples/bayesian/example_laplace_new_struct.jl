using MethodOfFundamentalSolutions
using Plots
using Distributions
using LinearAlgebra

# Assuming your module exports or includes these types:
# using MethodOfFundamentalSolutions 

# ==============================================================================
# 1. Geometry Setup & Spatial Uncertainty
# ==============================================================================
n_sensors = 100
θs = range(0, 2π, length=n_sensors+1)[1:end-1]

# Create structured SVectors for the true sensor positions (unit circle)
x0_sensors_true = [[cos(θ), sin(θ)] for θ in θs]

# Noise properties for sensor locations
σ_x = 0.01
Σ_x_block = σ_x^2 * I(2)
Σ_x = kron(I(n_sensors), Σ_x_block)

# Add noise to sensor positions to generate our experimental realization
noise_x_distribution = MvNormal(zeros(2 * n_sensors), Σ_x)
noise_x_flat = rand(noise_x_distribution)
noise_x_structured = [noise_x_flat[i:i+1] for i in 1:2:length(noise_x_flat)]
x0_sensors_noisy = x0_sensors_true .+ noise_x_structured

# Create the spatial distribution representing our knowledge of sensor positions
# (Centered at the noisy observed positions, with the known experimental covariance)
flat_noisy_sensors = vcat(x0_sensors_noisy...)
sensor_distribution = MvNormal(flat_noisy_sensors, Σ_x)

# ==============================================================================
# 2. Boundary Data & Measurement Uncertainty
# ==============================================================================
# Analytical boundary condition: u(x,y) = x
g_true = [p[1] for p in x0_sensors_true] 

σ_sensor = 0.05
Σ_sensor = σ_sensor^2 * I(n_sensors)

# Generate noisy measurements from the true field values
measurement_generator = MvNormal(g_true, Σ_sensor)
g_noisy = rand(measurement_generator)

# Create the field distribution representing our observed data + uncertainty
field_distribution = MvNormal(g_noisy, Σ_sensor)

# Create boundary data from the sensor distribution and field distribution
bd = BoundaryData(
    DirichletType();
    boundary_points = sensor_distribution,  # Literal AbstractMvNormal passed here
    fields = field_distribution              # Literal AbstractMvNormal passed here
)

# ==============================================================================
# 3. Solver Setup (Prior & Flags)
# ==============================================================================
n_sources = 10
σ_a = 1.0
Σ_a = σ_a^2 * I(n_sources)

# Create the Gaussian prior for the coefficients using standard Distributions.jl
prior_mean = zeros(n_sources) 
prior = MvNormal(prior_mean, Σ_a)

# Instantiate your Bayesian Solver configuration
solver = BayesianSolver(
    prior;
    optimise_source_positions_flag = true,
    use_greens_gradient_analytical_flag = true
)

 # Place initial guess for sources in a larger circle (r = 2.0)

θ_sources = range(0, 2π, length=n_sources+1)[1:end-1]

init_source_positions = [[2.0*cos(θ), 2.0*sin(θ)] for θ in θ_sources] 

# ==============================================================================
# 4. Simulation Assembly & Solving
# ==============================================================================
# Define the physical medium (Assuming standard acoustic/laplace dummy medium)
medium = LaplaceMedium{2, Float64}()

# Package into Simulation
sim = Simulation(
    medium, 
    bd; 
    solver =solver,
    source_positions = init_source_positions,
    particular_solution = NoParticularSolution(),
    ω = 2pi * 1.0
)

sol = solve(sim)

# ==============================================================================
# 5. Field Prediction 
# ==============================================================================
# Get points inside the unit disk for field prediction
using MultipleScattering
grid, idx = points_in_shape(Circle([0.0, 0.0], 1.0))
points = grid[idx]

fs = [
    field(DirichletType(), sol, x, x / norm(x)) 
for x in points];
    
field_mat = [[0.0] for x in grid]
field_mat[idx] = [[fs[i][1]] for i in eachindex(fs)];

field_predict = FieldResult(grid, [field_mat[i] for i in eachindex(field_mat)]);


p1 = plot(field_predict, field_apply = first, title = "Predicted Field")
plot(sol)
covs= [
    field_covariance(DirichletType(), sol, x, x / norm(x)) 
for x in points]

stds = [
    field_std(DirichletType(), sol, x, x / norm(x))
        for x in points]

std_mat = [[0.0] for x in grid]
std_mat[idx] = [[stds[i][1]] for i in eachindex(stds)];

std_predict = FieldResult(grid, [std_mat[i] for i in eachindex(std_mat)]);

p2 = plot(std_predict, field_apply = first, title = "Standard Deviation", colormap = :inferno)



plot(p1, p2, layout = (1, 2), size = (800, 400))

plot!(sol)


# ==============================================================================
# 6. If you want to optimise the prior and source positions simulataneously, you can call:
# ==============================================================================
solver_2 = BayesianSolver(
    prior;
    optimise_source_positions_flag = false,
    use_greens_gradient_analytical_flag = true
)

sim_2 = Simulation(
    medium, 
    bd; 
    solver =solver_2,
    source_positions = init_source_positions,
    particular_solution = NoParticularSolution(),
    ω = 2pi * 1.0
)


opt_sim= construct_prior(sim_2)

opt_sol = solve(opt_sim)


# ==============================================================================
# 5. Field Prediction 
# ==============================================================================
# Get points inside the unit disk for field prediction
using MultipleScattering
grid, idx = points_in_shape(Circle([0.0, 0.0], 1.0))
points = grid[idx]

fs = [
    field(DirichletType(), opt_sol, x, x / norm(x)) 
for x in points];
    
field_mat = [[0.0] for x in grid]
field_mat[idx] = [[fs[i][1]] for i in eachindex(fs)];

field_predict = FieldResult(grid, [field_mat[i] for i in eachindex(field_mat)]);


p1 = plot(field_predict, field_apply = first, title = "Predicted Field")
plot!(opt_sol)
covs= [
    field_covariance(DirichletType(), opt_sol, x, x / norm(x)) 
for x in points]

stds = [
    field_std(DirichletType(), opt_sol, x, x / norm(x))
        for x in points]

std_mat = [[0.0] for x in grid]
std_mat[idx] = [[stds[i][1]] for i in eachindex(stds)];

std_predict = FieldResult(grid, [std_mat[i] for i in eachindex(std_mat)]);

p2 = plot(std_predict, field_apply = first, title = "Standard Deviation", colormap = :inferno)



plot(p1, p2, layout = (1, 2), size = (800, 400))

plot(opt_sol)
