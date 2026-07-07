using MethodOfFundamentalSolutions
using LinearAlgebra
using Distributions


# ==============================================================================
# 1. Medium and Physical Domain Setup
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

# Flatten the deterministic traction vector for the statistical distribution
# Length will be: 4 sides * n points/side * 2 components = 8n
flat_points = vcat(vcat(points...)...)
flat_traction = vcat(vcat(traction_deterministic...)...)

# ==============================================================================
# 2. Bayesian Bounding Box: Large Source Square
# ==============================================================================
# The boundary box lives within x ∈ [-0.5, 0.5] and y ∈ [0.0, 0.5].
# We create a significantly larger source square spanning x ∈ [-2.0, 2.0], y ∈ [-1.75, 2.25]
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

σ_traction = 0.01*maximum(abs.(flat_traction))  # 1% of the maximum traction magnitude
Σ_traction = (σ_traction^2) * I(length(flat_traction))
field_distribution = MvNormal(flat_traction, Σ_traction)

# Define the prior distribution for the source charges
# Since it's Elastostatics (2D vector field), length = 2 * n_sources
σ_prior = 100.0
Σ_prior = (σ_prior^2) * I(2 * n_sources)
prior_distribution = MvNormal(zeros(2 * n_sources), Σ_prior)

σ_points=0.01
Σ_points = (σ_points^2) * I(length(flat_points))
points_distribution = MvNormal(flat_points, Σ_points)
points_flat = rand(points_distribution)
points_distribution = MvNormal(points_flat, Σ_points)
# Package everything cleanly into the Bayesian configurations

bd = BoundaryData(TractionType(); 
    boundary_points = points_distribution, 
    fields = field_distribution,  # Statistical distribution passed here
    #outward_normals = normals
)

solver_bayesian = BayesianSolver(
    prior_distribution;
    optimise_source_positions_flag = true, 
    use_greens_gradient_analytical_flag = true,
    max_iters = 5,
)

# ==============================================================================
# 4. Simulation and Execution
# ==============================================================================
sim = Simulation(medium, bd; 
    particular_solution = ParticularGravity(height = H),
    #particular_solution = NoParticularSolution(),
    solver = solver_bayesian,
    source_positions = source_positions
)

fsol = solve(sim)

# ==============================================================================
# 5. Field Prediction (Displacement)
# ==============================================================================
bd_plot = BoundaryData(TractionType(); 
    boundary_points = points, 
    fields = field_distribution,  # Statistical distribution passed here
    #outward_normals = normals
)


grid, idx = points_in_shape(bd_plot)
points = grid[idx]
# Define a dummy normal vector required for vector field evaluations
normal = [0.0, 1.0]

# Evaluate the displacement field at internal points (returns [u_x, u_y])
fs = [
    field(TractionType(), fsol, x, normal) 
    for x in points
]
    
# Extract the first component (u_x displacement) for visualization
field_mat = [[NaN] for x in grid]
field_mat[idx] = [[fs[i][2]] for i in eachindex(fs)]

field_predict = FieldResult(grid, [field_mat[i] for i in eachindex(field_mat)])

using Plots
# Plot predicted field
p1 = plot(field_predict, field_apply = first, title = "Predicted Field (σyy)",colormap = :inferno)
plot!(fsol)

# Calculate Bayesian field uncertainty configurations
covs = [
    field_covariance(TractionType(), fsol, x, normal) 
    for x in points
]

stds = [
    field_std(TractionType(), fsol, x, normal)
    for x in points
]

# Extract the standard deviation for the u_x component
std_mat = [[NaN] for x in grid]
std_mat[idx] = [[stds[i][2]] for i in eachindex(stds)]

std_predict = FieldResult(grid, [std_mat[i] for i in eachindex(std_mat)])

# Plot standard deviation field
p2 = plot(std_predict, field_apply = first, title = "Standard Deviation (σyy)", colormap = :inferno)

# Render side-by-side plots with source layout overlay
plot(fsol)
plot(p1, p2, layout = (1, 2), size = (800, 400))
plot(fsol, xlims = (-1.0, 1.0), ylims = (-0.5, 0.5), title = "Source Layout Overlay")


