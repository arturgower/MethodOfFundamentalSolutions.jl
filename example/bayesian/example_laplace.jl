using MultipleScattering
using MethodOfFundamentalSolutions
using Plots
using Distributions
using StaticArrays:SVector
using Accessors
using LinearAlgebra
using ForwardDiff
#include("../../src/bayesian.jl")
#include("../../src/physics/laplace.jl")
#include("../../plot/plot.jl")
#include("../../src/utils.jl")

# Example of Laplace equation with MFS and Bayesian inference, considering a unit disk. and boundary condition u=x in the boundary. Which gives the solution u(x,y)=x in the whole domain in the analytical case.

# Sensor postions and noise in the postions


θs = range(0, 2π, length=101)[1:end-1]
#θs = range(0, 2π, length=4)[1:end-1]
x0_sensors = [cos.(θs)'; sin.(θs)'] # 100 sensors in a circle

x0_sensors_flat=x0_sensors[:]

n_sensors=length(θs)

σ_x=0.01
#σ_x=0.00

 Σ_x_block = σ_x^2 * I(2) # 1 cm std dev in x and y, no correlation

Σ_x=kron(I(n_sensors), Σ_x_block) # Block diagonal covariance for all sensors

noise_x_distribution = MvNormal(zeros(length(x0_sensors_flat)), Σ_x)

noise_x = rand(noise_x_distribution)


x0_sensors_noisy_flat = x0_sensors_flat + noise_x

## boundary data and measurements noise
g_true= x0_sensors[1,:] # u=x in the boundary
σ=0.05
Σ_sensor= σ^2*I(n_sensors)
measurement_distribution = MvNormal(g_true, Σ_sensor)

g = rand(measurement_distribution)

## Initial guess for the source positions (hyperparameters)
n_sources=10

σ_a=1.0

Σ_a=σ_a^2*I(n_sources)

 init_source_positions_true = [2.0*cos.(range(0, 2π, length=n_sources+1)[1:end-1])'; 2.0*sin.(range(0, 2π, length=n_sources+1)[1:end-1])'] #10 sources in a larger circle 
 init_source_positions_flat=init_source_positions_true[:]
 σ_source=0.00
 Σ_source_block=σ_source^2*I(2) # 20 cm std dev in x and y, no correlation
 Σ_source=kron(I(n_sources), Σ_source_block) # Block diagonal covariance for all sources

 source_distribution = MvNormal(init_source_positions_flat, Σ_source)

 init_source_positions= rand(source_distribution)

compute_Cx(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, laplace_M)
compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, laplace_grad_M)

 norm(compute_Cx(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, laplace_M) - compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, laplace_grad_M))/norm(compute_Cx(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, laplace_M)) # Relative error between finite difference and analytical Cx, which should be small if the gradient is correct.

 ## Optimize the hyperparameters (source positions) using the analytical gradient of the log-marginal likelihood.

best_source_positions = optimize_hyperparameters(
    g, 
    x0_sensors_noisy_flat, 
    Σ_a, 
    Σ_sensor, 
    Σ_x_block, 
    init_source_positions, 
    laplace_M, 
    laplace_grad_M
)



best_source_positions_matrix = flat_to_pos_matrix(best_source_positions)


x_sources, y_sources = best_source_positions_matrix[1, :], best_source_positions_matrix[2, :]

#init_x_sources, init_y_sources = 
scatter(x_sources, y_sources, 
        label="Sources (x)", 
        color=:red, 
        markersize=6, 
        marker=:circle, 
        aspect_ratio=:equal, 
        grid=true) 

# 2. Calculate parametric coordinates for r = 1

plot!(Circle([0.0, 0.0], 1.0) )


#Compute posterior predictive distribution of the measurements at the sensor locations, given the optimized source positions.


μ_post, Σ_post = compute_coefficient_posterior(g, x0_sensors_noisy_flat, best_source_positions, Σ_a, Σ_sensor, Σ_x_block, laplace_M, laplace_grad_M)
#μ_post, Σ_post = compute_coefficient_posterior(g, x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_sensor, Σ_x_block, laplace_M, laplace_grad_M)


#Compute field predictions in the domain

grid,idx=points_in_shape(Circle([0.0, 0.0], 1.0)) |> collect

points=grid[idx] |> collect
points=Vector.(points)
points_matrix= hcat(points...)  |> collect

points_flat=points_matrix[:]

mean_field, variance_field=reconstruct_full_field(points_flat, best_source_positions, μ_post, Σ_post, laplace_M) |> collect
#mean_field, variance_field=reconstruct_full_field(points_flat, init_source_positions, μ_post, Σ_post, laplace_M) |> collect

fs=zeros(length(grid))
fs[idx]=mean_field


plot_reconstructed_fields_mean(points_flat, mean_field, x0_sensors_flat, source_positions_flat=best_source_positions)
#savefig("images/laplace_reconstructed_mean_optimized.pdf")
plot_reconstructed_fields_mean(points_flat, mean_field, x0_sensors_flat)
plot_reconstructed_fields_variance(points_flat, variance_field, x0_sensors_flat, 
source_positions_flat=best_source_positions)
#savefig("images/laplace_reconstructed_variance_optimized.pdf")

#plot_reconstructed_fields_mean(points_flat, mean_field, init_source_positions, x0_sensors_flat)
#plot_reconstructed_fields_variance(points_flat, variance_field, init_source_positions, x0_sensors_flat)

plot_reconstructed_fields_mean(points_flat, mean_field, x0_sensors_flat, source_positions_flat=best_source_positions)
plot_reconstructed_fields(points_flat, mean_field, variance_field, x0_sensors_flat, source_positions_flat=best_source_positions)
#savefig("images/laplace_reconstructed_fields_optimized.pdf")
#plot_reconstructed_fields(points_flat, mean_field, variance_field, x0_sensors_flat, source_positions_flat=init_source_positions)
#savefig("images/laplace_reconstructed_fields_initial.pdf")
norm(mean_field - points_matrix[1,:])/norm(points_matrix[1,:]) # Relative error of the mean since the true solution is u(x,y)=x, which corresponds to points_matrix[1,:] in the grid.