# ==============================================================================
# 1. Helper Functions to Wrap `greens` for the Matrix/Tensor Builders
# ==============================================================================

# Builds the N x K Matrix using the scalar greens function
# Builds the N x K Matrix using the scalar greens function
# Builds the N x K Matrix using the scalar greens function
function build_M(xb_flat::AbstractVector, chi_flat::AbstractVector)
    N = length(xb_flat) ÷ 2
    K = length(chi_flat) ÷ 2
    
    T = promote_type(eltype(xb_flat), eltype(chi_flat))
    M = zeros(T, N, K)
    
    dummy_normal = SVector(0.0, 1.0)
    
    for i in 1:N
        x_i = SVector(xb_flat[2i-1], xb_flat[2i])
        for j in 1:K
            s_j = SVector(chi_flat[2j-1], chi_flat[2j])
            
            # ---> FIXED: changed sim.medium to medium <---
            val = greens(DirichletType(), medium, x_i - s_j, dummy_normal)
            M[i, j] = val[1]
        end
    end
    return M
end

# Builds the N x K x 2 Tensor using the greens_gradient function
function build_grad_M(xb_flat::AbstractVector, chi_flat::AbstractVector)
    N = length(xb_flat) ÷ 2
    K = length(chi_flat) ÷ 2
    
    T = promote_type(eltype(xb_flat), eltype(chi_flat))
    grad_M = zeros(T, N, K, 2)
    
    dummy_normal = SVector(0.0, 1.0)
    
    for i in 1:N
        x_i = SVector(xb_flat[2i-1], xb_flat[2i])
        for j in 1:K
            s_j = SVector(chi_flat[2j-1], chi_flat[2j])
            
            # ---> FIXED: changed sim.medium to medium <---
            gval = greens_gradient(DirichletType(), medium, x_i - s_j, dummy_normal)
            
            grad_M[i, j, 1] = gval[1]
            grad_M[i, j, 2] = gval[2]
        end
    end
    return grad_M
end# ==============================================================================
# 2. Shared Setup Data
# ==============================================================================
medium = LaplaceMedium{2, Float64}()
θs = range(0, 2π, length=101)[1:end-1]
n_sensors = length(θs)
x0_sensors = [cos.(θs)'; sin.(θs)'] 
x0_sensors_flat = x0_sensors[:]

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

# Note: Keeping source positions completely flat to match the new builder logic
init_source_positions = init_source_positions_flat 

# ==============================================================================
# TEST SET 1: Consistency (ForwardDiff vs Analytical)
# ==============================================================================
@testset "calculate_Cx_laplace_consistency_test" begin
    
    # 1. Compute Exact Matrices via your wrapper functions
    Cx_forwarddiff = compute_Cx(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, build_M)
    Cx_analytical  = compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, build_grad_M)

    # Check 1: The two exact methods must match down to near machine precision
    @test isapprox(Cx_forwarddiff, Cx_analytical, rtol=1e-12, atol=1e-12)
end

# ==============================================================================
# TEST SET 2: Finite-Difference Convergence of the Gradient Tensor
# ==============================================================================
@testset "calculate_Cx_laplace_convergence_test" begin
    
    Cx_analytical = compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, build_grad_M)
    norm_analytical = norm(Cx_analytical)

    hs = [1e-2, 1e-3, 1e-4, 1e-5]
    relative_errors = Float64[]

    # Determine dimensions based on the nominal evaluation
    M_nom = build_M(x0_sensors_noisy_flat, init_source_positions)
    N, K = size(M_nom)

    for h in hs
        grad_M_num = zeros(N, K, 2)
        
        for s_idx in 1:n_sensors
            idx_x = 2s_idx - 1
            idx_y = 2s_idx
            
            # Helper to compactly perturb the flat array
            perturb(flat, i, val) = (out = copy(flat); out[i] += val; out)
            
            # Perturb X Coordinate
            M_plus_x  = build_M(perturb(x0_sensors_noisy_flat, idx_x, h), init_source_positions)
            M_minus_x = build_M(perturb(x0_sensors_noisy_flat, idx_x, -h), init_source_positions)
            
            # Perturb Y Coordinate
            M_plus_y  = build_M(perturb(x0_sensors_noisy_flat, idx_y, h), init_source_positions)
            M_minus_y = build_M(perturb(x0_sensors_noisy_flat, idx_y, -h), init_source_positions)
            
            # Because moving sensor `s_idx` ONLY changes row `s_idx` in the M matrix,
            # we can optimize this by only updating that specific row!
            grad_M_num[s_idx, :, 1] .= (M_plus_x[s_idx, :] .- M_minus_x[s_idx, :]) ./ (2h)
            grad_M_num[s_idx, :, 2] .= (M_plus_y[s_idx, :] .- M_minus_y[s_idx, :]) ./ (2h)
        end

        # Feed the numerical tensor into your unchanged analytical layout function
        Cx_num = compute_Cx_analytical(x0_sensors_noisy_flat, init_source_positions, Σ_a, Σ_x_block, (xb, chi) -> grad_M_num)

        # Track the tracking error compared to your exact analytical evaluation
        err = norm(Cx_num - Cx_analytical) / norm_analytical
        push!(relative_errors, err)
    end

    # Calculate empirical convergence rates: Δlog(Error) / Δlog(h)
    rates = [log(relative_errors[i+1] / relative_errors[i]) / log(hs[i+1] / hs[i]) for i in 1:(length(hs)-1)]

    # Assertions
    @test all(diff(relative_errors) .< 0) # Error must drop strictly

    expected_order = 2.0 
    for (j, rate) in enumerate(rates)
        @info "Testing tensor step h = $(hs[j]) -> $(hs[j+1]) | Observed Convergence Rate: $(round(rate, digits=4))"
        @test isapprox(rate, expected_order, atol=0.15)
    end
end

