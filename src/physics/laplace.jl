

abstract type FieldType end

struct DirichletType <: FieldType end

struct NeumannType <: FieldType end

function greens(field::DirichletType, x::SVector{2,T}, outward_normal::AbstractVector{T} = zeros(T,2)) where T
   
    G = zeros(Complex{T},1,1)
    G[1,1] = -1/(2π) * log(norm(x))
    return G
end

"""
Laplace M_func where hyperparameters (chi) are the 2D coordinates of the sources.
"""
function laplace_M(x_flat, chi)
    # CRUCIAL FIX: Automatically resolve whether to use Float64 or ForwardDiff.Dual
    NumType = promote_type(eltype(x_flat), eltype(chi))
    
    n_sensors = div(length(x_flat), 2)
    n_sources = div(length(chi), 2)
    
    # Preallocate using the dynamically resolved type
    M = zeros(NumType, n_sensors, n_sources)
    
    for i in 1:n_sensors
        # Use NumType instead of hardcoding Float64
        x_sens = SVector{2, NumType}(x_flat[2i-1], x_flat[2i])
        
        for j in 1:n_sources
            x_sour = SVector{2, NumType}(chi[2j-1], chi[2j])
            
            r_vec = x_sens - x_sour
            
            # This calls your greens signature, which handles NumType natively
            G_complex_mat = greens(DirichletType(), r_vec)
            
            M[i, j] = real(G_complex_mat[1, 1])
        end
    end
    return M
end

"""
Analytical gradient tensor builder for Laplace.
Computes ∂M/∂x_sensor, where source locations are read from hyperparameter chi.
Returns a 3D array of shape (n_sensors, n_sources, 2)
"""

function laplace_grad_M(x_flat, chi)
    # Apply the same type promotion logic here
    NumType = promote_type(eltype(x_flat), eltype(chi))
    
    n_sensors = div(length(x_flat), 2)
    n_sources = div(length(chi), 2)
    
    # Preallocate using NumType
    grad_M = zeros(NumType, n_sensors, n_sources, 2)
    
    for i in 1:n_sensors
        x_sens = SVector{2, NumType}(x_flat[2i-1], x_flat[2i])
        
        for j in 1:n_sources
            x_sour = SVector{2, NumType}(chi[2j-1], chi[2j])
            
            dx = x_sens[1] - x_sour[1]
            dy = x_sens[2] - x_sour[2]
            r2 = dx^2 + dy^2 + 1e-12 
            
            grad_M[i, j, 1] = -(1.0 / (2 * π)) * (dx / r2) 
            grad_M[i, j, 2] = -(1.0 / (2 * π)) * (dy / r2) 
        end
    end
    return grad_M
end