struct DirichletType <: FieldType end

struct NeumannType <: FieldType end

struct LaplaceMedium{Dim,T} <: PhysicalMedium{Dim, 1} end

function greens(
    field::DirichletType, 
    medium::LaplaceMedium{2, T1},
    x::SVector{2, T2}, 
    outward_normal::SVector
) where {T1, T2}
    
    # Calculate the purely real scalar value
    val = -1 / (2π) * log(norm(x))
    
    # Return it wrapped in a 1x1 Static Matrix of type T2
    return SMatrix{1, 1, T2}(val)
end

function greens_gradient(
    field::DirichletType, 
    medium::LaplaceMedium{2, T1}, 
    x::SVector{2, T2}, 
    outward_normal::SVector
) where {T1, T2}
    
    r2 = x[1]^2 + x[2]^2 + 1e-12 
    
    grad_x = -T2(1) / T2(2π) * (x[1] / r2)
    grad_y = -T2(1) / T2(2π) * (x[2] / r2)
    
   
    return SVector{2, T2}(grad_x, grad_y)
end


