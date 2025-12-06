
struct DirichletType <: FieldType end
#Dirchlet type of boundary conditions
struct NeumannType <: FieldType end
#Neumann type of boundary conditions


function greens(field::DirichletType, medium::Acoustic{T,2}, x::SVector{2,T}, outward_normal; ω::Float64= 2pi * 1.0) where T
   
    G = zeros(Complex{T},1,1)
    G[1,1] = -(im/4) * hankelh1(0, (ω / medium.c) * norm(x))
    return G
end
