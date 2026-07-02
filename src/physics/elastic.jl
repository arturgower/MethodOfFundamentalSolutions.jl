"""
    Elastostatic{Dim,T<:AbstractFloat}(ρ::T, c::T)
    Elastostatic(ρ::T, c::T, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium for static simulations.

Simulations in this medium produce scalar (Dim) fields in Dim dimensions.
"""
struct Elastostatic{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::T # Phase velocity of pressure wave
    cs::T # Phase velocity of shear wave
    
    # Constructor which supplies the dimension without explicitly mentioning type
    function Elastostatic(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
        # check Lame parameters are positive
        μ = ρ * cs^2  
        λ = ρ * cp^2 - 2μ

        if μ < 0 || λ < 0 
            @error "The Lame parameters need to be positive."
        end

        new{Dim,T}(ρ,T(cp),T(cs))
    end
end

struct DisplacementType <: FieldType end
# DisplacementType(field::AbstractVector) = DisplacementType{length(field)}(SVector(field...))

struct TractionType <: FieldType end
    # field::SVector{Dim,Float64}
    # outward_normal::SVector{Dim,Float64}

# function TractionType(field::AbstractVector, outward_normal::AbstractVector)
#     Dim = length(field)
#     return TractionType{Dim}(SVector(field...), SVector(outward_normal...) ./ norm(outward_normal))
# end
struct StrainType <: FieldType end
    # field::SVector{Dim,Float64}
    # strain_direction::SVector{Dim,Float64}

struct ParticularGravity{T} <: ParticularSolution 
    g::T
    height::T # the height times the average width should equal the volume.
    function ParticularGravity(; g::T = 9.81, height::T = zero(typeof(g))) where T
        new{T}(g, height)
    end
end

function greens(traction::TractionType, medium::Elastostatic{2,T}, x::SVector{2,T}, outward_normal::SVector) where T
    n = outward_normal
    μ = medium.ρ * medium.cs^2  
    λ = medium.ρ * medium.cp^2 - 2μ

    r2 = dot(x,x)
    xn = dot(x,n)

    Σn = [
        μ * (xn * (l == i) + x[i]*n[l] - x[l]*n[i]) + 2(λ+μ) * (x[l]*x[i] * xn / r2)
    for i = 1:2, l = 1:2]
    Σn = Σn ./ (2pi*(λ+2μ)*r2)

    return Σn 
end

function greens(displace::DisplacementType, medium::Elastostatic{2,T}, x::SVector{2,T}, outward_normal) where T
    μ = medium.ρ * medium.cs^2  
    λ = medium.ρ * medium.cp^2 - 2μ

    ν = λ / (2λ + 2μ) 
    r2 = dot(x,x)
    r = sqrt(r2)

    U = [
        (3-4ν) * log(r) * (i==j) - x[i] * x[j] / r2
    for i = 1:2, j = 1:2]

    U = U ./ (8pi*μ*(1-ν))
    return U
end

function greens(strain::StrainType, medium::Elastostatic{2,T}, x::SVector{2,T}, strain_direction::SVector) where T
    s = strain_direction
    μ = medium.ρ * medium.cs^2  
    λ = medium.ρ * medium.cp^2 - 2μ

    ν = λ / (2λ + 2μ)
    r2 = dot(x,x)
    xs = dot(x,s)

    Es = [
         (3-4ν)*( xs*(l==i) + x[i]*s[l] )- 2*x[l]*s[i]  - s[l]*x[i] - xs*(l==i) + 4*xs*x[i]*x[l]/r2
    for i = 1:2, l = 1:2]
    
    Es = Es ./ (16pi*μ*(1-ν)*r2)
    
    return Es

end

function field(traction::TractionType, medium::Elastostatic{Dim,T}, psol::ParticularGravity, x::AbstractVector, outward_normal::AbstractVector) where {T,Dim}
    σzz = medium.ρ*psol.g*(x[Dim] - psol.height)
    return SVector(zero(T), σzz * outward_normal[Dim])
end


function greens_gradient(
    traction::TractionType, 
    medium::Elastostatic{2,T}, 
    x::SVector{2,T}, 
    outward_normal::SVector
) where T

    n = outward_normal
    μ = medium.ρ * medium.cs^2  
    λ = medium.ρ * medium.cp^2 - 2μ

    r2 = dot(x,x)
    xn = dot(x,n)

    # Pre-calculate inverse powers of r to avoid repetitive divisions
    inv_r2 = one(T) / r2
    inv_r4 = inv_r2 * inv_r2
    inv_r6 = inv_r4 * inv_r2

    # Output tensor shape: (i, l, k) matching ∂(Σn_il)/∂x_k
    grad_Σn = [
        begin
            # Inline Kronecker Deltas converted to type T
            δ_li = T(l == i)
            δ_ik = T(i == k)
            δ_lk = T(l == k)

            # Term 1: Associated with 1/r^2 scaling
            term1 = μ * (n[k] * δ_li + δ_ik * n[l] - δ_lk * n[i]) * inv_r2

            # Term 2: Associated with 1/r^4 scaling
            term2 = 2 * (
                -μ * x[k] * (xn * δ_li + x[i] * n[l] - x[l] * n[i]) +
                (λ + μ) * (δ_lk * x[i] * xn + x[l] * δ_ik * xn + x[l] * x[i] * n[k])
            ) * inv_r4

            # Term 3: Associated with 1/r^6 scaling
            term3 = -8 * (λ + μ) * (x[l] * x[i] * xn * x[k]) * inv_r6

            # Apply global scalar normalization factor
            (term1 + term2 + term3) / (2 * pi * (λ + 2μ))
        end
    for i = 1:2, l = 1:2, k = 1:2]

    return grad_Σn 
end