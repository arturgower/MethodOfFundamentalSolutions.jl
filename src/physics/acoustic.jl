"""
    Acoustic{T<:AbstractFloat,Dim}(ρ::T, c::Complex{T})
    Acoustic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic acoustic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (1D) fields in Dim dimensions.
"""
struct Acoustic{T,Dim} <: PhysicalMedium{Dim,1}
    ω::T # Angular frequency
    ρ::T # Densityω::T
    c::Complex{T} # Phase velocity
end

# Constructor which supplies the dimension without explicitly mentioning type
Acoustic(ω::T, ρ::T, c::Union{T,Complex{T}}, Dim::Integer) where {T<:Number} =  Acoustic{T,Dim}(ω, ρ, Complex{T}(c))
Acoustic(Dim::Integer;  ω::T = 0.0, ρ::T = 0.0, c::Union{T,Complex{T}} = 0.0) where {T<:Number} =  Acoustic{T,Dim}(ω, ρ, Complex{T}(c))

name(a::Acoustic{T,Dim}) where {Dim,T} = "$(Dim)D Acoustic"

# The two acoustic boundary fields of the fundamental solution φ(x) = -(i/4) H₀¹(k |x|),
# k = ω/c. Across a transmission (penetrable) interface both are continuous:
# - TractionType: the field φ itself (the acoustic pressure);
# - DisplacementType: the normal displacement (1/ρ) ∂φ/∂n.

function greens(field::TractionType, medium::Acoustic{T,2}, x::SVector{2,T}, outward_normal::AbstractVector{T} = zeros(T,2)) where T

    G = zeros(Complex{T},1,1)
    G[1,1] = -(im/4) * hankelh1(0, (medium.ω / medium.c) * norm(x))
    return G
end

function greens(field::DisplacementType, medium::Acoustic{T,2}, x::SVector{2,T}, outward_normal::AbstractVector{T} = zeros(T,2)) where T

    k = medium.ω / medium.c
    r = norm(x)

    # ∂φ/∂n = φ'(r) ∂r/∂n with φ'(r) = (i k/4) H₁¹(k r) (since H₀¹' = -H₁¹) and
    # ∂r/∂n = x⋅n / r; the density weight 1/ρ makes it continuous across an interface
    G = zeros(Complex{T},1,1)
    G[1,1] = (im * k / 4) * hankelh1(1, k * r) * dot(x, outward_normal) / (r * medium.ρ)
    return G
end
