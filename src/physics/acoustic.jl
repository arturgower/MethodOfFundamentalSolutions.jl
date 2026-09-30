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

# The gradients of the two acoustic kernels with respect to the SENSOR position, i.e. with
# respect to `x = sensor - source` at a fixed outward normal — what the boundary-uncertainty
# machinery of the VariationalBayesianSolver needs (see `system_matrix_gradient`). Both are
# returned in the 1 × 1 × Dim layout of the elastic kernel, so that the single field
# component of the acoustic medium is indexed exactly like a multi-component one.

function greens_gradient(field::TractionType, medium::Acoustic{T,2}, x::AbstractVector, outward_normal::AbstractVector = zeros(T,2)) where T

    k = medium.ω / medium.c
    r = norm(x)

    # ∇φ = φ'(r) x/r with φ'(r) = (i k/4) H₁¹(k r)
    dφ = (im * k / 4) * hankelh1(1, k * r) / r

    G = zeros(Complex{T}, 1, 1, 2)
    G[1,1,1] = dφ * x[1]
    G[1,1,2] = dφ * x[2]
    return G
end

function greens_gradient(field::DisplacementType, medium::Acoustic{T,2}, x::AbstractVector, outward_normal::AbstractVector = zeros(T,2)) where T

    k = medium.ω / medium.c
    r = norm(x)
    xn = dot(x, outward_normal)

    # the kernel is A(r) (x⋅n) with A(r) = (i k)/(4ρ) H₁¹(k r)/r, so
    # ∂/∂x_d = A'(r) (x_d/r)(x⋅n) + A(r) n_d, and with H₁¹'(z) = H₀¹(z) - H₁¹(z)/z
    #     A'(r) = (i k)/(4ρ) [k H₀¹(k r)/r - 2 H₁¹(k r)/r²].
    A  = (im * k / (4 * medium.ρ)) * hankelh1(1, k * r) / r
    dA = (im * k / (4 * medium.ρ)) * (k * hankelh1(0, k * r) / r - 2 * hankelh1(1, k * r) / r^2)

    G = zeros(Complex{T}, 1, 1, 2)
    G[1,1,1] = dA * (x[1] / r) * xn + A * outward_normal[1]
    G[1,1,2] = dA * (x[2] / r) * xn + A * outward_normal[2]
    return G
end
