# Variational evidence maximization for the Method of Fundamental Solutions.
#
# Implements Algorithm 1 ("Variational evidence maximization: learning the basis χ and the
# prior {αᵢ}, with known measurement noise σ², and predicting the posterior of the
# coefficients a") of docs/theory/main-variational-evidence.tex.
#
# The priors of the coefficients a and of the boundary points x are diagonal:
#   p(a)  = ∏ᵢ N(aᵢ | 0, 1/αᵢ)          (ARD prior, precisions αᵢ learned by EM)
#   p(δx) = N(0, Σ_x),  Σ_x diagonal     (taken from the covariance of the boundary points)
# The measurement noise Σ = diag(σ²) is known and fixed throughout.
#
# Complex-valued problems (e.g. acoustics) are handled by stacking real and imaginary
# parts: g̃ = [Re g; Im g], M̃ = [Re M  -Im M; Im M  Re M], ã = [Re a; Im a], so that the
# whole algorithm runs on a real linear-Gaussian model.

"""
    VariationalBayesianSolver <: AbstractSolver

Solver implementing the variational evidence maximization of Algorithm 1 in
docs/theory/main-variational-evidence.tex: it infers the posterior of the source
coefficients `a`, learns the ARD prior precisions `αᵢ` (automatic relevance
determination, pruning superfluous sources), optionally learns the source positions `χ`,
and optionally infers the boundary perturbation `δx`. The measurement noise is known.

The measurement noise is specified through the boundary data, not through the solver: give
the boundary `fields` as a vector of `MvNormal` (one per boundary point) or as a single
`MvNormal` over the flattened fields, whose covariance is the sensor noise. For
complex-valued problems each point must be a `2FD`-dimensional `MvNormal` over the stacked
real and imaginary parts `[Re; Im]` of its `FD` field components.

The options shared with [`BayesianSolver`](@ref) live in the field `options::SolverOptions`
(see [`SolverOptions`](@ref)); the keyword constructor accepts them as keywords directly.

# Keyword arguments
- `priors`: optional vector of `MvNormal`, one per source; the (diagonal) covariance of
  `priors[j]` sets the initial prior variances `1/αᵢ` of the coefficients of source `j`
  (the means are ignored: the ARD prior is zero-mean).
- `prior`: optional `MvNormal`; its (diagonal) covariance sets the initial prior variances `1/αᵢ`.
- `prior_variance = 1.0`: initial prior variance of the coefficients. A scalar is shared by
  all coefficients; a vector must have one entry per coefficient (for complex problems either
  one entry per complex coefficient, or per real degree of freedom `[Re; Im]`).
- `optimise_source_positions_flag = false`: run the M-step over the source positions `χ`
  (a few L-BFGS steps per iteration, accepted only if the bound increases).
- `update_geometry_flag = false`: update the boundary factor `q(δx)` and re-center the
  boundary (E-step II). Requires the `boundary_points` to be an `MvNormal`; its diagonal
  covariance is the prior `Σ_x`. Not implemented for complex-valued problems.
- `learn_prior_flag = true`: run the EM update `αᵢ = 1/(μᵢ² + Σᵢᵢ)` of the prior precisions.
- `ard_prune_flag = true`: remove sources all of whose precisions exceed `ard_threshold`.
- `ard_threshold = 1e8`: precision above which a coefficient counts as switched off.
- `mackay_acceleration_flag = true`: use the MacKay/Tipping fixed-point update
  `αᵢ = γᵢ/μᵢ²` as an acceleration, falling back to the EM update whenever the bound decreases.
- `use_greens_gradient_analytical_flag = true`: use the analytic gradient of the expected
  misfit for the source-position M-step when `greens_gradient` is available (real problems,
  fixed geometry); otherwise finite differences are used.
- `elbo_tol = 1e-8`: relative tolerance on the increase of the evidence lower bound.
- `max_iters = 200`: maximum number of variational EM iterations.
- `source_position_iters = 5`: L-BFGS iterations per source-position M-step.
"""
struct VariationalBayesianSolver <: AbstractSolver
    options::SolverOptions
    prior_variance::Vector{Float64}
    ard_prune_flag::Bool
    ard_threshold::Float64
    mackay_acceleration_flag::Bool
    elbo_tol::Float64
end

function VariationalBayesianSolver(;
        priors::AbstractVector{<:AbstractMvNormal} = MvNormal[],
        prior::Union{Nothing, ContinuousMultivariateDistribution} = nothing,
        prior_variance::Union{Real, AbstractVector{<:Real}} = 1.0,
        optimise_source_positions_flag::Bool = false,
        update_geometry_flag::Bool = false,
        learn_prior_flag::Bool = true,
        ard_prune_flag::Bool = true,
        ard_threshold::Real = 1e8,
        mackay_acceleration_flag::Bool = true,
        use_greens_gradient_analytical_flag::Bool = true,
        elbo_tol::Real = 1e-8,
        max_iters::Int = 200,
        source_position_iters::Int = 5
    )

    # The diagonal ARD prior is fully characterized by its variances, so the `priors`,
    # `prior` and `prior_variance` keyword forms all reduce to the same variance vector
    # (the means of the MvNormals are ignored: the prior is zero-mean).
    pv = if !isempty(priors)
        vcat([diag(cov(p)) for p in priors]...)
    elseif prior !== nothing
        diag(cov(prior))
    else
        prior_variance
    end
    pvs = pv isa Real ? [Float64(pv)] : Vector{Float64}(pv)
    all(>(0), pvs) || throw(ArgumentError("all prior variances must be positive"))

    options = SolverOptions(;
        optimise_source_positions_flag = optimise_source_positions_flag,
        use_greens_gradient_analytical_flag = use_greens_gradient_analytical_flag,
        update_geometry_flag = update_geometry_flag,
        learn_prior_flag = learn_prior_flag,
        max_iters = max_iters,
        source_position_iters = source_position_iters
    )

    return VariationalBayesianSolver(
        options, pvs,
        ard_prune_flag, Float64(ard_threshold),
        mackay_acceleration_flag, Float64(elbo_tol)
    )
end

VariationalBayesianSolver(prior::ContinuousMultivariateDistribution; kws...) =
    VariationalBayesianSolver(; prior = prior, kws...)

"""
    VariationalSolution

Result of `solve` with a [`VariationalBayesianSolver`](@ref).

All the information about the coefficients lives in `fsol`, and all the information about
the boundary lives in `boundary_shape`.

# Fields
- `fsol::FundamentalSolution`: the coefficient posterior q(a) = N(μ_post, Σ_post) at the
  learned source positions: `fsol.coefficients` is the posterior mean and
  `fsol.coefficients_covariance` the posterior covariance over the retained degrees of
  freedom. Use it with `field`, `field_covariance` and `field_std` (which also accept the
  `VariationalSolution` directly). For complex problems the coefficients are complex and
  their posterior covariance is stored over the stacked real degrees of freedom
  `[Re a; Im a]`; `field_covariance` then returns the covariance of the stacked field
  `[Re f; Im f]`.
- `boundary_shape::BoundaryShape`: the posterior of the boundary. When the geometry was
  updated (`update_geometry_flag`) its `boundary_points` is an `MvNormal` whose mean is the
  re-centered boundary and whose covariance is the posterior covariance Σ_δx; otherwise it
  is the boundary shape of the input `BoundaryData`, unchanged.
- `prior_precisions`: the learned ARD precisions αᵢ of the retained coefficients.
- `elbo_history`: the evidence lower bound F after each iteration.
- `baseline_resets`: iterations at which the model changed significantly (source pruning or
  a large boundary re-centering), across which values of F are not comparable.
- `recenter_iterations`: every iteration at which the boundary was re-centered. F is only
  guaranteed to be monotone between, not across, re-centerings; for well-linearized
  problems the change of F across a small re-centering is negligible.
- `misfit_ratio`: the noise-weighted expected misfit R/N_g of eq. (R); ≈ 1 at convergence
  when the model is consistent with the known noise level.
"""
struct VariationalSolution{FS <: FundamentalSolution, BS <: BoundaryShape, T <: Real}
    fsol::FS
    boundary_shape::BS
    prior_precisions::Vector{T}
    elbo_history::Vector{T}
    baseline_resets::Vector{Int}
    recenter_iterations::Vector{Int}
    misfit_ratio::T
end

field(ft::FieldType, vsol::VariationalSolution, x::AbstractVector, outward_normal::AbstractVector = ones(x |> length)) =
    field(ft, vsol.fsol, x, outward_normal)
field_covariance(ft::FieldType, vsol::VariationalSolution, x::AbstractVector, outward_normal::AbstractVector = ones(x |> length)) =
    field_covariance(ft, vsol.fsol, x, outward_normal)
field_std(ft::FieldType, vsol::VariationalSolution, x::AbstractVector, outward_normal::AbstractVector = ones(x |> length)) =
    field_std(ft, vsol.fsol, x, outward_normal)

"""
    grid_source_positions(bd::BoundaryData; n = 15, scale = 2.0, clearance = 1.0)

Candidate MFS source positions "everywhere": a regular `n × n` grid covering the bounding
box of the boundary enlarged by `scale`, keeping only points outside the domain and further
than `clearance` times the average boundary spacing from the boundary. Intended as an
overcomplete initialization for a [`VariationalBayesianSolver`](@ref), whose automatic
relevance determination then switches off the unnecessary sources.
"""
function grid_source_positions(bd::BoundaryData{F, 2}; n::Int = 15, scale::Real = 2.0, clearance::Real = 1.0) where F
    pts = mean_points(bd)
    len = length(pts)

    xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
    centre = SVector((minimum(xs) + maximum(xs)) / 2, (minimum(ys) + maximum(ys)) / 2)
    halfwidth = SVector(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) ./ 2

    # average spacing between neighbouring boundary points, as in source_positions
    sampled_rng = LinRange(1, len, min(6, len)) .|> round .|> Int
    spacing = mean(map(pts[sampled_rng]) do p
        dists = [norm(p - q) for q in pts]
        idx = sortperm(dists)[2:min(3, len)]
        mean(dists[idx])
    end)

    grid = [
        centre + SVector(2u - 1, 2v - 1) .* (scale .* halfwidth)
    for u in LinRange(0, 1, n), v in LinRange(0, 1, n)]

    return filter(vec(grid)) do p
        p ∉ bd && minimum(norm(p - q) for q in pts) > clearance * spacing
    end
end

# ------------------------------------------------------------------------------------
# Internal helpers. All of them operate on the real (possibly [Re; Im]-stacked) model.
# ------------------------------------------------------------------------------------

_structured_positions(chiflat::AbstractVector, Dim::Int) =
    [SVector{Dim, Float64}(ntuple(k -> chiflat[i + k - 1], Dim)) for i in 1:Dim:length(chiflat)]

# The mean data and the known measurement-noise covariance, from the boundary data alone,
# in the ordering of the real working model: per-point stacking for real problems, and
# complex means with `[all Re; all Im]`-stacked variances for complex ones. The noise is
# returned as a vector of per-component variances when it is diagonal (independent sensor
# noise), or as a full covariance matrix when the boundary `fields` carry a correlated
# covariance. For complex problems each field point must be an `MvNormal` over the stacked
# real and imaginary parts `[Re; Im]` of its FD field components (diagonal noise only).
function _data_and_noise(bd::BoundaryData, FD::Int, iscomplex::Bool)
    fields = bd.fields
    noise_error = ArgumentError(
        "the measurement noise must be known: give the boundary `fields` as a vector of " *
        "`MvNormal` (one per boundary point), or a single `MvNormal` over the flattened " *
        "fields, whose covariance is the sensor noise"
    )

    if iscomplex
        fields isa AbstractVector{<:AbstractMvNormal} || throw(noise_error)
        all(length(d) == 2FD for d in fields) || throw(ArgumentError(
            "for complex-valued problems each field point must be an `MvNormal` of dimension " *
            "2 × $FD: the stacked real and imaginary parts `[Re; Im]` of its field"
        ))
        means = mean.(fields)
        vars = [diag(cov(d)) for d in fields]
        g0 = vcat([m[1:FD] .+ im .* m[(FD + 1):2FD] for m in means]...)
        s2 = vcat([v[1:FD] for v in vars]..., [v[(FD + 1):2FD] for v in vars]...)
        return g0, Vector{Float64}(s2)
    end

    Σf = cov(fields)
    Σf isa UniformScaling && throw(noise_error)
    g0 = Vector{Float64}(flat_fields(fields))
    Σf = Matrix{Float64}(Σf)
    # keep the fast diagonal path for independent sensor noise; otherwise use the full
    # (correlated) covariance
    noise = isdiag(Σf) ? diag(Σf) : Symmetric(Σf)
    return g0, noise
end

# The model matrix of the real working model: for complex physics the real/imaginary parts
# are stacked so that g̃ = M̃ ã with ã = [Re a; Im a].
function _stacked_system_matrix(source_positions, medium, bd, iscomplex::Bool)
    M0 = system_matrix(source_positions, medium, bd)
    if iscomplex
        Mr = real.(M0); Mi = imag.(M0)
        return [Mr -Mi; Mi Mr]
    else
        return Matrix{Float64}(M0)
    end
end

# Column indices, in the stacked coefficient vector, of the sources `keep`.
function _kept_columns(keep::AbstractVector{Int}, n_src::Int, FD::Int, iscomplex::Bool)
    re = reduce(vcat, [collect(((j - 1) * FD + 1):(j * FD)) for j in keep])
    return iscomplex ? vcat(re, re .+ n_src * FD) : re
end

_source_columns(j::Int, n_src::Int, FD::Int, iscomplex::Bool) = _kept_columns([j], n_src, FD, iscomplex)

function _initial_precisions(pv::Vector{Float64}, K::Int, Kh::Int, iscomplex::Bool)
    if length(pv) == 1
        return fill(1 / pv[1], K)
    elseif length(pv) == K
        return 1 ./ pv
    elseif iscomplex && length(pv) == Kh
        return vcat(1 ./ pv, 1 ./ pv)
    else
        throw(ArgumentError("prior_variance must be a scalar, or have one entry per coefficient ($Kh) or per real degree of freedom ($K)"))
    end
end

# The measurement-noise precision Σ⁻¹ applied to a vector or matrix `X`. The noise is stored
# either as a vector of per-component variances' reciprocals `w = 1/σ²` (independent sensor
# noise, the fast path) or as a full precision matrix `w = Σ⁻¹` (correlated noise).
_apply_precision(w::AbstractVector, X) = w .* X
_apply_precision(w::AbstractMatrix, X) = w * X

const NoisePrecision = Union{AbstractVector, AbstractMatrix}

# E-step I, eq. (vb_coefficient_update): q(a) = N(μ_post, Σ_post) with
# Σ_post = (M̄ᵀ Σ⁻¹ M̄ + Γ + diag α)⁻¹ and μ_post = Σ_post M̄ᵀ Σ⁻¹ g.
function _coefficient_posterior(M::AbstractMatrix, Γ, α::AbstractVector, w::NoisePrecision, g::AbstractVector)
    H = Matrix(M' * _apply_precision(w, M))
    Γ === nothing || (H .+= Γ)
    @inbounds for i in eachindex(α)
        H[i, i] += α[i]
    end

    C = cholesky(Symmetric(H); check = false)
    if !issuccess(C)
        jitter = 1e-12 * max(1.0, maximum(abs, diag(H)))
        C = cholesky(Symmetric(H + jitter * I))
    end
    logdetΣ = -logdet(C)
    Σpost = Matrix(Symmetric(inv(C)))
    μ = Σpost * (M' * _apply_precision(w, g))
    return μ, Σpost, logdetΣ
end

# The noise-weighted expected misfit R of eq. (R):
# R = E_q ‖g - M a - D(a) δx‖²_Σ⁻¹, evaluated after re-centering (μ_δx = 0), so that
# R = (g - Mμ)ᵀΣ⁻¹(g - Mμ) + tr(Σ⁻¹ M Σ_post Mᵀ) + μᵀΓμ + tr(Γ Σ_post).
function _expected_misfit(M::AbstractMatrix, Γ, μ::AbstractVector, Σpost::AbstractMatrix, w::NoisePrecision, g::AbstractVector)
    r = g - M * μ
    S = M * Σpost
    R = dot(r, _apply_precision(w, r)) + sum(_apply_precision(w, M) .* S)
    if Γ !== nothing
        R += dot(μ, Γ * μ) + dot(Γ, Σpost)
    end
    return R
end

# The evidence lower bound F of eq. (elbo_model), with the Gaussian KL of eq. (gaussian_kl).
# The prior of the boundary perturbation stays anchored at the original nominal boundary:
# relative to the current linearization point it is N(m0, Σx) with m0 = μ_x0 - x_lin.
function _elbo(R::Real, Nr::Int, noise_logdet::Real, α::AbstractVector, μ::AbstractVector,
        Σpost::AbstractMatrix, logdetΣ::Real;
        μδx = nothing, Σδx = nothing, logdetΣδx::Real = 0.0, sx2 = nothing, m0 = nothing
    )
    F = -0.5 * (Nr * log(2π) + noise_logdet) - 0.5 * R
    K = length(α)
    F -= 0.5 * (sum(α .* (abs2.(μ) .+ diag(Σpost))) - K - sum(log.(α)) - logdetΣ)
    if μδx !== nothing
        nx = length(μδx)
        F -= 0.5 * (sum(diag(Σδx) ./ sx2) + sum(abs2.(μδx .- m0) ./ sx2) - nx + sum(log.(sx2)) - logdetΣδx)
    end
    return F
end

# The extra coefficient precision Γ = Σ_kl Σδx^kl M_kᵀ Σ⁻¹ M_l of eq. (Gamma). Row r of M
# depends only on the boundary point of its own sensor, so only the Dim × Dim diagonal
# blocks of Σδx contribute.
function _gamma(gradM::AbstractArray, Σδx::AbstractMatrix, w::AbstractVector, d_m::Int, K::Int)
    N = size(gradM, 1)
    Dim = size(gradM, 3)
    Γ = zeros(K, K)
    for r in 1:N
        i = (r - 1) ÷ d_m + 1
        for d in 1:Dim, d2 in 1:Dim
            c = w[r] * Σδx[(i - 1) * Dim + d, (i - 1) * Dim + d2]
            c == 0 && continue
            vd = Vector(gradM[r, :, d])
            vd2 = Vector(gradM[r, :, d2])
            mul!(Γ, vd, transpose(vd2), c, 1.0)
        end
    end
    return Matrix(Symmetric(Γ))
end

# E-step II, eqs. (vb_boundary_precision)-(vb_boundary_update): q(δx) = N(μ_δx, Σ_δx).
# With a diagonal prior Σ_x and independent sensors both the precision E[DᵀΣ⁻¹D] and Σ_δx
# are block diagonal, one Dim × Dim block per sensor. The prior of δx relative to the
# current linearization point is N(m0, Σx) with m0 = μ_x0 - x_lin, so that re-centering
# does not move the prior: directions the data cannot inform are pulled back to the
# nominal boundary instead of drifting with each re-linearization.
function _boundary_posterior(M::AbstractMatrix, gradM::AbstractArray, μ::AbstractVector,
        Σpost::AbstractMatrix, w::AbstractVector, g::AbstractVector,
        sx2::AbstractVector, m0::AbstractVector, d_m::Int, Dim::Int
    )
    N = size(M, 1)
    n_sensors = N ÷ d_m
    nx = n_sensors * Dim
    μδx = zeros(nx)
    Σδx = zeros(nx, nx)
    logdetΣδx = 0.0
    Mμ = M * μ

    for i in 1:n_sensors
        A = zeros(Dim, Dim)
        b = zeros(Dim)
        for r in (((i - 1) * d_m + 1):(i * d_m))
            Mr = M[r, :]
            ΣMr = Σpost * Mr
            vs = [Vector(gradM[r, :, d]) for d in 1:Dim]
            Σvs = [Σpost * vs[d] for d in 1:Dim]
            for d in 1:Dim
                vμd = dot(vs[d], μ)
                b[d] += w[r] * (vμd * (g[r] - Mμ[r]) - dot(vs[d], ΣMr))
                for d2 in d:Dim
                    val = w[r] * (vμd * dot(vs[d2], μ) + dot(vs[d], Σvs[d2]))
                    A[d, d2] += val
                    d2 > d && (A[d2, d] += val)
                end
            end
        end
        ks = ((i - 1) * Dim + 1):(i * Dim)
        Sblock = inv(Symmetric(A + Diagonal(1 ./ sx2[ks])))
        Σδx[ks, ks] = Sblock
        μδx[ks] = Sblock * (b + m0[ks] ./ sx2[ks])
        logdetΣδx += logdet(Symmetric(Matrix(Sblock)))
    end
    return μδx, Σδx, logdetΣδx
end

# Re-center the boundary at the current estimate: μ_x ← μ_x + μ_δx (Algorithm 1, line 6).
# The outward normals are kept fixed, consistent with small perturbations.
function _recenter_boundary(bd::BoundaryData, μδx::AbstractVector, Dim::Int)
    pts = mean_points(bd)
    newpts = [
        pts[i] + SVector{Dim, Float64}(ntuple(d -> μδx[(i - 1) * Dim + d], Dim))
    for i in eachindex(pts)]

    shape = bd.boundary_shape
    return BoundaryData(bd.fieldtype,
        BoundaryShape(newpts, shape.normals, shape.interior_points),
        bd.fields
    )
end

# M-step objective over the source positions χ: F depends on χ only through the expected
# misfit R(χ) (eq. (R)), so minimizing R maximizes F at fixed q.
function _chi_misfit(chiflat::AbstractVector, medium, bd, iscomplex::Bool,
        μ::AbstractVector, Σpost::AbstractMatrix, w::NoisePrecision, g::AbstractVector, Dim::Int;
        Γ = nothing, Σδx = nothing, d_m::Int = 1
    )
    pos = _structured_positions(chiflat, Dim)
    M = _stacked_system_matrix(pos, medium, bd, iscomplex)
    Γx = if Σδx === nothing
        Γ
    else
        gradM = system_matrix_gradient(pos, medium, bd)
        _gamma(gradM, Σδx, w, d_m, size(M, 2))
    end
    return _expected_misfit(M, Γx, μ, Σpost, w, g)
end

# Analytic gradient of _chi_misfit for real problems with fixed geometry, eq. (chi_gradient).
# ∂M/∂(source j, coord d) is minus the derivative with respect to the boundary coordinate,
# which is what system_matrix_gradient returns, and only the columns of source j are nonzero.
function _chi_misfit_gradient!(G::AbstractVector, chiflat::AbstractVector, medium, bd,
        μ::AbstractVector, Σpost::AbstractMatrix, w::AbstractVector, g::AbstractVector,
        Dim::Int, FD::Int
    )
    pos = _structured_positions(chiflat, Dim)
    M = Matrix{Float64}(system_matrix(pos, medium, bd))
    gradM = system_matrix_gradient(pos, medium, bd)

    wr = w .* (g - M * μ)
    SM = M * Σpost

    for j in eachindex(pos), d in 1:Dim
        cols = ((j - 1) * FD + 1):(j * FD)
        P = @view gradM[:, cols, d]   # ∂M/∂χ_{j,d} = -P on these columns
        t1 = 2 * dot(P * view(μ, cols), wr)
        t2 = -2 * sum(w .* vec(sum(P .* view(SM, :, cols), dims = 2)))
        G[(j - 1) * Dim + d] = t1 + t2
    end
    return G
end

_has_greens_gradient(ft, medium, Dim) =
    hasmethod(greens_gradient, Tuple{typeof(ft), typeof(medium), SVector{Dim, Float64}, SVector{Dim, Float64}})

# MacKay/Tipping fixed-point update α = γ/μ², eq. (mackay_updates), with the EM update as a
# safe fallback for undetermined components.
function _mackay_precisions(α::AbstractVector, μ::AbstractVector, Σpost::AbstractMatrix, α_em::AbstractVector)
    d = diag(Σpost)
    return map(eachindex(α)) do i
        γ = 1 - α[i] * d[i]
        μ2 = abs2(μ[i])
        if γ <= 0
            α_em[i]
        elseif μ2 < 1e-300
            1e12
        else
            γ / μ2
        end
    end
end

# ------------------------------------------------------------------------------------
# Algorithm 1
# ------------------------------------------------------------------------------------

function solve(sim::Simulation{VariationalBayesianSolver, Dim}) where Dim

    solver = sim.solver
    medium = sim.medium
    FD = field_dimension(medium)
    bd = sim.boundary_data

    src_pos = Vector{SVector{Dim, Float64}}(sim.source_positions)
    M0 = system_matrix(src_pos, medium, bd)
    iscomplex = eltype(M0) <: Complex

    # --- data g (subtracting any particular solution) and its known noise variance,
    #     both taken from the boundary data ---
    g0, noise = _data_and_noise(bd, FD, iscomplex)
    g_particular = field(medium, bd, sim.particular_solution)
    g0 = g0 - vcat(g_particular...)

    g = iscomplex ? Vector{Float64}(vcat(real.(g0), imag.(g0))) : Vector{Float64}(g0)
    Nr = length(g)

    # measurement-noise precision w (= Σ⁻¹) and log|Σ|, from either a diagonal (vector of
    # variances) or a full (matrix) noise covariance
    w, noise_logdet = if noise isa AbstractVector
        all(>(0), noise) || throw(ArgumentError("all measurement noise variances must be positive"))
        (1 ./ noise, sum(log, noise))
    else
        Cf = cholesky(Symmetric(Matrix(noise)))
        (Matrix(Symmetric(inv(Cf))), logdet(Cf))
    end

    # --- geometry prior Σ_x (diagonal) ---
    Σx_raw = cov(bd.boundary_shape.boundary_points)
    do_geometry = solver.options.update_geometry_flag && !(Σx_raw isa UniformScaling)
    if do_geometry && iscomplex
        throw(ArgumentError("geometry updates are not implemented for complex-valued problems"))
    end
    if do_geometry && !(w isa AbstractVector)
        throw(ArgumentError("geometry updates require independent (diagonal) sensor noise; the boundary `fields` covariance must be diagonal"))
    end
    sx2 = do_geometry ? Vector{Float64}(diag(Σx_raw)) : nothing

    bd_current = bd
    n_sensors = length(mean_points(bd))
    d_m = size(M0, 1) ÷ n_sensors

    # --- initial prior precisions and model matrix ---
    n_src = length(src_pos)
    Kh = n_src * FD
    K = iscomplex ? 2Kh : Kh
    α = _initial_precisions(solver.prior_variance, K, Kh, iscomplex)
    α_cap = 1e12

    M = iscomplex ? _stacked_system_matrix(src_pos, medium, bd_current, true) : Matrix{Float64}(M0)

    gradM = do_geometry ? system_matrix_gradient(src_pos, medium, bd_current) : nothing
    Σδx = do_geometry ? Matrix(Diagonal(sx2)) : nothing            # initialize q(δx) at the prior
    logdetΣδx = do_geometry ? sum(log.(sx2)) : 0.0
    μδx = do_geometry ? zeros(length(sx2)) : nothing
    # prior mean of δx relative to the current linearization point: the prior of the
    # boundary stays anchored at the original nominal boundary μ_x0 across re-centerings
    x_nominal = do_geometry ? Vector{Float64}(vcat(mean_points(bd)...)) : nothing
    m0 = do_geometry ? zeros(length(sx2)) : nothing
    Γ = do_geometry ? _gamma(gradM, Σδx, w, d_m, K) : nothing

    elbo_history = Float64[]
    baseline_resets = Int[]
    recenter_iterations = Int[]
    F_prev = -Inf
    μ = zeros(K)
    Σpost = Matrix{Float64}(I, K, K)
    logdetΣ = 0.0

    for it in 1:solver.options.max_iters
        reset_baseline = false

        # --- E-step I: update q(a), eq. (vb_coefficient_update) ---
        μ, Σpost, logdetΣ = _coefficient_posterior(M, Γ, α, w, g)

        # --- E-step II: update q(δx) and re-center, eq. (vb_boundary_update) ---
        if do_geometry
            μδx, Σδx, logdetΣδx = _boundary_posterior(M, gradM, μ, Σpost, w, g, sx2, m0, d_m, Dim)
            step = norm(μδx, Inf)
            prior_scale = sqrt(maximum(sx2))
            if step > 1e-6 * prior_scale
                bd_current = _recenter_boundary(bd_current, μδx, Dim)
                M = _stacked_system_matrix(src_pos, medium, bd_current, iscomplex)
                gradM = system_matrix_gradient(src_pos, medium, bd_current)
                μδx = zero(μδx)
                m0 = x_nominal - Vector{Float64}(vcat(mean_points(bd_current)...))
                push!(recenter_iterations, it)
                # only a significant move of the linearization point counts as a model
                # change, across which values of the bound are not comparable
                reset_baseline = step > 1e-2 * prior_scale
            end
            Γ = _gamma(gradM, Σδx, w, d_m, K)
        end

        # --- M-step: prior precisions, EM update eq. (alpha_update) ---
        α_em = 1 ./ (abs2.(μ) .+ diag(Σpost))
        if solver.options.learn_prior_flag
            α_new = solver.mackay_acceleration_flag ? _mackay_precisions(α, μ, Σpost, α_em) : α_em
            α = min.(α_new, α_cap)
        end

        # --- ARD pruning: remove sources whose every precision has diverged ---
        if solver.options.learn_prior_flag && solver.ard_prune_flag
            n_active = length(src_pos)
            keep = [j for j in 1:n_active if any(α[_source_columns(j, n_active, FD, iscomplex)] .< solver.ard_threshold)]
            if isempty(keep)
                @warn "automatic relevance determination switched off every source; keeping the most relevant one"
                keep = [argmin([minimum(α[_source_columns(j, n_active, FD, iscomplex)]) for j in 1:n_active])]
            end
            if length(keep) < n_active
                cols = _kept_columns(keep, n_active, FD, iscomplex)
                α = α[cols]
                α_em = α_em[cols]
                μ = μ[cols]
                Σpost = Σpost[cols, cols]
                logdetΣ = logdet(cholesky(Symmetric(Σpost)))
                M = M[:, cols]
                src_pos = src_pos[keep]
                if do_geometry
                    gradM = gradM[:, _kept_columns(keep, n_active, FD, false), :]
                    Γ = Γ[cols, cols]
                end
                K = length(α)
                Kh = length(src_pos) * FD
                reset_baseline = true
            end
        end

        # --- M-step: source positions χ, a few L-BFGS steps on eq. (chi_gradient) ---
        if solver.options.optimise_source_positions_flag && !isempty(src_pos)
            chi0 = Vector{Float64}(vcat(src_pos...))
            R0 = _expected_misfit(M, Γ, μ, Σpost, w, g)

            obj = chi -> _chi_misfit(chi, medium, bd_current, iscomplex, μ, Σpost, w, g, Dim;
                Γ = Γ, Σδx = do_geometry ? Σδx : nothing, d_m = d_m)

            opts = Optim.Options(iterations = solver.options.source_position_iters)
            use_analytic = solver.options.use_greens_gradient_analytical_flag && !iscomplex &&
                !do_geometry && (w isa AbstractVector) && _has_greens_gradient(bd.fieldtype, medium, Dim)

            res = if use_analytic
                grad! = (Gv, chi) -> _chi_misfit_gradient!(Gv, chi, medium, bd_current, μ, Σpost, w, g, Dim, FD)
                optimize(obj, grad!, chi0, LBFGS(), opts)
            else
                optimize(obj, chi0, LBFGS(), opts)
            end

            # accept only steps that increase the bound, i.e. decrease the expected misfit
            if Optim.minimum(res) < R0
                src_pos = _structured_positions(Optim.minimizer(res), Dim)
                M = _stacked_system_matrix(src_pos, medium, bd_current, iscomplex)
                if do_geometry
                    gradM = system_matrix_gradient(src_pos, medium, bd_current)
                    Γ = _gamma(gradM, Σδx, w, d_m, K)
                end
            end
        end

        # --- monitor the bound, eq. (elbo_model) ---
        R = _expected_misfit(M, Γ, μ, Σpost, w, g)
        F = _elbo(R, Nr, noise_logdet, α, μ, Σpost, logdetΣ;
            μδx = μδx, Σδx = Σδx, logdetΣδx = logdetΣδx, sx2 = sx2, m0 = m0)

        # MacKay acceleration carries no monotonicity guarantee: fall back to the EM update
        # whenever the bound fails to increase (Section on ARD of the theory document).
        if solver.options.learn_prior_flag && solver.mackay_acceleration_flag && !reset_baseline && F < F_prev
            α = min.(α_em, α_cap)
            F = _elbo(R, Nr, noise_logdet, α, μ, Σpost, logdetΣ;
                μδx = μδx, Σδx = Σδx, logdetΣδx = logdetΣδx, sx2 = sx2, m0 = m0)
        end

        push!(elbo_history, F)
        reset_baseline && push!(baseline_resets, it)

        converged = !reset_baseline && abs(F - F_prev) <= solver.elbo_tol * (1 + abs(F))
        F_prev = F
        converged && break
    end

    # --- final inference at the learned hyperparameters ---
    μ, Σpost, logdetΣ = _coefficient_posterior(M, Γ, α, w, g)
    R = _expected_misfit(M, Γ, μ, Σpost, w, g)
    misfit_ratio = R / Nr

    relative_boundary_error = norm(M * μ - g) / norm(g)

    # complex coefficients keep their posterior covariance over the stacked real degrees
    # of freedom [Re a; Im a], the convention of FundamentalSolution
    coefficients = if iscomplex
        Khh = length(src_pos) * FD
        μ[1:Khh] .+ im .* μ[(Khh + 1):end]
    else
        copy(μ)
    end

    fsol = FundamentalSolution(medium;
        positions = collect(src_pos),
        coefficients = coefficients,
        coefficients_covariance = Σpost,
        particular_solution = sim.particular_solution,
        relative_boundary_error = relative_boundary_error
    )

    # all the boundary information is returned as a BoundaryShape: when the geometry was
    # updated, the posterior of the boundary points is an MvNormal with the re-centered
    # boundary as mean and Σ_δx as covariance; otherwise the input shape is passed through
    boundary_shape = if do_geometry
        shape = bd_current.boundary_shape
        posterior_points = MvNormal(
            Vector{Float64}(vcat(mean_points(bd_current)...)),
            Symmetric(Σδx)
        )
        BoundaryShape(posterior_points, shape.normals, shape.interior_points)
    else
        bd.boundary_shape
    end

    return VariationalSolution(
        fsol, boundary_shape, α, elbo_history, baseline_resets, recenter_iterations, misfit_ratio
    )
end
