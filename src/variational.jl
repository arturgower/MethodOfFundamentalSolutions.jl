# Variational evidence maximization for the Method of Fundamental Solutions.
#
# Implements the two-phase algorithm of docs/theory/main-variational-evidence.tex:
#
#   Phase 1 (select) — section "Phase 1: pruning a dense set of candidate sources": sparse
#   Bayesian learning at fixed source positions, run constructively as the fast marginal
#   likelihood algorithm of Tipping & Faul: starting from no sources, the single best
#   action — add a column at its exact optimal precision eq. (alpha_star), re-estimate an
#   active precision, or delete a column — is applied until no action increases the
#   evidence. Optionally, each newly added source's position is refined by a few
#   hard-coded gradient-descent steps.
#
#   Phase 2 (move) — Algorithm 1: variational EM on the survivors, learning the prior
#   precisions {αᵢ} by the EM update eq. (alpha_update) (with the MacKay acceleration),
#   optionally the source positions χ, and optionally the boundary perturbation δx.
#
# The priors of the coefficients a and of the boundary points x are diagonal:
#   p(a)  = ∏ᵢ N(aᵢ | 0, 1/αᵢ)          (ARD prior)
#   p(δx) = N(0, Σ_x),  Σ_x diagonal     (taken from the covariance of the boundary points)
# The measurement noise Σ = diag(σ²) is known and fixed throughout.
#
# Complex-valued problems (e.g. acoustics) are handled by stacking real and imaginary
# parts: g̃ = [Re g; Im g], M̃ = [Re M  -Im M; Im M  Re M], ã = [Re a; Im a], so that the
# whole algorithm runs on a real linear-Gaussian model.
#
# Everything is phrased on the whitened augmented model Φ = [W M̄; V], ĝ = [W g; 0], where
# W is a left square root of the noise precision (Σ⁻¹ = WᵀW) and Γ = VᵀV is the extra
# coefficient precision from the boundary uncertainty, eq. (Gamma): then the coefficient
# posterior, the expected misfit and the per-column evidence factors are all plain
# least-squares expressions in Φ and ĝ. Since the basis is typically overcomplete, every
# posterior is evaluated in the data-sized Woodbury form of
# eqs. (posterior_mean_woodbury)-(posterior_covariance_woodbury): the K × K precision
# ΦᵀΦ + diag α is never formed, and the only Cholesky is of a P × P matrix, P being the
# number of rows of Φ.

"""
    VariationalBayesianSolver <: AbstractSolver

Solver implementing the two-phase variational evidence maximization of
docs/theory/main-variational-evidence.tex. By default (Phase 1) it starts from no sources
and adds one basis column at a time by the exact evidence criterion of Tipping and Faul's
fast marginal likelihood algorithm — the constructive form of the section "Phase 1: pruning
a dense set of candidate sources" — and then (Phase 2) runs variational EM on the selected
sources: it infers the posterior of the source coefficients `a`, learns the ARD prior
precisions `αᵢ`, optionally learns the source positions `χ`, and optionally infers the
boundary perturbation `δx`. The measurement noise is known.

The measurement noise is specified through the boundary data, not through the solver: give
the boundary `fields` as a vector of `MvNormal` (one per boundary point) or as a single
`MvNormal` over the flattened fields, whose covariance is the sensor noise. For
complex-valued problems each point must be a `2FD`-dimensional `MvNormal` over the stacked
real and imaginary parts `[Re; Im]` of its `FD` field components.

The options shared with [`BayesianSolver`](@ref) live in the field `options::SolverOptions`
(see [`SolverOptions`](@ref)); the keyword constructor accepts them as keywords directly.

# Keyword arguments
- `ard_prune_flag = true`: select the sources sparsely. With `learn_prior_flag` this runs
  the one-source-at-a-time selection of Phase 1, and afterwards keeps pruning any source
  whose precisions all exceed `ard_threshold` during Phase 2. When `false` every given
  source is kept.
- `learn_prior_flag = true`: learn the prior precisions (Phase 1 exactly, then the EM
  update `αᵢ = 1/(μᵢ² + Σᵢᵢ)` in Phase 2). When `false` the precisions stay at their
  initial values, and Phase 1 is skipped.
- `optimise_source_positions_flag = false`: optimise the source positions by hard-coded
  gradient descent with backtracking on the expected misfit — `source_position_iters`
  steps for each newly added source in Phase 1, and per source per iteration in Phase 2.
- `source_position_iters = 5`: gradient-descent steps per source-position update.
- `update_geometry_flag = false`: update the boundary factor `q(δx)` and re-center the
  boundary (E-step II). Requires the `boundary_points` to be an `MvNormal`; its diagonal
  covariance is the prior `Σ_x`. Not implemented for complex-valued problems.
- `priors`, `prior`, `prior_variance = 1.0`: the initial prior variances `1/αᵢ` of the
  coefficients, used only when Phase 1 is skipped (Phase 1 sets each precision by its
  exact evidence maximizer). A vector of `MvNormal` (one per source), a single `MvNormal`,
  a scalar shared by all coefficients, or a vector with one entry per coefficient (for
  complex problems either per complex coefficient or per real degree of freedom `[Re; Im]`).
- `ard_threshold = 1e8`: precision above which a coefficient counts as switched off.
- `mackay_acceleration_flag = true`: use the MacKay/Tipping fixed-point update
  `αᵢ = γᵢ/μᵢ²` in Phase 2 as an acceleration, falling back to the EM update whenever the
  bound decreases.
- `elbo_tol = 1e-8`: relative tolerance on the increase of the evidence lower bound.
- `max_iters = 200`: maximum number of variational EM iterations of Phase 2, each of which
  updates the whole model.
- `max_select_iters = -1`: maximum number of actions of the Phase 1 selection, which adds,
  re-estimates or deletes ONE coefficient per action and so needs a budget of a different
  order from `max_iters` — typically a few times the number of data rows. Negative (the
  default) means automatic: four times the number of candidate coefficients. Selection stops
  by itself as soon as no action increases the evidence, so a generous budget costs nothing,
  whereas an under-sized one silently truncates the selection.
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
        elbo_tol::Real = 1e-8,
        max_iters::Int = 200,
        max_select_iters::Int = -1,
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
        update_geometry_flag = update_geometry_flag,
        learn_prior_flag = learn_prior_flag,
        max_iters = max_iters,
        max_select_iters = max_select_iters,
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
- `elbo_history`: the evidence lower bound F after each iteration, of Phase 1 (source
  selection) followed by Phase 2 (variational EM).
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
overcomplete set of candidates for a [`VariationalBayesianSolver`](@ref), which selects the
useful sources one at a time.
"""
function grid_source_positions(bd::BoundaryData{F, 2}; n::Int = 15, scale::Real = 2.0, clearance::Real = 1.0) where F
    pts = mean_points(bd)
    len = length(pts)

    xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
    centre = SVector((minimum(xs) + maximum(xs)) / 2, (minimum(ys) + maximum(ys)) / 2)
    halfwidth = SVector(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) ./ 2

    spacing = _boundary_spacing(pts)

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

# precisions are capped here rather than at Inf, so that a switched-off coefficient still
# has a well-defined (tiny) prior variance in Phase 2
const ALPHA_CAP = 1e12

# Average spacing between neighbouring boundary points, as in source_positions.
function _boundary_spacing(pts)
    len = length(pts)
    sampled_rng = LinRange(1, len, min(6, len)) .|> round .|> Int
    return mean(map(pts[sampled_rng]) do p
        dists = [norm(p - q) for q in pts]
        idx = sortperm(dists)[2:min(3, len)]
        mean(dists[idx])
    end)
end

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

# Cholesky of a symmetric positive definite matrix, retried with a small diagonal jitter when
# round-off makes the factorization fail.
function _safe_cholesky(A::AbstractMatrix)
    S = Symmetric(Matrix(A))
    C = cholesky(S; check = false)
    issuccess(C) && return C
    jitter = 1e-12 * max(1.0, maximum(abs, diag(S)))
    return cholesky(Symmetric(S + jitter * I))
end

# A left square root W of the measurement-noise precision, so that Σ⁻¹ = WᵀW: the whitened
# model matrix W M is what every data-sized posterior is built from. The noise precision is
# stored either as a vector of reciprocal variances w = 1/σ² (independent sensor noise, the
# fast path) or as a full precision matrix w = Σ⁻¹ (correlated noise).
_whitener(w::AbstractVector) = Diagonal(sqrt.(w))
_whitener(w::AbstractMatrix) = Matrix(transpose(_safe_cholesky(w).L))

# The augmented design Φ = [W M̄; V] of the whitened model: with Γ = VᵀV the extra
# coefficient precision from the boundary uncertainty, ΦᵀΦ = M̄ᵀ Σ⁻¹ M̄ + Γ.
_phi(Mw::AbstractMatrix, ::Nothing) = Mw
_phi(Mw::AbstractMatrix, Vγ::AbstractMatrix) = vcat(Mw, Vγ)

# E-step I, eq. (vb_coefficient_update), on the whitened augmented model:
# Σ_post = (ΦᵀΦ + diag α)⁻¹ and μ_post = Σ_post Φᵀ ĝ, evaluated in the data-sized Woodbury
# form of eq. (posterior_covariance_woodbury). With the diagonal prior A = diag α,
#   Σ_post = A⁻¹ - A⁻¹Φᵀ (I + Φ A⁻¹ Φᵀ)⁻¹ Φ A⁻¹,   log|Σ_post| = -log|I + Φ A⁻¹ Φᵀ| - Σᵢ log αᵢ
# (the log-determinant by the matrix determinant lemma). The only factorization is of the
# P × P matrix I + Φ A⁻¹ Φᵀ, P being the number of rows of Φ.
function _coefficient_posterior(Φ::AbstractMatrix, α::AbstractVector, ĝ::AbstractVector)
    αinv = 1 ./ α
    ΦA = Φ .* transpose(αinv)                  # Φ A⁻¹, P × K
    C = _safe_cholesky(ΦA * transpose(Φ) + I)  # I + Φ A⁻¹ Φᵀ, P × P
    Σpost = Matrix(Symmetric(Diagonal(αinv) - transpose(ΦA) * (C \ ΦA)))
    logdetΣ = -logdet(C) - sum(log, α)
    μ = Σpost * (transpose(Φ) * ĝ)
    return μ, Σpost, logdetΣ
end

# The noise-weighted expected misfit R of eq. (R) on the whitened augmented model:
# R = ‖ĝ - Φμ‖² + tr(Φ Σ_post Φᵀ). Through the rows [W M̄; V] of Φ this expands into the
# misfit of the mean solution, the extra misfit from the coefficient uncertainty, and the
# boundary-uncertainty terms μᵀΓμ + tr(Γ Σ_post).
function _expected_misfit(Φ::AbstractMatrix, μ::AbstractVector, Σpost::AbstractMatrix, ĝ::AbstractVector)
    r = ĝ - Φ * μ
    return dot(r, r) + sum((Φ * Σpost) .* Φ)
end

# The evidence lower bound F of eq. (elbo_model), with the Gaussian KL of eq. (gaussian_kl).
# The prior of the boundary perturbation stays anchored at the original nominal boundary:
# relative to the current linearization point it is N(m0, Σx) with m0 = μ_x0 - x_lin.
function _elbo(R::Real, Nr::Int, noise_logdet::Real, α::AbstractVector, μ::AbstractVector,
        Σpost::AbstractMatrix, logdetΣ::Real;
        μδx = nothing, Σδx = nothing,
        logdetΣδx::Real = 0.0,
        sx2 = nothing, m0 = nothing
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

# The factor V of the extra coefficient precision Γ = VᵀV = Σ_kl Σδx^kl M_kᵀ Σ⁻¹ M_l of
# eq. (Gamma). Row r of M depends only on the boundary point of its own sensor, so only the
# Dim × Dim diagonal blocks Σδx^(i) of Σδx contribute and, with the K × Dim gradient block
# G_r = [∂M_r/∂x_1 … ∂M_r/∂x_Dim] of data row r and its sensor i(r),
#   Γ = Σ_r w_r G_r Σδx^(i(r)) G_rᵀ.
# Each Σδx^(i) is a covariance, so factoring it as L_i L_iᵀ gives Γ = VᵀV with the Dim rows
# √w_r L_i(r)ᵀ G_rᵀ of V per data row r. Γ is only ever needed through this factor, so it is
# never formed: that keeps the coefficient posterior in its data-sized form, and costs
# O(Nr Dim² K) instead of the O(Nr Dim² K²) of accumulating the K × K matrix.
function _gamma_factor(gradM::AbstractArray, Σδx::AbstractMatrix, w::AbstractVector, d_m::Int)
    N, K = size(gradM, 1), size(gradM, 2)
    Dim = size(gradM, 3)
    Vγ = zeros(N * Dim, K)
    for i in 1:(N ÷ d_m)
        ks = ((i - 1) * Dim + 1):(i * Dim)
        Lt = transpose(_safe_cholesky(Σδx[ks, ks]).L)
        for r in (((i - 1) * d_m + 1):(i * d_m))
            Vγ[((r - 1) * Dim + 1):(r * Dim), :] = sqrt(w[r]) .* (Lt * transpose(Matrix(gradM[r, :, :])))
        end
    end
    return Vγ
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
            ALPHA_CAP
        else
            γ / μ2
        end
    end
end

# ------------------------------------------------------------------------------------
# Source positions: hard-coded gradient descent on the expected misfit
# ------------------------------------------------------------------------------------

# The whitened columns that one source at position χ contributes to the augmented design
# Φ = [W M; V]: with boundary uncertainty (Σδx given) the V rows come from the Γ factor of
# the source's own gradient columns.
function _source_basis(χ::SVector, medium, bd, iscomplex::Bool, W, Σδx, w, d_m::Int)
    Φj = W * _stacked_system_matrix([χ], medium, bd, iscomplex)
    if Σδx !== nothing
        Φj = vcat(Φj, _gamma_factor(system_matrix_gradient([χ], medium, bd), Σδx, w, d_m))
    end
    return Φj
end

# One gradient-descent step with backtracking on the expected misfit R of eq. (R) as a
# function of the position χ of a single source, holding the coefficient posterior q(a) and
# every other source fixed. The bound F depends on χ only through R, so an accepted step
# increases F (a generalized EM step). `basis(χ)` returns the columns of Φ that move with
# the source (the columns `cols`), and the gradient is by central finite differences on
# those columns alone. The trial step starts at a quarter of the source's distance to the
# boundary — so no step can ever cross the boundary — and is halved until R decreases;
# steps that bring the source closer than `clearance` to a boundary point are rejected,
# since such near-singular columns just chase individual sensors.
function _descend_step(χ::SVector{Dim, Float64}, basis, Φ::AbstractMatrix,
        cols::AbstractVector{Int}, μ::AbstractVector, Σpost::AbstractMatrix,
        ĝ::AbstractVector, pts, clearance::Real) where Dim

    others = setdiff(axes(Φ, 2), cols)
    r = ĝ - Φ[:, others] * μ[others]              # residual without this source
    μj = μ[cols]
    Σjj = Σpost[cols, cols]
    B = Φ[:, others] * Σpost[others, cols]        # posterior coupling to the other sources
    # R(χ) up to terms independent of χ: ‖r - Φⱼμⱼ‖² + 2 tr(ΦⱼᵀB) + tr(ΦⱼᵀΦⱼ Σⱼⱼ)
    R(Φj) = sum(abs2, r - Φj * μj) + 2 * sum(Φj .* B) + sum((transpose(Φj) * Φj) .* Σjj)
    boundary_distance(x) = minimum(norm(x - p) for p in pts)

    Rχ = R(basis(χ))
    dist = boundary_distance(χ)
    h = 1e-6 * dist
    G = SVector{Dim}(ntuple(Dim) do d
        χp = Base.setindex(χ, χ[d] + h, d)
        χm = Base.setindex(χ, χ[d] - h, d)
        (R(basis(χp)) - R(basis(χm))) / (2h)
    end)
    Gn = norm(G)
    Gn > 0 || return χ, false

    t = dist / (4 * Gn)
    for _ in 1:8
        χt = χ - t * G
        if R(basis(χt)) < Rχ && boundary_distance(χt) >= clearance
            return χt, true
        end
        t /= 2
    end
    return χ, false
end

# `source_position_iters` rounds of: one descent step on the position of the source with
# columns `cols` of Φ, then a re-solve of the coefficient posterior q(a) — alternating the
# two lets the source travel far while its amplitude follows, where descent at fixed q(a)
# stalls after a short move. Both parts increase the bound F, so the alternation does too.
# The moved columns are written into Φ; returns the position and the refreshed posterior.
function _optimise_source_position!(Φ::AbstractMatrix, cols::AbstractVector{Int}, basis,
        χ0::SVector, α::AbstractVector, ĝ::AbstractVector, iters::Int, pts, clearance::Real)

    χ = χ0
    μ, Σpost, logdetΣ = _coefficient_posterior(Φ, α, ĝ)
    for _ in 1:iters
        χnew, moved = _descend_step(χ, basis, Φ, cols, μ, Σpost, ĝ, pts, clearance)
        moved || break
        χ = χnew
        Φ[:, cols] = basis(χ)
        μ, Σpost, logdetΣ = _coefficient_posterior(Φ, α, ĝ)
    end
    return χ, μ, Σpost, logdetΣ, χ != χ0
end

# ------------------------------------------------------------------------------------
# Phase 1: add one source at a time (fast marginal likelihood, Tipping & Faul)
# ------------------------------------------------------------------------------------

# The log evidence as a function of a single precision α, eq. (sq_factors), up to terms
# independent of α; the limit ℓ(∞) = 0 is the column switched off.
_log_evidence_1(α::Real, s::Real, q::Real) = isinf(α) ? 0.0 : (log(α / (α + s)) + q^2 / (α + s)) / 2

# Sparse Bayesian learning at fixed candidate positions, run constructively: starting from
# no sources, apply the single best action — add a column at its optimal precision
# eq. (alpha_star), re-estimate the precision of an active column, or delete one — until no
# action increases the evidence. Each action is the exact coordinate maximizer of the
# evidence, so the bound increases monotonically. After each addition the position of the
# added source is optionally refined by a few gradient-descent steps.
# Mutates src_pos (positions may move); returns the indices of the surviving sources, the
# precisions of their columns, and the evidence-bound history.
function _add_sources!(src_pos, medium, bd, iscomplex::Bool, W, w, d_m::Int, Σδx,
        ĝ::AbstractVector, Nr::Int, noise_logdet::Real, FD::Int, solver)

    n_src = length(src_pos)
    K = (iscomplex ? 2 : 1) * n_src * FD

    Vγ = Σδx === nothing ? nothing :
        _gamma_factor(system_matrix_gradient(src_pos, medium, bd), Σδx, w, d_m)
    Φ = _phi(W * _stacked_system_matrix(src_pos, medium, bd, iscomplex), Vγ)

    α = fill(Inf, K)
    elbo = Float64[]
    F_prev = 0.0
    pts = mean_points(bd)
    clearance = _boundary_spacing(pts) / 2

    # one action per iteration, so this phase needs its own budget: a negative
    # `max_select_iters` asks for the automatic four-per-candidate-coefficient safety net
    budget = solver.options.max_select_iters
    budget < 0 && (budget = 4K)

    for _ in 1:budget
        active = findall(isfinite, α)

        # the evidence covariance C = I + Φ_A A⁻¹ Φ_Aᵀ of the active model, and the factors
        # S = φᵢᵀC⁻¹φᵢ (sparsity) and Q = φᵢᵀC⁻¹ĝ (quality) of every candidate column
        ΦA = Φ[:, active]
        C = _safe_cholesky(ΦA * Diagonal(1 ./ α[active]) * transpose(ΦA) + I)
        S = vec(sum(Φ .* (C \ Φ); dims = 1))
        Q = transpose(Φ) * (C \ ĝ)

        # the single best action: for each column the leave-one-out factors s, q of
        # eq. (sq_factors), its optimal precision eq. (alpha_star), and the evidence gain
        best_i, best_Δ, best_α = 0, solver.elbo_tol * (1 + abs(F_prev)), Inf
        for i in 1:K
            s, q = if isinf(α[i])
                S[i], Q[i]
            else            # S, Q include column i itself; remove its own contribution
                den = max(α[i] - S[i], 1e-12 * α[i])
                (α[i] * S[i] / den, α[i] * Q[i] / den)
            end
            α_i = q^2 > s ? s^2 / (q^2 - s) : Inf
            Δ = _log_evidence_1(α_i, s, q) - _log_evidence_1(α[i], s, q)
            if Δ > best_Δ
                best_i, best_Δ, best_α = i, Δ, α_i
            end
        end
        best_i == 0 && break        # no action improves the evidence: converged

        added = isinf(α[best_i])
        α[best_i] = best_α

        # refine the position of the newly added source by a few gradient-descent steps
        if added && solver.options.optimise_source_positions_flag && solver.options.source_position_iters > 0
            j = mod(best_i - 1, n_src * FD) ÷ FD + 1
            cols = _source_columns(j, n_src, FD, iscomplex)
            sub = findall(c -> isfinite(α[c]), cols)      # this source's active columns
            keep = findall(isfinite, α)
            basis = χ -> _source_basis(χ, medium, bd, iscomplex, W, Σδx, w, d_m)
            χ, _, _, _, moved = _optimise_source_position!(
                Φ[:, keep], findall(in(cols[sub]), keep), χc -> basis(χc)[:, sub],
                src_pos[j], α[keep], ĝ, solver.options.source_position_iters, pts, clearance
            )
            if moved
                src_pos[j] = χ
                Φ[:, cols] = basis(χ)
            end
        end

        # monitor the bound: with q(a) the exact posterior it equals the log evidence
        keep = findall(isfinite, α)
        μ, Σpost, logdetΣ = _coefficient_posterior(Φ[:, keep], α[keep], ĝ)
        R = _expected_misfit(Φ[:, keep], μ, Σpost, ĝ)
        F = _elbo(R, Nr, noise_logdet, α[keep], μ, Σpost, logdetΣ)
        push!(elbo, F)
        F_prev = F
    end

    if !any(isfinite, α)
        @warn "no candidate source explains the data beyond the noise; keeping the best one"
        i = argmax(abs2.(transpose(Φ) * ĝ) ./ vec(sum(abs2, Φ; dims = 1)))
        α[i] = ALPHA_CAP
    end

    keep_src = [j for j in 1:n_src if any(isfinite, α[_source_columns(j, n_src, FD, iscomplex)])]
    cols = _kept_columns(keep_src, n_src, FD, iscomplex)
    return keep_src, min.(α[cols], ALPHA_CAP), elbo
end

# ------------------------------------------------------------------------------------
# solve: Phase 1 (select the sources), then Phase 2 (variational EM, Algorithm 1)
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
    W = _whitener(w)

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

    Σδx = do_geometry ? Matrix(Diagonal(sx2)) : nothing            # initialize q(δx) at the prior
    logdetΣδx = do_geometry ? sum(log.(sx2)) : 0.0
    μδx = do_geometry ? zeros(length(sx2)) : nothing
    # prior mean of δx relative to the current linearization point: the prior of the
    # boundary stays anchored at the original nominal boundary μ_x0 across re-centerings
    x_nominal = do_geometry ? Vector{Float64}(vcat(mean_points(bd)...)) : nothing
    m0 = do_geometry ? zeros(length(sx2)) : nothing

    # whitened data of the augmented model: the Γ ridge rows carry zero data
    ĝ = do_geometry ? vcat(W * g, zeros(Nr * Dim)) : W * g

    elbo_history = Float64[]

    # --- Phase 1: add one source at a time by the exact evidence criterion ---
    if solver.ard_prune_flag && solver.options.learn_prior_flag
        keep_src, α, elbo1 = _add_sources!(src_pos, medium, bd_current, iscomplex,
            W, w, d_m, Σδx, ĝ, Nr, noise_logdet, FD, solver)
        src_pos = src_pos[keep_src]
        append!(elbo_history, elbo1)
    else
        n_src = length(src_pos)
        Kh = n_src * FD
        α = _initial_precisions(solver.prior_variance, iscomplex ? 2Kh : Kh, Kh, iscomplex)
    end

    # --- Phase 2: variational EM on the selected sources ---
    M = _stacked_system_matrix(src_pos, medium, bd_current, iscomplex)
    gradM = do_geometry ? system_matrix_gradient(src_pos, medium, bd_current) : nothing
    Vγ = do_geometry ? _gamma_factor(gradM, Σδx, w, d_m) : nothing
    Φ = _phi(W * M, Vγ)

    baseline_resets = Int[]
    recenter_iterations = Int[]
    # continue monitoring from the last Phase-1 value: the bound is comparable across the
    # transition, so the MacKay fallback also guards the first Phase-2 iteration
    F_prev = isempty(elbo_history) ? -Inf : elbo_history[end]
    μ = zeros(length(α))
    Σpost = Matrix{Float64}(I, length(α), length(α))
    logdetΣ = 0.0

    for _ in 1:solver.options.max_iters
        reset_baseline = false

        # --- E-step I: update q(a), eq. (vb_coefficient_update) ---
        μ, Σpost, logdetΣ = _coefficient_posterior(Φ, α, ĝ)

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
                push!(recenter_iterations, length(elbo_history) + 1)
                # only a significant move of the linearization point counts as a model
                # change, across which values of the bound are not comparable
                reset_baseline = step > 1e-2 * prior_scale
            end
            Vγ = _gamma_factor(gradM, Σδx, w, d_m)
            Φ = _phi(W * M, Vγ)
        end

        # --- M-step: prior precisions, EM update eq. (alpha_update) ---
        α_em = 1 ./ (abs2.(μ) .+ diag(Σpost))
        if solver.options.learn_prior_flag
            α_new = solver.mackay_acceleration_flag ? _mackay_precisions(α, μ, Σpost, α_em) : α_em
            α = min.(α_new, ALPHA_CAP)
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
                Φ = Φ[:, cols]
                src_pos = src_pos[keep]
                if do_geometry
                    gradM = gradM[:, _kept_columns(keep, n_active, FD, false), :]
                    Vγ = Vγ[:, cols]
                end
                reset_baseline = true
            end
        end

        # --- M-step: source positions χ, a few gradient-descent steps per source ---
        if solver.options.optimise_source_positions_flag && solver.options.source_position_iters > 0
            pts = mean_points(bd_current)
            clearance = _boundary_spacing(pts) / 2
            moved_any = false
            for j in eachindex(src_pos)
                cols = _source_columns(j, length(src_pos), FD, iscomplex)
                basis = χ -> _source_basis(χ, medium, bd_current, iscomplex, W, Σδx, w, d_m)
                χ, μ, Σpost, logdetΣ, moved = _optimise_source_position!(Φ, cols, basis,
                    src_pos[j], α, ĝ, solver.options.source_position_iters, pts, clearance)
                if moved
                    src_pos[j] = χ
                    moved_any = true
                end
            end
            if moved_any
                M = _stacked_system_matrix(src_pos, medium, bd_current, iscomplex)
                do_geometry && (gradM = system_matrix_gradient(src_pos, medium, bd_current))
            end
        end

        # --- monitor the bound, eq. (elbo_model) ---
        R = _expected_misfit(Φ, μ, Σpost, ĝ)
        F = _elbo(R, Nr, noise_logdet, α, μ, Σpost, logdetΣ;
            μδx = μδx, Σδx = Σδx, logdetΣδx = logdetΣδx, sx2 = sx2, m0 = m0)

        # MacKay acceleration carries no monotonicity guarantee: fall back to the EM update
        # whenever the bound fails to increase (Section on ARD of the theory document).
        if solver.options.learn_prior_flag && solver.mackay_acceleration_flag && !reset_baseline && F < F_prev
            α = min.(α_em, ALPHA_CAP)
            F = _elbo(R, Nr, noise_logdet, α, μ, Σpost, logdetΣ;
                μδx = μδx, Σδx = Σδx, logdetΣδx = logdetΣδx, sx2 = sx2, m0 = m0)
        end

        push!(elbo_history, F)
        reset_baseline && push!(baseline_resets, length(elbo_history))

        converged = !reset_baseline && abs(F - F_prev) <= solver.elbo_tol * (1 + abs(F))
        F_prev = F
        converged && break
    end

    # --- final inference at the learned hyperparameters ---
    μ, Σpost, logdetΣ = _coefficient_posterior(Φ, α, ĝ)
    R = _expected_misfit(Φ, μ, Σpost, ĝ)
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
