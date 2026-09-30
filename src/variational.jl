# Variational evidence maximization for the Method of Fundamental Solutions.
#
# Implements the two-phase algorithm of docs/theory/main-variational-evidence.tex:
#
#   Phase 1 (select) — section "Phase 1: pruning a dense set of candidate sources": sparse
#   Bayesian learning at FIXED source positions, run constructively as the fast marginal
#   likelihood algorithm of Tipping & Faul. Starting from no sources, the single best action —
#   add a column at its exact optimal precision eq. (alpha_star), re-estimate an active
#   precision, or delete a column — is applied until no action increases the evidence.
#
#   Phase 2 (move) — Algorithm 1: variational EM on the survivors, learning the prior
#   precisions {αᵢ} by the EM update eq. (alpha_update) (with the MacKay acceleration),
#   optionally the source positions χ by gradient descent, and optionally the boundary
#   perturbation δx.
#
# The priors of the coefficients a and of the boundary points x are diagonal:
#   p(a)  = ∏ᵢ N(aᵢ | 0, 1/αᵢ)          (ARD prior)
#   p(δx) = N(0, Σ_x),  Σ_x diagonal     (taken from the covariance of the boundary points)
# The measurement noise Σ = diag(σ²) is known and fixed throughout.
#
# Complex-valued problems (e.g. acoustics) are handled by stacking real and imaginary parts:
# g̃ = [Re g; Im g], M̃ = [Re M  -Im M; Im M  Re M], ã = [Re a; Im a], so that the whole
# algorithm runs on a real linear-Gaussian model — the `WorkingModel` below.
#
# Everything is phrased on the whitened augmented design Φ = [W M̄; V], with data ĝ = [W g; 0],
# where W is a left square root of the noise precision (Σ⁻¹ = WᵀW) and Γ = VᵀV is the extra
# coefficient precision from the boundary uncertainty, eq. (Gamma): the coefficient posterior,
# the expected misfit and the per-column evidence factors are then plain least-squares
# expressions in Φ and ĝ. Since the basis is typically overcomplete, Phase 2 evaluates every
# posterior in the data-sized Woodbury form of eqs. (posterior_mean_woodbury)-
# (posterior_covariance_woodbury): the K × K precision ΦᵀΦ + diag α is never formed, and the
# only Cholesky is of a P × P matrix, P being the number of rows of Φ. Phase 1 goes the other
# way round — its active set is far smaller than P — and works with the k × k precision of the
# active columns alone; see `_select_columns`.

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

# Keyword arguments
- `select_sources_flag = ard_prune_flag`: run the Phase 1 selection, which starts from an
  empty basis and adds one candidate coefficient at a time by the exact evidence criterion.
  Set to `false` (with `ard_prune_flag` left on) to enter Algorithm 1 directly on a set of
  sources that Phase 1 has already chosen — the two-stage form of the theory, where Phase 2
  is `Require`d to start from the survivors of Phase 1 and their precisions (pass those back
  as `prior_variance = 1 ./ vsol.prior_precisions`).
- `ard_prune_flag = true`: delete any source whose precisions all exceed `ard_threshold`
  during Phase 2 (Algorithm 1, line 8), and, unless `select_sources_flag` says otherwise,
  run the Phase 1 selection first. When `false` every given source is kept.
- `learn_prior_flag = true`: learn the prior precisions (Phase 1 exactly, then the EM
  update `αᵢ = 1/(μᵢ² + Σᵢᵢ)` in Phase 2). When `false` the precisions stay at their
  initial values, and Phase 1 is skipped.
- `optimise_source_positions_flag = false`: optimise the source positions in Phase 2 by
  gradient descent with backtracking on the expected misfit, `source_position_iters` steps
  per source per iteration. The positions are fixed throughout Phase 1, as the theory
  prescribes: all its guarantees hold only for a fixed basis.
- `source_position_iters = 5`: gradient-descent steps per source-position update.
- `boundary_ridge_flag = true`: when the `boundary_points` carry a covariance `Σ_x`, add
  the ridge `Γ` of eq. (Gamma_x) to the coefficient precision — the closed-form evidence
  bound of the section "A closed-form evidence bound: boundary variance as pseudo-data".
  It tempers the coefficients in the directions in which the boundary trace is most
  sensitive to where the sensors actually are, without inferring the boundary. This is what
  Phase 1 always uses; with `update_geometry_flag = false` it is also all that Phase 2 uses,
  which is the regime to pick when the sensor positions are uncertain but not identifiable
  from the data. Set to `false` to ignore `Σ_x` entirely. Requires independent (diagonal)
  sensor noise.
- `update_geometry_flag = false`: additionally update the boundary factor `q(δx)` and
  re-center the boundary (E-step II), i.e. learn the boundary rather than only allow for
  it. Requires the `boundary_points` to be an `MvNormal` (or a vector of them, one per
  point), whose diagonal covariance is `Σ_x`.
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
    prior_variance::Vector{Float64}
    select_sources_flag::Bool
    ard_prune_flag::Bool
    ard_threshold::Float64
    learn_prior_flag::Bool
    mackay_acceleration_flag::Bool
    optimise_source_positions_flag::Bool
    source_position_iters::Int
    boundary_ridge_flag::Bool
    update_geometry_flag::Bool
    elbo_tol::Float64
    max_iters::Int
    max_select_iters::Int
end

function VariationalBayesianSolver(;
        priors::AbstractVector{<:AbstractMvNormal} = MvNormal[],
        prior::Union{Nothing, ContinuousMultivariateDistribution} = nothing,
        prior_variance::Union{Real, AbstractVector{<:Real}} = 1.0,
        ard_prune_flag::Bool = true,
        select_sources_flag::Bool = ard_prune_flag,
        ard_threshold::Real = 1e8,
        learn_prior_flag::Bool = true,
        mackay_acceleration_flag::Bool = true,
        optimise_source_positions_flag::Bool = false,
        source_position_iters::Int = 5,
        boundary_ridge_flag::Bool = true,
        update_geometry_flag::Bool = false,
        elbo_tol::Real = 1e-8,
        max_iters::Int = 200,
        max_select_iters::Int = -1
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

    return VariationalBayesianSolver(
        pvs, select_sources_flag, ard_prune_flag, Float64(ard_threshold),
        learn_prior_flag, mackay_acceleration_flag,
        optimise_source_positions_flag, source_position_iters,
        boundary_ridge_flag, update_geometry_flag,
        Float64(elbo_tol), max_iters, max_select_iters
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

# ------------------------------------------------------------------------------------
# The working model and the two Gaussian factors of the variational posterior
# ------------------------------------------------------------------------------------

# precisions are capped here rather than at Inf, so that a switched-off coefficient still
# has a well-defined (tiny) prior variance in Phase 2
const ALPHA_CAP = 1e12

"""
    WorkingModel

Everything both phases share about the real (possibly `[Re; Im]`-stacked) working model: the
physics, the data with the particular solution subtracted, and the known measurement noise.
Only the boundary points, the source positions and `q(δx)` change during a solve, so those
are passed around separately.

# Fields
- `medium`, `field_dim`, `iscomplex`: the physics, and how it maps onto the real model.
- `g`: the data of the real model, i.e. the boundary fields minus the particular solution,
  `[Re g; Im g]`-stacked for complex physics.
- `ĝ`: the data of the whitened augmented model, `[W g; 0]` — the Γ rows carry no data.
- `W`, `w`, `noise_logdet`: a left square root of the noise precision (Σ⁻¹ = WᵀW), the
  precision Σ⁻¹ itself (a vector of reciprocal variances for independent sensor noise, a
  matrix for correlated noise), and log|Σ|.
- `row_sensor`: the boundary sensor each row of the real model belongs to. Real problems
  stack `d_m` rows per sensor in sensor order; complex ones stack that block twice, once for
  the real parts and once for the imaginary parts. Every quantity that groups rows by sensor
  — the Γ factor and E-step II — goes through this map rather than through row arithmetic.
"""
struct WorkingModel{Dim, P <: PhysicalMedium{Dim}, WT <: AbstractMatrix, NT}
    medium::P
    field_dim::Int
    iscomplex::Bool
    g::Vector{Float64}
    ĝ::Vector{Float64}
    W::WT
    w::NT
    noise_logdet::Float64
    row_sensor::Vector{Int}
end

"""
    BoundaryFactor

The Gaussian factor q(δx) = N(`mean`, `cov`) over the boundary perturbation, present
whenever the boundary points carry a covariance Σ_x that is not being ignored. It is frozen
at its prior N(0, Σ_x) unless E-step II updates it (`update_geometry_flag`), and in that
frozen state its contribution to the bound is identically zero.

`prior_mean` is the prior mean of δx *relative to the current linearization point*, i.e. the
nominal boundary seen from the re-centered one: the prior stays anchored at the original
nominal boundary, so directions the data cannot inform are pulled back to it instead of
drifting with each re-linearization. `prior_var` is the diagonal of Σ_x.
"""
mutable struct BoundaryFactor
    mean::Vector{Float64}
    cov::Matrix{Float64}
    logdet_cov::Float64
    prior_mean::Vector{Float64}
    prior_var::Vector{Float64}
end

BoundaryFactor(sx2::Vector{Float64}) = BoundaryFactor(
    zeros(length(sx2)), Matrix(Diagonal(sx2)), sum(log, sx2), zeros(length(sx2)), sx2
)

# The coefficient factor q(a) = N(mean, cov), with the log-determinant of the covariance the
# bound needs (which the Woodbury form of `_coefficient_posterior` returns for free).
struct CoefficientPosterior
    mean::Vector{Float64}
    cov::Matrix{Float64}
    logdet_cov::Float64
end

_empty_posterior() = CoefficientPosterior(Float64[], zeros(0, 0), 0.0)

# The posterior restricted to the coefficients `cols`, after ARD switched off the others.
_restrict(post::CoefficientPosterior, cols::AbstractVector{Int}) =
    CoefficientPosterior(post.mean[cols], post.cov[cols, cols],
        logdet(_safe_cholesky(post.cov[cols, cols])))

# ------------------------------------------------------------------------------------
# Assembling the working model from the boundary data
# ------------------------------------------------------------------------------------

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

# Build the working model of `sim`: subtract the particular solution, stack [Re; Im] for
# complex physics, and whiten by the known measurement noise. `use_ridge` says whether the
# augmented design carries the Γ rows, which the whitened data pads with zeros.
function _working_model(sim::Simulation, use_ridge::Bool)
    medium = sim.medium
    bd = sim.boundary_data
    FD = field_dimension(medium)
    Dim = spatial_dimension(medium)

    M0 = system_matrix(sim.source_positions, medium, bd)
    iscomplex = eltype(M0) <: Complex

    # A fundamental solution is singular at its own source, so a source sitting on a boundary
    # point makes a column infinite — and the failure would otherwise surface much later, as an
    # unreadable LAPACK error inside one of the factorizations.
    if !all(isfinite, M0)
        bad = findall(!isfinite, M0)
        j = min(first(bad)[2] ÷ max(1, FD) + 1, length(sim.source_positions))
        throw(ArgumentError(
            "the system matrix is not finite in $(length(bad)) of its entries: source " *
            "$j at $(sim.source_positions[j]) lies on (or unusably close to) a boundary " *
            "point. Keep every source clear of the boundary."
        ))
    end

    g0, noise = _data_and_noise(bd, FD, iscomplex)
    g0 = g0 - vcat(field(medium, bd, sim.particular_solution)...)
    g = iscomplex ? Vector{Float64}(vcat(real.(g0), imag.(g0))) : Vector{Float64}(g0)

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
    ĝ = use_ridge ? vcat(W * g, zeros(length(g) * Dim)) : W * g

    n_sensors = length(mean_points(bd))
    d_m = size(M0, 1) ÷ n_sensors
    row_sensor = repeat(repeat(1:n_sensors, inner = d_m), iscomplex ? 2 : 1)

    return WorkingModel{Dim, typeof(medium), typeof(W), typeof(w)}(
        medium, FD, iscomplex, g, ĝ, W, w, noise_logdet, row_sensor
    )
end

# The diagonal Σ_x of the boundary prior, or `nothing` when the boundary is deterministic or
# its covariance is to be ignored. Its presence is what turns the ridge Γ on.
function _boundary_variances(sim::Simulation, solver::VariationalBayesianSolver)
    Σx = cov(sim.boundary_data.boundary_shape.boundary_points)
    (Σx isa UniformScaling) && return nothing
    (solver.update_geometry_flag || solver.boundary_ridge_flag) || return nothing
    return Vector{Float64}(diag(Σx))
end

# The model matrix of the real working model: for complex physics the real/imaginary parts
# are stacked so that g̃ = M̃ ã with ã = [Re a; Im a].
function _stacked_system_matrix(model::WorkingModel, bd, src_pos)
    M0 = system_matrix(src_pos, model.medium, bd)
    model.iscomplex || return Matrix{Float64}(M0)
    Mr = real.(M0); Mi = imag.(M0)
    return [Mr -Mi; Mi Mr]
end

# ∂M̃/∂x of the same real working model, as an (Nr × K × Dim) array: the derivative of every
# row of M̃ with respect to the position of its own boundary sensor, stacked exactly like
# `_stacked_system_matrix` so that rows and columns keep the [Re; Im] ordering.
function _stacked_gradient(model::WorkingModel, bd, src_pos)
    G = system_matrix_gradient(src_pos, model.medium, bd)
    model.iscomplex || return Array{Float64}(G)

    N, K, Dim = size(G)
    S = zeros(Float64, 2N, 2K, Dim)
    for d in 1:Dim
        Gr = real.(view(G, :, :, d)); Gi = imag.(view(G, :, :, d))
        S[1:N, 1:K, d]               =  Gr
        S[1:N, (K + 1):2K, d]        = -Gi
        S[(N + 1):2N, 1:K, d]        =  Gi
        S[(N + 1):2N, (K + 1):2K, d] =  Gr
    end
    return S
end

# The factor V of the extra coefficient precision Γ = VᵀV = Σ_kl Σδx^kl M_kᵀ Σ⁻¹ M_l of
# eq. (Gamma). Row r of M depends only on the boundary point of its own sensor, so only the
# Dim × Dim diagonal blocks Σδx^(i) of Σδx contribute and, with the K × Dim gradient block
# G_r of data row r and its sensor i(r), Γ = Σ_r w_r G_r Σδx^(i(r)) G_rᵀ. Factoring each
# block as L_i L_iᵀ gives Γ = VᵀV with the Dim rows √w_r L_i(r)ᵀ G_rᵀ of V per data row r.
# Γ is only ever needed through this factor, so it is never formed: that keeps the
# coefficient posterior in its data-sized form, and costs O(Nr Dim² K) instead of the
# O(Nr Dim² K²) of accumulating the K × K matrix.
function _gamma_factor(model::WorkingModel, gradM::AbstractArray, boundary::BoundaryFactor)
    N, K, Dim = size(gradM)
    w, row_sensor, Σδx = model.w, model.row_sensor, boundary.cov

    Lts = map(1:maximum(row_sensor)) do i
        ks = ((i - 1) * Dim + 1):(i * Dim)
        transpose(_safe_cholesky(Σδx[ks, ks]).L)
    end
    Vγ = zeros(N * Dim, K)
    for r in 1:N
        Vγ[((r - 1) * Dim + 1):(r * Dim), :] =
            sqrt(w[r]) .* (Lts[row_sensor[r]] * transpose(Matrix(gradM[r, :, :])))
    end
    return Vγ
end

# The pieces of the model that depend on where the boundary and the sources are: the
# unwhitened stacked system matrix M (which E-step II and the boundary error need), its
# gradient, and the whitened augmented design Φ = [W M; V].
function _assemble(model::WorkingModel, bd, src_pos, boundary)
    M = _stacked_system_matrix(model, bd, src_pos)
    boundary === nothing && return M, nothing, model.W * M
    gradM = _stacked_gradient(model, bd, src_pos)
    return M, gradM, vcat(model.W * M, _gamma_factor(model, gradM, boundary))
end

# The augmented design alone, used for a single source when its position moves.
_design(model::WorkingModel, bd, src_pos, boundary) = last(_assemble(model, bd, src_pos, boundary))

# ------------------------------------------------------------------------------------
# Bookkeeping of the coefficient columns
# ------------------------------------------------------------------------------------

# Column indices, in the stacked coefficient vector, of the sources `keep` out of `n_src`.
function _kept_columns(model::WorkingModel, keep::AbstractVector{Int}, n_src::Int)
    FD = model.field_dim
    re = reduce(vcat, [collect(((j - 1) * FD + 1):(j * FD)) for j in keep]; init = Int[])
    return model.iscomplex ? vcat(re, re .+ n_src * FD) : re
end

_source_columns(model::WorkingModel, j::Int, n_src::Int) = _kept_columns(model, [j], n_src)

# The sources with at least one coefficient still switched on, and the columns they own.
function _active_sources(model::WorkingModel, α::AbstractVector, n_src::Int, threshold::Real)
    keep = [j for j in 1:n_src if any(<(threshold), α[_source_columns(model, j, n_src)])]
    return keep, _kept_columns(model, keep, n_src)
end

function _initial_precisions(pv::Vector{Float64}, model::WorkingModel, n_src::Int)
    Kh = n_src * model.field_dim
    K = model.iscomplex ? 2Kh : Kh
    length(pv) == 1 && return fill(1 / pv[1], K)
    length(pv) == K && return 1 ./ pv
    (model.iscomplex && length(pv) == Kh) && return vcat(1 ./ pv, 1 ./ pv)
    throw(ArgumentError("prior_variance must be a scalar, or have one entry per coefficient ($Kh) or per real degree of freedom ($K)"))
end

# ------------------------------------------------------------------------------------
# The variational updates
# ------------------------------------------------------------------------------------

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
    return CoefficientPosterior(Σpost * (transpose(Φ) * ĝ), Σpost, logdetΣ)
end

# The noise-weighted expected misfit R of eq. (R) on the whitened augmented model:
# R = ‖ĝ - Φμ‖² + tr(Φ Σ_post Φᵀ). Through the rows [W M̄; V] of Φ this expands into the
# misfit of the mean solution, the extra misfit from the coefficient uncertainty, and the
# boundary-uncertainty terms μᵀΓμ + tr(Γ Σ_post).
function _expected_misfit(Φ::AbstractMatrix, post::CoefficientPosterior, ĝ::AbstractVector)
    r = ĝ - Φ * post.mean
    return dot(r, r) + sum((Φ * post.cov) .* Φ)
end

# The evidence lower bound F of eq. (elbo_model), with the Gaussian KL of eq. (gaussian_kl).
# The boundary factor contributes its own KL, which vanishes identically while q(δx) is
# frozen at the prior — that is, everywhere except after E-step II has moved it.
function _elbo(model::WorkingModel, R::Real, α::AbstractVector,
        post::CoefficientPosterior, boundary::Union{Nothing, BoundaryFactor})

    F = -0.5 * (length(model.g) * log(2π) + model.noise_logdet) - 0.5 * R
    # `init` for the empty model of the first Phase-1 action, which has no coefficients
    F -= 0.5 * (sum(α .* (abs2.(post.mean) .+ diag(post.cov)); init = 0.0) - length(α)
                - sum(log, α; init = 0.0) - post.logdet_cov)
    if boundary !== nothing
        sx2 = boundary.prior_var
        F -= 0.5 * (sum(diag(boundary.cov) ./ sx2)
                    + sum(abs2.(boundary.mean .- boundary.prior_mean) ./ sx2)
                    - length(sx2) + sum(log, sx2) - boundary.logdet_cov)
    end
    return F
end

# E-step II, eqs. (vb_boundary_precision)-(vb_boundary_update): update q(δx) in place. With a
# diagonal prior Σ_x and independent sensors both the precision E[DᵀΣ⁻¹D] and Σ_δx are block
# diagonal, one Dim × Dim block per sensor.
function _update_boundary!(boundary::BoundaryFactor, model::WorkingModel,
        M::AbstractMatrix, gradM::AbstractArray, post::CoefficientPosterior)

    w, g, row_sensor = model.w, model.g, model.row_sensor
    μ, Σpost = post.mean, post.cov
    sx2, m0 = boundary.prior_var, boundary.prior_mean
    Dim = size(gradM, 3)

    n_sensors = maximum(row_sensor)
    rows_of = [Int[] for _ in 1:n_sensors]
    for r in eachindex(row_sensor)
        push!(rows_of[row_sensor[r]], r)
    end

    Mμ = M * μ
    boundary.logdet_cov = 0.0
    for i in 1:n_sensors
        A = zeros(Dim, Dim)
        rhs = zeros(Dim)
        for r in rows_of[i]
            ΣMr = Σpost * M[r, :]
            vs = [Vector(gradM[r, :, d]) for d in 1:Dim]
            Σvs = [Σpost * vs[d] for d in 1:Dim]
            for d in 1:Dim
                vμd = dot(vs[d], μ)
                rhs[d] += w[r] * (vμd * (g[r] - Mμ[r]) - dot(vs[d], ΣMr))
                for d2 in d:Dim
                    val = w[r] * (vμd * dot(vs[d2], μ) + dot(vs[d], Σvs[d2]))
                    A[d, d2] += val
                    d2 > d && (A[d2, d] += val)
                end
            end
        end
        ks = ((i - 1) * Dim + 1):(i * Dim)
        Sblock = inv(Symmetric(A + Diagonal(1 ./ sx2[ks])))
        boundary.cov[ks, ks] = Sblock
        boundary.mean[ks] = Sblock * (rhs + m0[ks] ./ sx2[ks])
        boundary.logdet_cov += logdet(Symmetric(Matrix(Sblock)))
    end
    return boundary
end

# Re-center the boundary at the current estimate: μ_x ← μ_x + μ_δx (Algorithm 1, line 6).
# The outward normals are kept fixed, consistent with small perturbations.
function _recenter_boundary(bd::BoundaryData{F, Dim}, μδx::AbstractVector) where {F, Dim}
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
function _mackay_precisions(α::AbstractVector, post::CoefficientPosterior, α_em::AbstractVector)
    d = diag(post.cov)
    return map(eachindex(α)) do i
        γ = 1 - α[i] * d[i]
        μ2 = abs2(post.mean[i])
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
# Phase 2: moving the source positions by gradient descent
# ------------------------------------------------------------------------------------

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
        cols::AbstractVector{Int}, post::CoefficientPosterior,
        ĝ::AbstractVector, pts, clearance::Real) where Dim

    μ, Σpost = post.mean, post.cov
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

# `iters` rounds of: one descent step on the position of the source with columns `cols` of Φ,
# then a re-solve of the coefficient posterior q(a) — alternating the two lets the source
# travel far while its amplitude follows, where descent at fixed q(a) stalls after a short
# move. Both parts increase the bound F, so the alternation does too. The moved columns are
# written into Φ; returns the position and the refreshed posterior.
function _optimise_source_position!(Φ::AbstractMatrix, cols::AbstractVector{Int}, basis,
        χ0::SVector, α::AbstractVector, ĝ::AbstractVector, iters::Int, pts, clearance::Real)

    χ = χ0
    post = _coefficient_posterior(Φ, α, ĝ)
    for _ in 1:iters
        χnew, moved = _descend_step(χ, basis, Φ, cols, post, ĝ, pts, clearance)
        moved || break
        χ = χnew
        Φ[:, cols] = basis(χ)
        post = _coefficient_posterior(Φ, α, ĝ)
    end
    return χ, post, χ != χ0
end

# ------------------------------------------------------------------------------------
# Phase 1: add one source at a time (fast marginal likelihood, Tipping & Faul)
# ------------------------------------------------------------------------------------

# The log evidence as a function of a single precision α, eq. (sq_factors), up to terms
# independent of α; the limit ℓ(∞) = 0 is the column switched off.
_log_evidence_1(α::Real, s::Real, q::Real) = isinf(α) ? 0.0 : (log(α / (α + s)) + q^2 / (α + s)) / 2

# Relative drift of the S, Q recurrence of `_select_columns`, measured against ‖φᵢ‖², at
# which the factors are re-evaluated exactly.
const SQ_REFRESH_TOL = 1e-9

# The sparsity and quality factors S = diag(ΦᵀC⁻¹Φ) and Q = ΦᵀC⁻¹ĝ of every candidate column,
# evaluated from scratch with C = I + Φ_A A⁻¹ Φ_Aᵀ the evidence covariance of the active model.
# Through Woodbury, C⁻¹ = I - Φ_A Σ Φ_Aᵀ with Σ = (Φ_AᵀΦ_A + A)⁻¹ = (L Lᵀ)⁻¹, both follow from
# Y = L⁻¹B alone, given the Gram rows B = Φ_Aᵀ Φ of the active columns `bcols`:
#     S = diag(ΦᵀΦ) - diag(YᵀY),   Q = Φᵀĝ - Yᵀ(L⁻¹ Φ_Aᵀĝ).
# This is the O(k²K) anchor of the recurrence in `_select_columns`, not its per-action cost.
function _sq_factors(Φnorm2::AbstractVector, Φtĝ::AbstractVector, B::AbstractMatrix,
        bcols::AbstractVector{Int}, α::AbstractVector)

    L = _safe_cholesky(B[:, bcols] + Diagonal(α[bcols])).L
    Y = L \ B
    S = Φnorm2 .- vec(sum(abs2, Y; dims = 1))
    Q = Φtĝ .- transpose(Y) * (L \ Φtĝ[bcols])
    return S, Q
end

# Sparse Bayesian learning at fixed source positions, run constructively: starting from no
# columns, apply the single best action — add a column at its optimal precision
# eq. (alpha_star), re-estimate the precision of an active column, or delete one — until no
# action increases the evidence. Each action is the exact coordinate maximizer of the
# evidence, so the bound increases monotonically.
#
# What an action costs is decided by the sparsity and quality factors S, Q of ALL K candidate
# columns. Read literally they ask for C⁻¹ = (I + Φ_A A⁻¹ Φ_Aᵀ)⁻¹ applied to the whole basis,
# O(P²K) per action, which is what makes a dense candidate set expensive. Two things avoid it:
#
#   * Woodbury the other way round, C⁻¹ = I - Φ_A Σ Φ_Aᵀ with the k × k posterior covariance
#     Σ = (Φ_AᵀΦ_A + A)⁻¹ of the k active columns, whose precision the Gram rows B = Φ_Aᵀ Φ
#     already hold. Selection keeps k ≪ P, the opposite of the regime the data-sized form of
#     `_coefficient_posterior` is written for. This alone brings S, Q down to the O(k²K) of
#     `_sq_factors`;
#   * S and Q are then not rebuilt at all but carried by the rank-1 recurrence of Tipping &
#     Faul — the step the "fast" in fast marginal likelihood refers to. An action changes C
#     by the rank-1 term (1/α_new - 1/α_old) φᵢφᵢᵀ, so by Sherman-Morrison every factor
#     shifts by a multiple of e = ΦᵀC⁻¹φᵢ = φᵢᵀΦ - (Φ_Aᵀφᵢ)ᵀ Σ B, which costs O(kK):
#         S ← S - e²/d,   Q ← Q - e Qᵢ/d,   d = 1/(1/α_new - 1/α_old) + Sᵢ,
#     one expression covering all three actions (an addition has 1/α_old = 0, a deletion
#     1/α_new = 0). No solve against the whole basis survives;
#   * B is carried across actions too: an action changes the active set by a single column,
#     so B gains or loses one row, O(PK) instead of O(PkK).
#
# A recurrence drifts, so it is anchored: eᵢ is φᵢᵀC⁻¹φᵢ evaluated afresh from B and the
# active precision, which the carried Sᵢ must reproduce, and their disagreement measures the
# drift accumulated since the last exact evaluation. Past `SQ_REFRESH_TOL` the action ends
# with an exact `_sq_factors`, which in practice happens a handful of times in a whole run.
#
# The bound needs no re-solve either. With q(a) the exact posterior F is the log evidence, and
# each action is its exact coordinate maximizer, so the evidence gain Δ of the chosen action
# IS the increase of F: the history is accumulated from the bound of the empty model.
#
# Returns the precisions of all K candidate columns — Inf for those never switched on — and
# the evidence-bound history.
function _select_columns(Φ::AbstractMatrix, model::WorkingModel, solver::VariationalBayesianSolver)

    ĝ = model.ĝ
    K = size(Φ, 2)
    α = fill(Inf, K)
    elbo = Float64[]

    # quantities of the whole candidate basis: the positions are fixed throughout Phase 1,
    # so Φ never changes and neither do these
    Φnorm2 = vec(sum(abs2, Φ; dims = 1))
    Φtĝ = transpose(Φ) * ĝ

    # The active columns `bcols`, in activation order, and their rows B = Φ_Aᵀ Φ of the Gram
    # matrix: an addition appends a row and a deletion swaps in the last one, so each row is
    # formed once per activation. `B` is a buffer grown geometrically, of which only the
    # first `length(bcols)` rows are in use.
    bcols = Int[]
    B = Matrix{Float64}(undef, 0, K)

    # the sparsity and quality factors, carried across actions by the recurrence below. With
    # no active column C = I, so they start at their anchor values diag(ΦᵀΦ) and Φᵀĝ.
    S, Q = _sq_factors(Φnorm2, Φtĝ, view(B, eachindex(bcols), :), bcols, α)

    # the bound of the empty model, from which the history is accumulated: no coefficients,
    # so the KL term vanishes and the expected misfit is ‖ĝ‖²
    F_empty = _elbo(model, dot(ĝ, ĝ), Float64[], _empty_posterior(), nothing)
    F_prev = F_empty

    # one action per iteration, so this phase needs its own budget: a negative
    # `max_select_iters` asks for the automatic four-per-candidate-coefficient safety net
    budget = solver.max_select_iters < 0 ? 4K : solver.max_select_iters

    for _ in 1:budget
        # the single best action: for each column the leave-one-out factors s, q of
        # eq. (sq_factors), its optimal precision eq. (alpha_star), and the evidence gain.
        # S is a Schur complement and so non-negative; rounding can push it just below zero
        # for a column that the active set already explains to machine precision, and the
        # `s > 0` guard drops such a column as uninformative.
        best_i, best_Δ, best_α = 0, solver.elbo_tol * (1 + abs(F_prev)), Inf
        for i in 1:K
            s, q = if isinf(α[i])
                S[i], Q[i]
            else            # S, Q include column i itself; remove its own contribution
                den = max(α[i] - S[i], 1e-12 * α[i])
                (α[i] * S[i] / den, α[i] * Q[i] / den)
            end
            α_i = (s > 0 && q^2 > s) ? s^2 / (q^2 - s) : Inf
            Δ = _log_evidence_1(α_i, s, q) - _log_evidence_1(α[i], s, q)
            if Δ > best_Δ
                best_i, best_Δ, best_α = i, Δ, α_i
            end
        end
        best_i == 0 && break        # no action improves the evidence: converged

        added = isinf(α[best_i])
        deleted = !added && isinf(best_α)
        r = added ? 0 : findfirst(==(best_i), bcols)     # the row of column best_i in B

        # e = ΦᵀC⁻¹φᵢ = φᵢᵀΦ - (Φ_Aᵀφᵢ)ᵀ Σ B, the direction along which this action shifts
        # every S and Q, all of it from the CURRENT active model. The Gram row φᵢᵀΦ is in B
        # already unless the column is being added, in which case it is the row B is about to
        # gain anyway; Σ acts through the Cholesky of the active precision Φ_AᵀΦ_A + A.
        Bv = view(B, eachindex(bcols), :)
        gram_i = added ? transpose(Φ) * view(Φ, :, best_i) : B[r, :]
        L = _safe_cholesky(Bv[:, bcols] + Diagonal(α[bcols])).L
        e = gram_i .- transpose(Bv) * (transpose(L) \ (L \ Bv[:, best_i]))

        # eᵢ is φᵢᵀC⁻¹φᵢ freshly evaluated, which the carried Sᵢ must reproduce: how far apart
        # they are is the drift the recurrence has accumulated
        drifted = abs(e[best_i] - S[best_i]) > SQ_REFRESH_TOL * Φnorm2[best_i]

        # Sherman-Morrison, one expression for all three actions: an addition comes from
        # α_old = ∞ and a deletion from α_new = ∞
        κinv = added ? best_α :
            deleted ? -α[best_i] : α[best_i] * best_α / (α[best_i] - best_α)
        d = κinv + S[best_i]
        Qi = Q[best_i]
        @. S -= e^2 / d
        @. Q -= e * Qi / d

        α[best_i] = best_α

        # keep B in step with the active set: one row appended, or one swapped out
        if added
            nb = length(bcols) + 1
            if nb > size(B, 1)
                Bnew = Matrix{Float64}(undef, max(8, 2 * size(B, 1)), K)
                copyto!(view(Bnew, 1:(nb - 1), :), view(B, 1:(nb - 1), :))
                B = Bnew
            end
            copyto!(view(B, nb, :), gram_i)
            push!(bcols, best_i)
        elseif deleted
            copyto!(view(B, r, :), view(B, length(bcols), :))
            bcols[r] = bcols[end]
            pop!(bcols)
        end

        # the chosen action is the exact coordinate maximizer of the evidence, so its gain is
        # the increase of F, and no re-solve is needed to monitor the bound
        F_prev += best_Δ
        push!(elbo, F_prev)

        drifted && ((S, Q) = _sq_factors(Φnorm2, Φtĝ, view(B, eachindex(bcols), :), bcols, α))
    end

    if !any(isfinite, α)
        @warn "no candidate source explains the data beyond the noise; keeping the best one"
        α[argmax(abs2.(Φtĝ) ./ Φnorm2)] = ALPHA_CAP
    end

    # The history above was accumulated from per-action evidence gains, never re-derived. Close
    # the loop the way the section "Monitoring the bound" prescribes: evaluate F exactly at the
    # selected model and compare. Agreement certifies the whole S, Q recurrence; disagreement is
    # the only symptom a drifted recurrence has. The comparison is against the evidence the
    # phase claims to have GAINED, not against F itself: the gain is the sum of a few thousand
    # exact per-action increments, so a small relative slip in each is expected and harmless,
    # whereas a discrepancy of the order of the gain means the recurrence — and with it the
    # choice of sources — has come apart.
    if !isempty(bcols) && !isempty(elbo)
        post = _coefficient_posterior(Φ[:, bcols], α[bcols], ĝ)
        F_exact = _elbo(model, _expected_misfit(Φ[:, bcols], post, ĝ), α[bcols], post, nothing)
        gain = elbo[end] - F_empty
        if abs(F_exact - elbo[end]) > 1e-3 * max(abs(gain), 1 + abs(F_exact))
            @warn(
                "the evidence recurrence of the source selection has drifted from the bound it " *
                "claims to maximize: the selected set is not trustworthy",
                accumulated = elbo[end], exact = F_exact, gained = gain, selected = length(bcols)
            )
        end
        elbo[end] = F_exact
    end

    return α, elbo
end

# ------------------------------------------------------------------------------------
# solve: Phase 1 (select the sources), then Phase 2 (variational EM, Algorithm 1)
# ------------------------------------------------------------------------------------

function solve(sim::Simulation{VariationalBayesianSolver, Dim}) where Dim

    solver = sim.solver
    bd = sim.boundary_data
    src_pos = Vector{SVector{Dim, Float64}}(sim.source_positions)

    # An uncertain boundary enters at two levels, controlled separately. The boundary factor
    # q(δx) exists as soon as Σ_x is to be used at all: frozen at its prior it contributes
    # only the ridge Γ of eq. (Gamma_x) — the closed-form evidence bound of the section
    # "A closed-form evidence bound: boundary variance as pseudo-data" — and E-step II then
    # additionally learns δx. Phase 1 always runs in the first regime, as the theory prescribes.
    sx2 = _boundary_variances(sim, solver)
    boundary = sx2 === nothing ? nothing : BoundaryFactor(sx2)
    learn_geometry = solver.update_geometry_flag && boundary !== nothing

    model = _working_model(sim, boundary !== nothing)
    if boundary !== nothing && !(model.w isa AbstractVector)
        throw(ArgumentError("an uncertain boundary requires independent (diagonal) sensor noise; the boundary `fields` covariance must be diagonal"))
    end

    # INFERRING the boundary rests on the linearization eq. (linearization_M), and M varies on
    # the scale of the distance from a source to the boundary: a source closer to the boundary
    # than a few standard deviations of Σ_x has no meaningful first-order expansion, so any δx
    # estimated through it is meaningless rather than merely uncertain. The ridge alone is not
    # subject to this — Γ is a well-defined quadratic penalty on the coefficients whatever Σ_x
    # is taken to mean — so this warns only when E-step II is going to move the boundary.
    prior_scale = learn_geometry ? sqrt(maximum(sx2)) : 0.0
    if learn_geometry
        pts0 = mean_points(bd)
        too_close = count(χ -> minimum(norm(χ - p) for p in pts0) < 3 * prior_scale, src_pos)
        too_close > 0 && @warn(
            "some sources lie closer to the boundary than the boundary is uncertain, where " *
            "the linearization of M in δx — and with it the ridge Γ — has no meaning; move " *
            "them further in, or reduce Σ_x",
            boundary_prior_std = prior_scale,
            sources_within_3_std = too_close,
            of = length(src_pos)
        )
    end

    M, gradM, Φ = _assemble(model, bd, src_pos, boundary)
    elbo_history = Float64[]

    # --- Phase 1: add one source at a time by the exact evidence criterion ---
    if solver.select_sources_flag && solver.learn_prior_flag
        α_all, elbo1 = _select_columns(Φ, model, solver)
        keep, cols = _active_sources(model, α_all, length(src_pos), Inf)
        α = min.(α_all[cols], ALPHA_CAP)
        src_pos = src_pos[keep]
        M, Φ = M[:, cols], Φ[:, cols]
        gradM === nothing || (gradM = gradM[:, cols, :])
        append!(elbo_history, elbo1)
    else
        α = _initial_precisions(solver.prior_variance, model, length(src_pos))
    end

    # --- Phase 2: variational EM on the selected sources (Algorithm 1) ---
    x_nominal = Vector{Float64}(vcat(mean_points(bd)...))
    baseline_resets = Int[]
    recenter_iterations = Int[]
    # continue monitoring from the last Phase-1 value: the bound is comparable across the
    # transition, so the MacKay fallback also guards the first Phase-2 iteration
    F_prev = isempty(elbo_history) ? -Inf : elbo_history[end]
    post = _empty_posterior()

    for _ in 1:solver.max_iters
        reset_baseline = false

        # --- E-step I: update q(a), eq. (vb_coefficient_update) ---
        post = _coefficient_posterior(Φ, α, model.ĝ)

        # --- E-step II: update q(δx) and re-center, eq. (vb_boundary_update) ---
        if learn_geometry
            _update_boundary!(boundary, model, M, gradM, post)
            step = norm(boundary.mean, Inf)
            if step > 1e-6 * prior_scale
                bd = _recenter_boundary(bd, boundary.mean)
                fill!(boundary.mean, 0.0)
                boundary.prior_mean = x_nominal - Vector{Float64}(vcat(mean_points(bd)...))
                push!(recenter_iterations, length(elbo_history) + 1)
                # only a significant move of the linearization point counts as a model
                # change, across which values of the bound are not comparable
                reset_baseline = step > 1e-2 * prior_scale
            end
            M, gradM, Φ = _assemble(model, bd, src_pos, boundary)
        end

        # --- M-step: prior precisions, EM update eq. (alpha_update) ---
        α_em = 1 ./ (abs2.(post.mean) .+ diag(post.cov))
        if solver.learn_prior_flag
            α_new = solver.mackay_acceleration_flag ? _mackay_precisions(α, post, α_em) : α_em
            α = min.(α_new, ALPHA_CAP)
        end

        # --- ARD pruning: remove sources whose every precision has diverged ---
        if solver.learn_prior_flag && solver.ard_prune_flag
            n_active = length(src_pos)
            keep, cols = _active_sources(model, α, n_active, solver.ard_threshold)
            if isempty(keep)
                @warn "automatic relevance determination switched off every source; keeping the most relevant one"
                keep = [argmin([minimum(α[_source_columns(model, j, n_active)]) for j in 1:n_active])]
                cols = _kept_columns(model, keep, n_active)
            end
            if length(keep) < n_active
                α, α_em = α[cols], α_em[cols]
                post = _restrict(post, cols)
                M, Φ = M[:, cols], Φ[:, cols]
                gradM === nothing || (gradM = gradM[:, cols, :])
                src_pos = src_pos[keep]
                reset_baseline = true
            end
        end

        # --- M-step: source positions χ, a few gradient-descent steps per source ---
        if solver.optimise_source_positions_flag && solver.source_position_iters > 0
            pts = mean_points(bd)
            clearance = _boundary_spacing(pts) / 2
            moved_any = false
            for j in eachindex(src_pos)
                cols = _source_columns(model, j, length(src_pos))
                basis = χ -> _design(model, bd, [χ], boundary)
                χ, post, moved = _optimise_source_position!(Φ, cols, basis, src_pos[j], α,
                    model.ĝ, solver.source_position_iters, pts, clearance)
                if moved
                    src_pos[j] = χ
                    moved_any = true
                end
            end
            moved_any && ((M, gradM, Φ) = _assemble(model, bd, src_pos, boundary))
        end

        # --- monitor the bound, eq. (elbo_model) ---
        R = _expected_misfit(Φ, post, model.ĝ)
        F = _elbo(model, R, α, post, boundary)

        # MacKay acceleration carries no monotonicity guarantee: fall back to the EM update
        # whenever the bound fails to increase (Section on ARD of the theory document).
        if solver.learn_prior_flag && solver.mackay_acceleration_flag && !reset_baseline && F < F_prev
            α = min.(α_em, ALPHA_CAP)
            F = _elbo(model, R, α, post, boundary)
        end

        push!(elbo_history, F)
        reset_baseline && push!(baseline_resets, length(elbo_history))

        converged = !reset_baseline && abs(F - F_prev) <= solver.elbo_tol * (1 + abs(F))
        F_prev = F
        converged && break
    end

    # --- final inference at the learned hyperparameters ---
    post = _coefficient_posterior(Φ, α, model.ĝ)
    misfit_ratio = _expected_misfit(Φ, post, model.ĝ) / length(model.g)
    relative_boundary_error = norm(M * post.mean - model.g) / norm(model.g)

    # complex coefficients keep their posterior covariance over the stacked real degrees
    # of freedom [Re a; Im a], the convention of FundamentalSolution
    coefficients = if model.iscomplex
        Khh = length(src_pos) * model.field_dim
        post.mean[1:Khh] .+ im .* post.mean[(Khh + 1):end]
    else
        copy(post.mean)
    end

    fsol = FundamentalSolution(model.medium;
        positions = collect(src_pos),
        coefficients = coefficients,
        coefficients_covariance = post.cov,
        particular_solution = sim.particular_solution,
        relative_boundary_error = relative_boundary_error
    )

    # all the boundary information is returned as a BoundaryShape: when the geometry was
    # updated, the posterior of the boundary points is an MvNormal with the re-centered
    # boundary as mean and Σ_δx as covariance; otherwise the input shape is passed through
    boundary_shape = if learn_geometry
        shape = bd.boundary_shape
        posterior_points = MvNormal(
            Vector{Float64}(vcat(mean_points(bd)...)),
            Symmetric(boundary.cov)
        )
        BoundaryShape(posterior_points, shape.normals, shape.interior_points)
    else
        sim.boundary_data.boundary_shape
    end

    return VariationalSolution(
        fsol, boundary_shape, α, elbo_history, baseline_resets, recenter_iterations, misfit_ratio
    )
end
