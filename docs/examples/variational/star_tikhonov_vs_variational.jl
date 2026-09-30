# # Where do you put the sources? Tikhonov MFS against the variational solver
#
# One obstacle — a smooth five-lobed star, r(θ) = 0.75(1 + 0.2 sin 5θ), ANALYTIC everywhere —
# with a Dirichlet (sound-soft) condition, lit symmetrically by five incident point sources, one
# on each of its five symmetry axes. The question the example asks is the only question
# classical MFS really has: WHERE DO THE SOURCES GO?
#
# Classical MFS answers it by hand: put one source behind every boundary point, a fixed
# offset in from the boundary. The offset is a free parameter, and the whole method lives or
# dies by it. The variational solver answers it from the data: lay candidates down
# EVERYWHERE inside the obstacle, let sparse Bayesian learning SELECT the useful ones
# (Phase 1), then let variational EM MOVE the survivors (Phase 2).
#
# Everything here is judged by ONE number, measured the same way for both methods: the
# residual of the boundary condition at 2000 points spread along the WHOLE boundary — ten
# times denser than the 200 sensors either method is allowed to fit, and interleaved with
# them, so no test point is a point that was fitted. The Dirichlet condition says the total
# pressure vanishes there, so |total field| at a test point IS the error, and we report it
# relative to the strongest incident signal on the boundary.
#
# BOTH SOLVERS ARE RUN ON THE SAME DATA. On an analytic boundary classical MFS with EXACT
# boundary values is superb — rms 3e-5 at the right offset — which is far below any noise floor,
# so comparing that against a solver handed 0.3% noisy data would measure the data and not the
# method. Section 1 therefore sweeps the offset BOTH ways, and section 2 compares like with
# like. (On the corner version of this obstacle this distinction was invisible: the corner error
# was ~100× the noise floor and dominated every number in the example.)
#
# The five things it shows:
#
#   1. TIKHONOV CANNOT SEE ITS OWN ERROR. With exact data the offset spans four decades of
#      true error, and the residual at the fitted sensors points the WRONG WAY across the whole
#      useful range: at a quarter step the residual is 2.6e-9 and the error 185× worse than at
#      one and a half steps, where the residual is 38000× LARGER. With noisy data the offset
#      barely matters until the system falls apart — the offset sensitivity everyone worries
#      about is a clean-data phenomenon.
#   2. VARIATIONAL, sources everywhere → ARD selects → Phase 2 moves them. Nobody chooses an
#      offset. Like for like on the same noisy data it TIES the oracle-tuned Tikhonov (rms
#      0.0035 against 0.0032) with fewer than half the sources, and its own misfit diagnostic
#      R/N_g ≈ 0.6 says, correctly, that it has fitted the data to the noise and no further.
#   3. UNCERTAIN SENSOR POSITIONS, AND WHICH DIRECTION THEY ARE UNCERTAIN IN. A boundary
#      covariance Σx enters as the ridge Γ = Σᵣ wᵣ Gᵣ Σδx Gᵣᵀ, and its DIRECTION decides what
#      that ridge penalises. Taken ALONG THE BOUNDARY TANGENT, Σδx = σ²ttᵀ makes
#          aᵀΓa = σ² Σᵣ wᵣ |∂(Ma)ᵣ/∂s|²,
#      a penalty on the arc-length derivative of the fitted boundary trace — a smoothness
#      prior, and precisely the thing that suppresses oscillation between collocation points.
#      Taken along the NORMAL it penalises ∂/∂n instead, which does neither. On a SMOOTH
#      boundary with the sensors correctly placed the fit is already smooth and the prior has
#      almost nothing to do; it earns its keep when the sensors are genuinely MISPLACED, and
#      there the tangential direction is the better of the two on error and on roughness alike.
#   4. THE SOURCES MOVE WITH THE DATA ERROR. Sweep the sensor noise σ and watch the selected
#      set change: from 0.1% to 10% noise it shrinks by a factor of three, and the DEEP sources
#      go first — the 90th-percentile depth falls from 4.6 steps to 1.5. Deep sources overlap
#      each other on the boundary, so their individual amplitudes are the first thing noise
#      makes unidentifiable. Classical MFS has no such mechanism: same layout at any noise.
#   5. THE LAYOUT IS NOT A CLASSICAL ONE. The selected sources do not sit on an offset curve.
#      Their depths spread over two decades where classical MFS uses one, Phase 2 pushes a few
#      of them deep into the interior, and there are half as many of them as boundary points.
#      Curvature still attracts them, but far more weakly than on the corner version: the ten
#      curvature extrema cover 30% of the boundary and hold 45% of the sources, a 1.5×
#      enrichment where the corner version piled them into the five tips outright.
#
# Needs Plots.jl (not a dependency of this package): run from an environment where both
# `using MethodOfFundamentalSolutions` and `using Plots` work.

using MethodOfFundamentalSolutions
using Distributions, LinearAlgebra, StaticArrays, Statistics, Random
using Plots
gr()

const FIGDIR = joinpath(@__DIR__, "figures")
mkpath(FIGDIR)
Random.seed!(2026)

# ---------------------------------------------------------------------------------------------
# The obstacle: a smooth five-lobed star, r(θ) = 0.75(1 + 0.2 sin 5θ), ANALYTIC everywhere.
#
# The smoothness is deliberate. An earlier version used a triangle wave, giving ten genuine
# corners, and they dominated everything: at a convex tip the EXTERIOR wedge angle is 259°, so
# the local solution goes like r^(π/ω) = r^0.69 and its gradient is unbounded. The analytic
# continuation of the scattered field into the obstacle then has a BRANCH POINT sitting exactly
# on the tip — not a pole, so no finite sum of fundamental solutions can represent it, and
# every method below committed its largest error at the same five points. Measured there: the
# error in the tip band ran 18× the median for Tikhonov at half a step and 460× at one step.
# That is a statement about corner singularities, not about where sources belong, and it
# drowned out the comparison this example is for. With an analytic boundary the continuation
# reaches a finite depth, the sources have something they can represent, and what is left is
# the question of WHERE.
# ---------------------------------------------------------------------------------------------
star_r(θ) = 0.75 * (1 + 0.2 * sin(5θ))
boundary_point(θ) = star_r(θ) * SVector(cos(θ), sin(θ))

function boundary_normal(θ; h = 1e-5)
    t = boundary_point(θ + h) - boundary_point(θ - h)
    n = SVector(t[2], -t[1]) / norm(SVector(t[2], -t[1]))
    return dot(n, boundary_point(θ)) < 0 ? -n : n
end

inside(x) = norm(x) < star_r(atan(x[2], x[1]))

# Points equally spaced in ARC LENGTH. Uniform θ would be a trap: on the steep lobe flanks the
# curve races through a lot of arc per unit θ, so uniform-θ points leave exactly the parts with
# the most curvature under-sampled. The half-step offset keeps the sensors interleaved with the
# test points below.
function arc_length_thetas(N; dense = 20000)
    θd = LinRange(0, 2π, dense + 1)
    seg = [norm(boundary_point(θd[i + 1]) - boundary_point(θd[i])) for i in 1:dense]
    cum = cumsum(seg)
    return [θd[searchsortedfirst(cum, s)] for s in ((1:N) .- 0.5) .* (cum[end] / N)], cum[end]
end

N_fit = 200
θ_fit, perimeter = arc_length_thetas(N_fit)
pts = [boundary_point(θ) for θ in θ_fit]
normals = [boundary_normal(θ) for θ in θ_fit]

# the inter-boundary step: the classical MFS offset, and the yardstick for every distance below
step = perimeter / N_fit

# Is the boundary actually resolved? For a smooth curve the two scales that matter are the
# wavelength and the radius of curvature; the sampling has to beat both.
function curvature(θ; h = 1e-4)
    p0, pm, pp = boundary_point(θ), boundary_point(θ - h), boundary_point(θ + h)
    d1 = (pp - pm) / (2h); d2 = (pp - 2p0 + pm) / h^2
    return (d1[1] * d2[2] - d1[2] * d2[1]) / norm(d1)^3
end
κ = curvature.(LinRange(0, 2π, 4001)[1:4000])
min_radius = 1 / maximum(abs, κ)

# The test set: ten times denser, and interleaved with the sensors by construction — the k-th
# sensor sits at arc length (k - 1/2)·L/200 and the j-th test point at (j - 1/2)·L/2000, and
# those never coincide. This is where every error below is measured.
N_test = 2000
θ_test, _ = arc_length_thetas(N_test)
test_pts = [boundary_point(θ) for θ in θ_test]
test_arc = ((1:N_test) .- 0.5) .* (perimeter / N_test)

# ---------------------------------------------------------------------------------------------
# Physics: a sound-soft star, wavelength 1, so the obstacle spans about 1.8 wavelengths. The
# tightest radius of curvature is in the concave valleys, about a tenth of a wavelength and
# still several inter-boundary steps, so the boundary is resolved everywhere (printed below).
#
# The illumination is SYMMETRIC: five incident point sources on a ring, one on each of the
# star's five symmetry axes. That matters for both of the things this example measures.
# It makes the results readable — the whole problem now has the star's own five-fold symmetry,
# so every spike is lit identically, the boundary error repeats five times, and any structure
# left in a source layout is structure the method chose rather than a shadow. And it is the
# fairest possible setting for classical MFS: a single off-centre source leaves half the
# boundary in shadow, where a uniform offset curve spends sources on a part of the boundary
# that carries almost no signal, while the variational solver simply declines to.
#
# Exact radial symmetry — one source at the star's centre, or a dense ring — is deliberately
# NOT used: an incident trace that is exactly radially symmetric on the boundary is matched
# exactly by a SINGLE monopole at the centre, since that monopole is regular and radiating
# everywhere outside. The exact answer would then be one source, the total field would vanish
# identically, and there would be nothing left to compare. Five sources give the symmetry
# without the degeneracy: the trace keeps its m = ±5, ±10, … angular content, which no single
# source can reproduce.
# ---------------------------------------------------------------------------------------------
ω = 2π
medium = Acoustic(2; ω = ω, ρ = 1.0, c = 1.0)

incident_radius = 3.0
incident_positions = [incident_radius * SVector(cos(θ), sin(θ))
    for θ in (π / 10) .+ (2π / 5) .* (0:4)]      # π/10: the angles of the five lobe tips
incident = PointSource(incident_positions, ones(ComplexF64, 5))

incident_trace(p) = field(TractionType(), medium, incident, p, SVector(1.0, 0.0))[1]
incident_scale = maximum(abs, incident_trace.(pts))

@info "the boundary is smooth and resolved" perimeter = round(perimeter, digits = 3) step = round(step, digits = 4) min_radius_of_curvature = round(min_radius, digits = 3) points_per_min_radius = round(min_radius / step, digits = 1) points_per_wavelength = round((2π / ω) / step, digits = 1)

# The sensor noise, fixed once here because BOTH solvers have to see it. σ is a fraction of the
# strongest incident signal on the boundary.
σ_rel = 0.003         # 0.3%

# THE error measure of this example: |total pressure| along the whole boundary, relative to
# the strongest incident signal on it. `where_pts = pts` restricts it to the fitted sensors.
boundary_error(sol; where_pts = test_pts) =
    [abs(field(TractionType(), sol, p)[1]) for p in where_pts] ./ incident_scale

# How ROUGH the fitted trace is along the boundary: the rms of its arc-length derivative,
# expressed as the variation across one inter-boundary step relative to the incident amplitude.
# This is the quantity a tangential Σx penalises, so it is what section 3 has to move.
function roughness(sol)
    f = [field(TractionType(), sol, p)[1] for p in test_pts]
    ds = perimeter / N_test
    return sqrt(mean(abs2, (circshift(f, -1) .- f) ./ ds)) * step / incident_scale
end

# a dense sampling of the true curve, to measure how deep inside a source really sits
curve = [boundary_point(θ) for θ in LinRange(0, 2π, 6001)[1:6000]]
depth(x) = minimum(norm(x - q) for q in curve)

# =============================================================================================
# 1. TIKHONOV: one source per boundary point, a fixed offset in from the boundary
# =============================================================================================
# BOTH solvers see the SAME data: the Dirichlet condition as a NOISY measurement, the identical
# realisation from the identical seed, written here as complex scalars because that is the form
# the Tikhonov solver takes and as `MvNormal`s in section 2 because that is what the variational
# solver takes. Exact boundary values are not used anywhere in this example. They would flatter
# Tikhonov out of all proportion — on an analytic boundary it converges to rms 3e-5 with exact
# data, far below any noise floor — and no measured boundary condition is ever exact.
noisy_zeros(σ; seed = 4242) = (Random.seed!(seed); [σ .* randn(2) for _ in pts])
bd_noisy_tik = BoundaryData(TractionType();
    boundary_points = pts,
    fields = [SVector(0.0 .* complex(z[1], z[2])) for z in noisy_zeros(σ_rel * incident_scale)],
    normals = normals,
    interior_points = [SVector(0.0, 0.0)]
)

# sources at p - m·step·n: INWARD, because for an exterior scattering problem every source
# must lie inside the obstacle, where it may be singular
tikhonov_sources(m) = [pts[i] - m * step * normals[i] for i in eachindex(pts)]

function solve_tikhonov(m)
    src = tikhonov_sources(m)
    all(inside, src) || @warn "offset $(m) step: some sources fall outside the star" n = count(!inside, src)
    return solve(Simulation(medium, bd_noisy_tik;
        source_positions = src, particular_solution = incident, solver = TikhonovSolver()))
end

offsets = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
tik_sweep = map(offsets) do m
    f = solve_tikhonov(m)
    (offset = m, fit = maximum(boundary_error(f; where_pts = pts)), all = maximum(boundary_error(f)),
        rms = sqrt(mean(boundary_error(f) .^ 2)))
end

for r in tik_sweep
    @info "Tikhonov" offset = "$(r.offset) steps" sensor_residual = round(r.fit, sigdigits = 3) whole_boundary = round(r.all, sigdigits = 3) rms = round(r.rms, sigdigits = 3)
end

# the best offset, chosen with knowledge no user of the method has: this is what section 2 has
# to beat, and it is a deliberately generous baseline
best_tik = tik_sweep[argmin([r.rms for r in tik_sweep])]
@info "the best Tikhonov offset, picked by an oracle" offset = best_tik.offset rms = round(best_tik.rms, sigdigits = 3) whole_boundary_max = round(best_tik.all, sigdigits = 3) spread_over_the_sweep = round(maximum(r.rms for r in tik_sweep if r.rms < 0.1) / best_tik.rms, digits = 2)

fsol_half = solve_tikhonov(0.25)          # crowded against the boundary: the worst usable offset
fsol_tik  = solve_tikhonov(best_tik.offset)   # the oracle's choice, the baseline for section 2

err_half = boundary_error(fsol_half)
err_tik  = boundary_error(fsol_tik)
err_tik_fit = boundary_error(fsol_tik; where_pts = pts)

# =============================================================================================
# 2. VARIATIONAL: candidates everywhere, ARD selects (Phase 1), variational EM moves (Phase 2)
# =============================================================================================
# A blind regular grid carpeting the inside of the star. It does not know where the lobes
# are; the clearance is well under the lobe half-width so that candidates DO reach into the
# tips — no solver can keep a source that was never offered.
grid_spacing, grid_clearance = 0.02, 0.006
xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
candidates = filter(vec([SVector(x, y)
        for x in minimum(xs):grid_spacing:maximum(xs), y in minimum(ys):grid_spacing:maximum(ys)])) do p
    inside(p) && depth(p) > grid_clearance
end

@info "an overcomplete basis" candidates = length(candidates) coefficients = 2 * length(candidates) data_rows = 2 * N_fit

# The Dirichlet condition as NOISY data: every sensor reports zero total pressure to within a
# known σ. This is how the variational solver is told the measurement noise, and it is the
# only place a noise level is ever specified.
dirichlet_data(σ) = [MvNormal(σ .* randn(2), σ^2 * I(2)) for _ in pts]

function star_data(σ_rel; points = pts, Σx = nothing, seed = 4242)
    Random.seed!(seed)
    fields = dirichlet_data(σ_rel * incident_scale)
    bpts = Σx === nothing ? points : MvNormal(vcat(Vector.(points)...), Σx)
    return BoundaryData(TractionType(); boundary_points = bpts, fields = fields,
        normals = normals, interior_points = [SVector(0.0, 0.0)])
end

bd_noisy = star_data(σ_rel)     # σ_rel was fixed above, so both solvers see the same noise

# --- Phase 1: SELECT. Starting from an EMPTY basis the solver repeatedly adds the single
#     candidate coefficient that most increases the evidence, each at its exact optimal prior
#     precision, and stops when no candidate is worth adding. No offset is chosen anywhere.
solver_select = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 200, elbo_tol = 1e-9)
t_select = @elapsed v_ard = solve(Simulation(medium, bd_noisy;
    source_positions = candidates, particular_solution = incident, solver = solver_select))

pos_ard = copy(v_ard.fsol.positions)
err_ard = boundary_error(v_ard)

@info "Phase 1 (select)" kept = length(pos_ard) of = length(candidates) seconds = round(t_select, digits = 1) misfit_ratio = round(v_ard.misfit_ratio, digits = 3) whole_boundary_max = round(maximum(err_ard), sigdigits = 3) rms = round(sqrt(mean(err_ard .^ 2)), sigdigits = 3)

# --- Phase 2: MOVE. Algorithm 1 on the survivors, entered exactly as the theory requires —
#     at the positions AND the precisions Phase 1 finished with (`select_sources_flag = false`
#     skips a second selection; `prior_variance` hands the learned precisions back). Now the
#     source positions are hyperparameters too, updated by gradient steps that are only
#     accepted when they increase the bound.
solver_move = VariationalBayesianSolver(
    prior_variance = 1 ./ v_ard.prior_precisions,
    select_sources_flag = false,           # Phase 1 is done: do not select again
    ard_prune_flag = true,                 # but Algorithm 1 line 8 still prunes
    optimise_source_positions_flag = true,
    source_position_iters = 3,
    ard_threshold = 1e10, max_iters = 60, elbo_tol = 1e-9
)
t_move = @elapsed v_move = solve(Simulation(medium, bd_noisy;
    source_positions = pos_ard, particular_solution = incident, solver = solver_move))

pos_move = copy(v_move.fsol.positions)
err_move = boundary_error(v_move)
err_move_fit = boundary_error(v_move; where_pts = pts)

# How far did they actually go? Phase 2 both moves and prunes, so the lists are not index-aligned;
# a source can move at most a quarter of its distance to the boundary per accepted step, so the
# nearest Phase-1 position is where each survivor came from.
travel_phase2 = [minimum(norm(b - a) for a in pos_ard) / step for b in pos_move]

@info "Phase 2 (move)" sources = length(pos_move) seconds = round(t_move, digits = 1) misfit_ratio = round(v_move.misfit_ratio, digits = 3) whole_boundary_max = round(maximum(err_move), sigdigits = 3) rms = round(sqrt(mean(err_move .^ 2)), sigdigits = 3) median_travel_in_steps = round(median(travel_phase2), digits = 3) max_travel_in_steps = round(maximum(travel_phase2), digits = 3)

# --- Phase 2 again, with a TANGENTIAL smoothing prior. Nothing about the sensors has changed
#     — they are exactly on the boundary here — but asserting that they are only located to
#     ±σ ALONG the boundary turns Γ into σ² Σᵣ wᵣ |∂(Ma)ᵣ/∂s|², a penalty on the arc-length
#     derivative of the fitted trace. It is a smoothness prior wearing a geometry prior's
#     clothes, aimed at a fit that honours its 200 collocation points and wanders between them.
#     On THIS obstacle it has very little to bite on — the boundary is analytic, the sources are
#     placed by evidence, and the trace is already smooth, so the prior changes the roughness by
#     a few percent and the error not at all. It is kept here because it is the same object
#     section 3 uses, where the sensors are genuinely wrong and it does real work.
σ_smooth = 0.25 * step
smoothing_prior = Matrix(Diagonal(vcat([
    [σ_smooth^2 * n[2]^2 + (0.1σ_smooth)^2, σ_smooth^2 * n[1]^2 + (0.1σ_smooth)^2]
for n in normals]...)))          # the tangent is the normal rotated by π/2: t = (-n₂, n₁)

#     The positions are held FIXED for this one, so that the only difference from `v_ard` is
#     the prior. Combining the smoothing prior WITH the position M-step is not done: the descent
#     then minimises the PENALISED misfit, so the sources chase the smoothness term instead of
#     the data, which costs both accuracy and a great deal of run time. Use one or the other.
solver_smooth = VariationalBayesianSolver(
    prior_variance = 1 ./ v_ard.prior_precisions,
    select_sources_flag = false,
    ard_prune_flag = true,
    ard_threshold = 1e10, max_iters = 100, elbo_tol = 1e-10
)
bd_smooth = star_data(σ_rel; Σx = smoothing_prior)
t_smooth = @elapsed v_smooth = solve(Simulation(medium, bd_smooth;
    source_positions = pos_ard, particular_solution = incident, solver = solver_smooth))

pos_smooth = copy(v_smooth.fsol.positions)
err_smooth = boundary_error(v_smooth)

@info "Phase 2 with a tangential smoothing prior (positions fixed)" sources = length(pos_smooth) seconds = round(t_smooth, digits = 1) whole_boundary_max = round(maximum(err_smooth), sigdigits = 3) rms = round(sqrt(mean(err_smooth .^ 2)), sigdigits = 3) rms_without_prior = round(sqrt(mean(err_ard .^ 2)), sigdigits = 3) roughness_without = round(roughness(v_ard), sigdigits = 3) roughness_with = round(roughness(v_smooth), sigdigits = 3)

# Does the posterior know where it is wrong? The solution is a distribution, so every boundary
# point carries a standard deviation s(x) = √(Var Re f + Var Im f) of the scattered field.
post_std = [sqrt(sum(abs2, field_std(TractionType(), v_move, p))) for p in test_pts] ./ incident_scale
@info "the posterior as an error indicator" correlation_of_logs = round(cor(log.(err_move .+ 1e-12), log.(post_std .+ 1e-12)), digits = 3) median_std = round(median(post_std), sigdigits = 3) median_error = round(median(err_move), sigdigits = 3)

# =============================================================================================
# 3. UNCERTAIN SENSOR POSITIONS: THE TANGENTIAL RIDGE
# =============================================================================================
# Now the sensors are genuinely misplaced: each nominal position is the true one pushed off
# the curve by a normal offset of standard deviation σx. The displacement HAS to be normal —
# sliding a sensor ALONG the boundary leaves it on the boundary, where it reports the same
# value, so a tangential displacement is not an error at all. The Dirichlet condition "total
# pressure = 0" is true on the TRUE curve, so imposing it at the nominal points imposes a wrong
# condition — the realistic geometry error, and the error is still measured on the TRUE curve.
#
# What the solver does with that is decided by the covariance we hand it, because the ridge it
# produces,
#     Γ = Σᵣ wᵣ Gᵣ Σδx Gᵣᵀ,   Gᵣ = ∂Mᵣ/∂x,
# is a directional derivative penalty. Along the boundary TANGENT, Σδx = σ² t tᵀ gives
#     aᵀΓa = σ² Σᵣ wᵣ |∂(M a)ᵣ/∂s|²,
# the squared arc-length derivative of the fitted trace: a smoothness penalty along the
# boundary, which is exactly what an MFS fit that interpolates its collocation points and
# oscillates between them needs. That is the prior used here, and it is compared against the
# only honest alternative — ignoring the uncertainty altogether, which is what a deterministic
# solver does. (A NORMAL Σx penalises ∂/∂n instead. It is a real constraint, but not a smoothing
# one, and it is not run: this example asserts a tangential uncertainty throughout.)
#
# The solver's boundary prior is diagonal per sensor, so a rank-one σ² d dᵀ enters through its
# diagonal (σ²d₁², σ²d₂²) — exact when d is axis-aligned, a bounding approximation otherwise —
# plus a small isotropic floor that keeps each 2 × 2 block non-degenerate.
σx = 0.010                              # ≈ 0.35 of an inter-boundary step
prior_floor = 0.1
tangents = [SVector(-n[2], n[1]) for n in normals]

Random.seed!(99)
offsets_true = σx .* randn(N_fit)
pts_nominal = [pts[i] + offsets_true[i] * normals[i] for i in eachindex(pts)]

direction_prior(σ, dirs) = Matrix(Diagonal(vcat([
    [σ^2 * dirs[i][1]^2 + (prior_floor * σ)^2,
     σ^2 * dirs[i][2]^2 + (prior_floor * σ)^2] for i in eachindex(pts)]...)))

# A clearance is still needed: half the nominal sensors have been pushed INSIDE the star, and a
# fundamental solution is singular at its own source, so a source may not be offered a position
# on a sensor. The same figure keeps candidates far enough from the boundary that the
# linearization M(x + δx) ≈ M + Σₖ Mₖ δxₖ behind Γ is meaningful — M varies on the scale of a
# source's distance to the boundary, so a candidate closer than a few σx has no expansion.
geometry_clearance = 3 * σx
candidates_nominal = filter(candidates) do p
    depth(p) > geometry_clearance && minimum(norm(p - q) for q in pts_nominal) > geometry_clearance
end
@info "candidates usable when the boundary is uncertain" kept = length(candidates_nominal) of = length(candidates) clearance_in_steps = round(geometry_clearance / step, digits = 2)

# Pick the basis ONCE, from the nominal sensors, and then treat Σx both ways on that SAME
# basis — so the comparison is about the regularisation and nothing else. This is the two-phase
# split of the theory read literally: Phase 1 chooses which sources, Algorithm 1 (Phase 2) then
# does the inference, and it is Phase 2 that carries Σx.
v_select = solve(Simulation(medium, star_data(σ_rel; points = pts_nominal);
    source_positions = candidates_nominal, particular_solution = incident,
    solver = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 200, elbo_tol = 1e-9)))
pos_geo = copy(v_select.fsol.positions)

# Two ways to treat the same misplaced sensors, both on `pos_geo`:
#   (a) ignore   — pretend the nominal positions are exact (what a deterministic solver does)
#   (b) tangent  — Σx = σx² t tᵀ: the ridge penalises ∂/∂s, i.e. it SMOOTHS the trace
function solve_geometry(mode)
    Σx = mode === :ignore ? nothing : direction_prior(σx, tangents)
    bd = star_data(σ_rel; points = pts_nominal, Σx = Σx)
    solver = VariationalBayesianSolver(
        prior_variance = 1 ./ v_select.prior_precisions,
        select_sources_flag = false,          # the basis is already chosen
        ard_prune_flag = true,
        ard_threshold = 1e10, max_iters = 100, elbo_tol = 1e-10,
        boundary_ridge_flag = mode !== :ignore
    )
    return solve(Simulation(medium, bd;
        source_positions = pos_geo, particular_solution = incident, solver = solver))
end

v_ignore  = solve_geometry(:ignore)
v_tangent = solve_geometry(:tangent)

# the error is still measured on the TRUE boundary: that is the object the solution has to
# get right, whatever the solver was told about the sensors
err_ignore  = boundary_error(v_ignore)
err_tangent = boundary_error(v_tangent)

for (tag, v, e) in (("Σx ignored", v_ignore, err_ignore), ("Σx tangential", v_tangent, err_tangent))
    @info "misplaced sensors: $tag" sources = length(v.fsol.positions) coefficient_norm = round(norm(v.fsol.coefficients), sigdigits = 3) true_boundary_max = round(maximum(e), sigdigits = 3) rms = round(sqrt(mean(e .^ 2)), sigdigits = 3) trace_roughness = round(roughness(v), sigdigits = 3) misfit_ratio = round(v.misfit_ratio, digits = 3)
end

# =============================================================================================
# 4. THE SOURCES FOLLOW THE DATA ERROR
# =============================================================================================
# The depth of the selected sources FALLS as the noise rises, which is the opposite of what a
# "noisier data need smoother solutions" intuition predicts. The reason is conditioning, and the
# sweep measures it: DEPTH IS THE PHYSICAL PARAMETERISATION OF THE SINGULAR SPECTRUM. A deep
# source has a broad, smooth boundary trace, so deep sources are mutually near-parallel — a
# dictionary drawn only from below 8 steps has condition number ~1e9 and only 14 of its 40
# singular values above 1e-2 of the largest, against cond 15 and all 40 for a dictionary inside
# half a step. The evidence keeps adding directions until the next one is buried in the noise,
# so the condition number of what it selects tracks the data precision, and dropping the
# unresolvable directions reads geometrically as dropping the DEEP sources.
condition_number(P, bd) = (M = MethodOfFundamentalSolutions.system_matrix(P, medium, bd);
    cond(M ./ [norm(view(M, :, j)) for j in axes(M, 2)]'))

noise_levels = [0.001, 0.003, 0.01, 0.03, 0.1]
noise_runs = map(noise_levels) do σ
    bd = star_data(σ)
    v = solve(Simulation(medium, bd;
        source_positions = candidates, particular_solution = incident,
        solver = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 200, elbo_tol = 1e-9)))
    d = depth.(v.fsol.positions)
    e = boundary_error(v)
    κ = condition_number(v.fsol.positions, bd)
    @info "noise sweep" σ = σ sources = length(d) median_depth_in_steps = round(median(d) / step, digits = 2) q90_depth_in_steps = round(quantile(d, 0.9) / step, digits = 2) condition_number = round(κ, sigdigits = 3) rms_error = round(sqrt(mean(e .^ 2)), sigdigits = 3)
    (σ = σ, vsol = v, depths = d, rms = sqrt(mean(e .^ 2)), cond = κ)
end

# ...and a caveat the sweep is honest enough to print. What the evidence prefers is not what
# minimises the boundary error: restricting the SAME candidate grid to sources deeper than two
# steps beats the free choice at BOTH ends of the noise range. The evidence maximises marginal
# likelihood AT THE 200 SENSORS, and a shallow, sharply peaked column is better at explaining
# individual sensor values — noise included — than at representing the field between them.
deep_only = filter(p -> depth(p) > 2 * step, candidates)
for σ in (first(noise_levels), last(noise_levels))
    v = solve(Simulation(medium, star_data(σ); source_positions = deep_only,
        particular_solution = incident,
        solver = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 200, elbo_tol = 1e-9)))
    free = noise_runs[findfirst(r -> r.σ == σ, noise_runs)]
    @info "restricting the candidates to depth > 2 steps" σ = σ deep_only_rms = round(sqrt(mean(boundary_error(v) .^ 2)), sigdigits = 3) deep_only_sources = length(v.fsol.positions) free_choice_rms = round(free.rms, sigdigits = 3) free_choice_sources = length(free.depths)
end

# =============================================================================================
# FIGURES
# =============================================================================================
θ_dense = LinRange(0, 2π, 1200)
outline_x = [boundary_point(θ)[1] for θ in θ_dense]
outline_y = [boundary_point(θ)[2] for θ in θ_dense]
const OBSTACLE_COLOR = RGB{Float64}(0.22, 0.22, 0.25)

function star_panel(; title = "", lims = 1.15)
    plt = plot(; aspect_ratio = 1, xlims = (-lims, lims), ylims = (-lims, lims),
        axis = false, grid = false, title = title, titlefontsize = 10, legend = false)
    plot!(plt, outline_x, outline_y; lc = :black, lw = 1.2)
    return plt
end

# --- Figure 1: the three source layouts ---
p1 = star_panel(title = "Tikhonov: $(N_fit) sources, one step in")
scatter!(p1, [s[1] for s in tikhonov_sources(1.0)], [s[2] for s in tikhonov_sources(1.0)];
    mc = :steelblue, ms = 1.8, msw = 0)

p2 = star_panel(title = "Phase 1: $(length(pos_ard)) selected from $(length(candidates)) candidates")
scatter!(p2, [c[1] for c in candidates], [c[2] for c in candidates]; mc = :gray85, ms = 1.0, msw = 0)
scatter!(p2, [s[1] for s in pos_ard], [s[2] for s in pos_ard]; mc = :crimson, ms = 2.4, msw = 0)

p3 = star_panel(title = "Phase 2: $(length(pos_move)) sources, moved")
# Phase 2 both moves and prunes, so the two lists are not index-aligned; a move is bounded by a
# quarter of the source's distance to the boundary per step, so the nearest Phase-1 source is
# where each survivor came from.
for b in pos_move
    a = pos_ard[argmin([norm(b - q) for q in pos_ard])]
    plot!(p3, [a[1], b[1]], [a[2], b[2]]; lc = :gray55, lw = 0.9)
end
scatter!(p3, [s[1] for s in pos_ard], [s[2] for s in pos_ard]; mc = :gray70, ms = 1.6, msw = 0)
scatter!(p3, [s[1] for s in pos_move], [s[2] for s in pos_move]; mc = :crimson, ms = 2.4, msw = 0)

savefig(plot(p1, p2, p3; layout = (1, 3), size = (1350, 470)),
    joinpath(FIGDIR, "star_source_layouts.png"))

# --- Figure 2: the error along the whole boundary, and the Tikhonov offset sweep ---
pe = plot(; yscale = :log10, xlabel = "arc length along the boundary",
    ylabel = "|total pressure| / |incident|", legend = :topright, legendfontsize = 7,
    title = "boundary error at 2000 points, 200 of which were fitted", titlefontsize = 10,
    ylims = (1e-5, 1.0), bottom_margin = 7Plots.mm, left_margin = 5Plots.mm, top_margin = 4Plots.mm)
plot!(pe, test_arc, max.(err_half, 1e-12); lc = :orange, lw = 1.0, label = "Tikhonov, ¼ step in")
plot!(pe, test_arc, max.(err_tik, 1e-12);  lc = :steelblue, lw = 1.0,
    label = "Tikhonov, $(best_tik.offset) steps in (the oracle's offset)")
plot!(pe, test_arc, max.(err_move, 1e-12); lc = :gray35, lw = 1.0, label = "variational (Phase 1 + 2, positions moved)")
plot!(pe, test_arc, max.(err_smooth, 1e-12); lc = :crimson, lw = 1.3, label = "variational + tangential smoothing prior (positions fixed)")
plot!(pe, test_arc, max.(post_std, 1e-12); lc = :crimson, lw = 1.0, ls = :dot, label = "variational posterior s(x)")

ps = plot([r.offset for r in tik_sweep], [r.all for r in tik_sweep];
    yscale = :log10, marker = :circle, ms = 3, lc = :steelblue, mc = :steelblue,
    label = "whole boundary (max)", xlabel = "Tikhonov source offset / inter-boundary step",
    ylabel = "relative error", legend = :bottomright, legendfontsize = 7,
    title = "the residual you can see, and the error you cannot", titlefontsize = 10,
    bottom_margin = 7Plots.mm, left_margin = 5Plots.mm)
plot!(ps, [r.offset for r in tik_sweep], [r.fit for r in tik_sweep];
    marker = :diamond, ms = 3, lc = :gray50, mc = :gray50, ls = :dash,
    label = "at the fitted sensors (max)")
plot!(ps, [r.offset for r in tik_sweep], [r.rms for r in tik_sweep];
    marker = :square, ms = 3, lc = :seagreen, mc = :seagreen, label = "whole boundary (rms)")
hline!(ps, [maximum(err_move)]; lc = :crimson, lw = 1.5, label = "variational, no offset to choose")

savefig(plot(pe, ps; layout = (1, 2), size = (1350, 480)),
    joinpath(FIGDIR, "star_boundary_error.png"))

# --- Figure 3: uncertain sensor positions ---
pg = star_panel(title = "sensors misplaced by σx = $(round(σx / step, digits = 2)) steps, along the normal")
for i in 1:4:N_fit
    plot!(pg, [pts[i][1], pts_nominal[i][1]], [pts[i][2], pts_nominal[i][2]]; lc = :gray50, lw = 0.7)
end
scatter!(pg, [p[1] for p in pts_nominal], [p[2] for p in pts_nominal]; mc = :darkorange, ms = 1.8, msw = 0)

pgi = plot(; yscale = :log10, xlabel = "arc length along the TRUE boundary",
    ylabel = "|total pressure| / |incident|", legend = :topright, legendfontsize = 7,
    title = "error on the true boundary, with misplaced sensors", titlefontsize = 10,
    ylims = (1e-4, 10.0), bottom_margin = 7Plots.mm, left_margin = 5Plots.mm)
lab(tag, v) = "$tag  (‖a‖ = $(round(norm(v.fsol.coefficients), sigdigits = 3)))"
plot!(pgi, test_arc, max.(err_ignore, 1e-12);  lc = :gray40, lw = 1.0, label = lab("Σx ignored", v_ignore))
plot!(pgi, test_arc, max.(err_tangent, 1e-12); lc = :crimson, lw = 1.4, label = lab("Σx tangential", v_tangent))

savefig(plot(pg, pgi; layout = (1, 2), size = (1350, 480)),
    joinpath(FIGDIR, "star_sensor_uncertainty.png"))

# --- Figure 4: the selected sources follow the data error ---
noise_panels = map(noise_runs) do r
    p = star_panel(title = "σ = $(round(100 * r.σ, digits = 1))%  ·  $(length(r.depths)) sources")
    scatter!(p, [s[1] for s in r.vsol.fsol.positions], [s[2] for s in r.vsol.fsol.positions];
        mc = :crimson, ms = 2.4, msw = 0)
    p
end

pdepth = plot(; xscale = :log10, xlabel = "sensor noise σ / |incident|",
    ylabel = "source depth / inter-boundary step", yscale = :log10,
    legend = :topright, legendfontsize = 7, title = "how deep the sources sit", titlefontsize = 10,
    bottom_margin = 7Plots.mm, left_margin = 5Plots.mm)
plot!(pdepth, noise_levels, [median(r.depths) / step for r in noise_runs];
    marker = :circle, ms = 4, lc = :crimson, mc = :crimson, label = "median depth", lw = 2)
plot!(pdepth, noise_levels, [quantile(r.depths, 0.1) / step for r in noise_runs];
    lc = :crimson, ls = :dot, label = "10th / 90th percentile")
plot!(pdepth, noise_levels, [quantile(r.depths, 0.9) / step for r in noise_runs];
    lc = :crimson, ls = :dot, label = "")
hline!(pdepth, [1.0]; lc = :steelblue, lw = 2, label = "classical MFS: one step, always")

savefig(plot(noise_panels..., pdepth; layout = (2, 3), size = (1350, 850)),
    joinpath(FIGDIR, "star_noise_sources.png"))

# --- Figure 5: the layout is not a classical one ---
depths_move = depth.(pos_move) ./ step
ph = histogram(depths_move; bins = 30, fc = :crimson, lc = :white, lw = 0.4,
    legend = :topright, legendfontsize = 7, label = "variational ($(length(pos_move)) sources)",
    xlabel = "source depth / inter-boundary step", ylabel = "sources",
    title = "classical MFS uses one depth; the evidence uses many", titlefontsize = 10,
    bottom_margin = 6Plots.mm, left_margin = 4Plots.mm)
vline!(ph, [1.0]; lc = :steelblue, lw = 2.5, label = "classical MFS (all $(N_fit), one depth)")

# where along the boundary the sources crowd: nearest boundary arc length of each source
function nearest_arc(x)
    j = argmin([norm(x - q) for q in curve])
    return (j - 0.5) * perimeter / length(curve)
end
pc = histogram(nearest_arc.(pos_move); bins = 40, fc = :crimson, lc = :white, lw = 0.4,
    xlabel = "arc length of the nearest boundary point", ylabel = "sources",
    legend = false, title = "where along the boundary the sources sit (curvature extrema dashed)",
    titlefontsize = 10, bottom_margin = 6Plots.mm, left_margin = 4Plots.mm)
# mark the ten curvature extrema: the convex tips (5θ = π/2) and the concave valleys (5θ = 3π/2)
extremum_θs = vcat([π / 10 + 2π * k / 5 for k in 0:4], [3π / 10 + 2π * k / 5 for k in 0:4])
vline!(pc, sort(nearest_arc.(boundary_point.(extremum_θs))); lc = :black, ls = :dash, lw = 0.8)

savefig(plot(ph, pc; layout = (1, 2), size = (1350, 470)),
    joinpath(FIGDIR, "star_source_depths.png"))

# --- Figure 6: the field, and what the posterior does not know ---
xlims, ylims = (-3.6, 3.6), (-3.6, 3.6)
nx, ny = 210, 210
xg = LinRange(xlims..., nx); yg = LinRange(ylims..., ny)
vals = fill(NaN * im, ny, nx)
unc = fill(NaN, ny, nx)
for (iy, y) in enumerate(yg), (ix, x) in enumerate(xg)
    p = SVector(x, y)
    inside(p) && continue
    vals[iy, ix] = field(TractionType(), v_move, p)[1]
    s = field_std(TractionType(), v_move, p)
    unc[iy, ix] = sqrt(s[1]^2 + s[2]^2)
end

mags = [abs(vals[iy, ix]) for iy in 1:ny for ix in 1:nx
    if isfinite(abs(vals[iy, ix])) &&
        minimum(norm(SVector(xg[ix], yg[iy]) - q) for q in incident_positions) > 0.3]
cmax = quantile(mags, 0.99)
cg = cgrad(:balance)
img = map(CartesianIndices(vals)) do I
    v = vals[I]
    isfinite(abs(v)) || return OBSTACLE_COLOR
    c = get(cg, clamp((real(v) + cmax) / (2cmax), 0, 1))
    RGB{Float64}(c.r, c.g, c.b)
end

pf = plot(xg, yg, img; aspect_ratio = 1, axis = false, grid = false,
    title = "total field (posterior mean)", titlefontsize = 10, legend = false)
plot!(pf, outline_x, outline_y; lc = :black, lw = 1.2)
scatter!(pf, [s[1] for s in pos_move], [s[2] for s in pos_move]; mc = :crimson, ms = 1.8, msw = 0)
scatter!(pf, [q[1] for q in incident_positions], [q[2] for q in incident_positions];
    mc = :lime, ms = 6, msw = 1)

pu = heatmap(xg, yg, 100 .* unc ./ mean(mags); c = :viridis, aspect_ratio = 1,
    xlims = xlims, ylims = ylims, axis = false, grid = false, colorbar = true,
    title = "posterior uncertainty, % of the mean |field|", titlefontsize = 10)
plot!(pu, Shape(outline_x, outline_y); fillcolor = OBSTACLE_COLOR, lc = :black, lw = 1.2, label = "")

savefig(plot(pf, pu; layout = (1, 2), size = (1300, 640)),
    joinpath(FIGDIR, "star_field_uncertainty.png"))

# =============================================================================================
# THE NUMBERS THE README QUOTES
# =============================================================================================
@info "1. Tikhonov cannot see its own error" sensor_residual_quarter_step = round(maximum(boundary_error(fsol_half; where_pts = pts)), sigdigits = 3) whole_boundary_quarter_step = round(maximum(err_half), sigdigits = 3) sensor_residual_best = round(maximum(err_tik_fit), sigdigits = 3) whole_boundary_best = round(maximum(err_tik), sigdigits = 3) best_offset = best_tik.offset
@info "2. the like-for-like comparison, both solvers on the SAME 0.3% noisy data" tikhonov_offset = best_tik.offset tikhonov_rms = round(sqrt(mean(err_tik .^ 2)), sigdigits = 3) tikhonov_max = round(maximum(err_tik), sigdigits = 3) tikhonov_sources = N_fit variational_rms = round(sqrt(mean(err_move .^ 2)), sigdigits = 3) variational_max = round(maximum(err_move), sigdigits = 3) variational_sources = length(pos_move) rms_with_smoothing = round(sqrt(mean(err_smooth .^ 2)), sigdigits = 3)
@info "2b. how far did the Phase 2 M-step actually move the sources?" median_travel_in_steps = round(median(travel_phase2), digits = 3) mean_travel = round(mean(travel_phase2), digits = 3) q90_travel = round(quantile(travel_phase2, 0.9), digits = 3) max_travel = round(maximum(travel_phase2), digits = 3) candidate_grid_spacing_in_steps = round(grid_spacing / step, digits = 2) rms_before = round(sqrt(mean(err_ard .^ 2)), sigdigits = 3) rms_after = round(sqrt(mean(err_move .^ 2)), sigdigits = 3)
@info "3. the tangential ridge against ignoring the sensor error" ignored = round(sqrt(mean(err_ignore .^ 2)), sigdigits = 3) tangential = round(sqrt(mean(err_tangent .^ 2)), sigdigits = 3) roughness_ignored = round(roughness(v_ignore), sigdigits = 3) roughness_tangential = round(roughness(v_tangent), sigdigits = 3) coefficient_norm_ignored = round(norm(v_ignore.fsol.coefficients), sigdigits = 3) coefficient_norm_tangential = round(norm(v_tangent.fsol.coefficients), sigdigits = 3)
@info "4. the sources follow the noise" sources = [length(r.depths) for r in noise_runs] median_depth_steps = round.([median(r.depths) / step for r in noise_runs], digits = 2)
@info "5. not a classical layout" depth_decades = round(log10(maximum(depths_move) / minimum(depths_move)), digits = 2) shallowest = round(minimum(depths_move), digits = 2) deepest = round(maximum(depths_move), digits = 2) median_depth = round(median(depths_move), digits = 2) quartiles = round.(quantile(depths_move, [0.25, 0.5, 0.75, 0.95]), digits = 2) within_one_step = count(<(1), depths_move) deeper_than_five = count(>(5), depths_move)
# is the layout actually correlated with curvature, as it was on the corner version? Compare the
# density of sources near the ten curvature extrema against the rest of the boundary.
let arcs = nearest_arc.(pos_move), ex = sort(nearest_arc.(boundary_point.(extremum_θs))),
    band = 3 * step, gap(s) = minimum(min(abs(s - c), perimeter - abs(s - c)) for c in ex)
    near = count(a -> gap(a) <= band, arcs)
    @info "5b. does curvature attract sources?" fraction_of_boundary_near_an_extremum = round(10 * 2band / perimeter, digits = 3) fraction_of_sources_near_one = round(near / length(arcs), digits = 3)
end
@info "wrote figures to" FIGDIR
