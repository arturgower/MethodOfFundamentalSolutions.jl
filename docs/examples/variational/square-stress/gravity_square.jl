# # A block under its own weight: where do you put the sources when the answer is known?
#
# The `Gravity rectangle` testset of `test/benchmarks_test.jl` is one of the few MFS problems
# with a CLOSED-FORM answer. A 2L × H block of density ρ sits on the ground under gravity. Its
# stress is uniaxial and linear in height,
#
#     σ_yy = ρ g (y - H),      σ_xx = σ_xy = 0,
#
# which is an exact solution: it satisfies equilibrium with the body force (∂σ_yy/∂y = ρg
# balances f_y = -ρg), it is compatible, and it produces exactly the tractions imposed on the
# four faces — the ground pushes up with the weight per unit length W/(2L) on the bottom face,
# and the top and the two sides are traction free. So everything below can be checked against
# the truth INSIDE the block, not just against the residual on its boundary.
#
# The particular solution `ParticularGravity()` carries the body force with σ_yy = ρ g y. It is
# the right forcing but the wrong constant: it is the stress of a block whose TOP is at y = 0,
# so it leaves the MFS sources a homogeneous problem to solve,
#
#     σ_yy = -ρ g H  uniformly,
#
# a state of uniform vertical compression. That is an awkward thing to build out of point
# forces: every Kelvin source decays like 1/r, so a nearly CONSTANT stress across the block can
# only come from sources far away, or from several near ones cancelling. Where those sources
# should go is not obvious, and it is what this example is about. (`ParticularGravity(height =
# H)` would be the exact solution outright and the coefficients would come out zero — that is
# the second half of the benchmark testset, and there is nothing to learn from it.)
#
# The benchmark's own comment says it all: "need some fiddling to choose the best source
# distance and tolerance". Classical MFS gets there with `relative_source_distance = 1.24` and
# `tolerance = 2e-7`, both found by hand, both without ever consulting a quantity you would have
# in a real problem. This example hands the same data to the `VariationalBayesianSolver` and
# tunes NOTHING.
#
# THE DATA ARE EXACT, AND THE SOLVER IS TOLD THEY ARE NOISY. Every sensor reports the exact
# traction, but the boundary `fields` carry a covariance of σ = 0.5% of ρgH, so the solver
# believes there is a 0.5% measurement error and regularises accordingly. This is deliberate.
# With a noise REALISATION in the data, the last few sources go where they go partly because of
# that particular draw, and the answer to "where do sources WANT to be?" is contaminated by an
# accident of the random seed — on an earlier version of this example the selected layout came
# out visibly asymmetric although the problem is not. Deleting the realisation while keeping the
# assumed noise leaves the regularisation exactly as it was and makes the placement question
# clean: what follows is where the sources go because of the DATA, and nothing else. The price
# is that the misfit diagnostic R/N_g now reads 0.11 instead of 1 — correctly, since the
# residual really is far below the noise the solver was told about.
#
# The five things it shows:
#
#   1. CANDIDATES EVERYWHERE. A blind 41 × 41 grid over a box five times the size of the block
#      supplies 1600 candidate sources — 3200 coefficients against 160 data rows. Nobody chooses
#      an offset, a depth, or a regularisation parameter.
#   2. PHASE 1 (SELECT), run ON ITS OWN with `max_iters = 0`, starts from an EMPTY basis and adds
#      one coefficient at a time by the exact evidence criterion. It stops after 63 actions with
#      12 sources. Measured against the exact stress on a grid inside the block that is already
#      rms 0.18% against Tikhonov's 9.9%, worst point 0.47% against 38%.
#   3. PHASE 2 re-enters Algorithm 1 at the positions AND precisions Phase 1 finished with, and
#      takes the rms down another quarter to 0.14%. Most of that is the EM re-estimation of the
#      precisions rather than the moving; the M-step over the positions then slides the sources a
#      median 0.19 and at most 1.3 boundary steps off the grid and takes the worst point down a
#      further 12%. The interesting part is that THE EVIDENCE CANNOT SEE ANY OF IT: the bound
#      ends within 0.01 nat of where Phase 1 left it, because Phase 1's precisions already
#      maximize that same bound coordinate-wise. What improves is the field between the sensors,
#      which no criterion evaluated ON the sensors can measure.
#   4. THE LAYOUT IS PHYSICS, NOT AN OFFSET CURVE, AND IT IS SYMMETRIC. Nothing sits on a curve
#      one step out from the boundary. The evidence puts eight sources ON THE VERTICAL AXIS, in
#      stacks far above and far below the block — the far-field squeeze that puts a slab into
#      uniform compression — plus two symmetric lateral PAIRS, at x = ±1.125 and ±2.5, which trim
#      the σ_xx the squeeze would otherwise leave behind. The problem is symmetric about x = 0
#      and, with the noise realisation gone, so is the layout: Phase 1 selects a set that is
#      symmetric to MACHINE PRECISION. It is a property of the problem rather than of the grid,
#      too — candidate grids of 416, 1600 and 3552 points all select the same twelve. Phase 2's
#      position M-step is the one part of this that does not respect the symmetry: it descends on
#      one source at a time and drifts 0.7 steps out of symmetry in the flattest direction it
#      has, while still improving the field.
#   5. THE ANSWER IS A POSTERIOR, and it stays calibrated even though the noise it is sized for
#      is not in the data: the truth is inside the 2σ band at every point of the x = 0 slice, at
#      a median z-score of 0.78. That is less of a coincidence than it looks. The assumed σ is
#      what tells the evidence when to stop adding sources, so it sets the accuracy the fit
#      settles at — and it also sets the width of the posterior. Both come out near a fifth of
#      the assumed 0.5%, and they track each other because they have the same cause.
#
# Needs Plots.jl (not a dependency of this package): run from an environment where both
# `using MethodOfFundamentalSolutions` and `using Plots` work.

using MethodOfFundamentalSolutions
using Distributions, LinearAlgebra, StaticArrays, Statistics
using Plots
gr()

# the figures are written next to this script
const FIGDIR = @__DIR__

# ---------------------------------------------------------------------------------------------
# The block, its boundary sensors, and the exact stress — all exactly as the benchmark sets them
# ---------------------------------------------------------------------------------------------
medium = Elastostatic(2; ρ = 1.0, cp = 3.0, cs = 2.0)
g_acc = 9.81

n = 20; L = 0.5; H = 0.5          # the block is 2L wide and H tall, resting on y = 0
xb = [LinRange(-L, L, n+2)[2:end-1]; zeros(n) .+ L; LinRange(L, -L, n+2)[2:end-1]; zeros(n) .- L]
yb = [zeros(n); LinRange(0, H, n+2)[2:end-1]; zeros(n) .+ H; LinRange(H, 0, n+2)[2:end-1]]

points = [SVector(xb[i], yb[i]) for i in eachindex(xb)]
normals = [
    [SVector( 0.0, -1.0) for _ in 1:n];      # bottom face, on the ground
    [SVector( 1.0,  0.0) for _ in 1:n];      # right face, free
    [SVector( 0.0,  1.0) for _ in 1:n];      # top face, free
    [SVector(-1.0,  0.0) for _ in 1:n]       # left face, free
]

weight = 2L * H * medium.ρ * g_acc
σ_scale = weight / (2L)                      # = ρ g H, the largest stress in the block
step = 2L / (n + 1)                          # sensor spacing along the top and bottom faces

# how far a source sits from the block, the yardstick every distance below is quoted in
depth(p) = minimum(norm(p - q) for q in points)

# the exact solution, used ONLY to generate the data and to score the answers
σyy_exact(y) = medium.ρ * g_acc * (y - H)
traction_exact(p, nv) = SVector(0.0, σyy_exact(p[2]) * nv[2])

# ---------------------------------------------------------------------------------------------
# The data: the exact tractions, with an ASSUMED measurement noise
# ---------------------------------------------------------------------------------------------
# The MEAN of each `MvNormal` is the exact traction — no realisation is drawn — and its
# covariance is the noise the solver is asked to assume. That covariance is the only place a
# noise level is ever specified, and it is what replaces Tikhonov's hand-chosen `tolerance`: it
# sets how hard the evidence pushes to explain the last fraction of a percent of the data, and
# so how many sources get selected. Passing exact values with a nonzero assumed σ is therefore
# not a contradiction — it says "regularise as if the data were 0.5% noisy", which is exactly
# the knob, with none of the placement noise that a realisation would inject.
σ_noise = 0.005 * σ_scale                    # 0.5% of the largest boundary traction
tractions = [traction_exact(points[i], normals[i]) for i in eachindex(points)]

bd = BoundaryData(TractionType(); boundary_points = points, normals = normals,
    fields = [MvNormal(t, σ_noise^2 * I(2)) for t in tractions])

# the same exact numbers, as plain vectors, for the deterministic solver
bd_tik = BoundaryData(TractionType(); boundary_points = points, normals = normals, fields = tractions)

# ---------------------------------------------------------------------------------------------
# How every answer below is scored: the whole stress tensor, inside the block
# ---------------------------------------------------------------------------------------------
# A traction is σ·n, so evaluating the field with n = (0,1) returns (σ_xy, σ_yy) and with
# n = (1,0) returns (σ_xx, σ_xy). Two evaluations therefore give the complete stress state, and
# the error is its distance from the exact one, relative to ρ g H. The test grid is inset by a
# small margin: MFS is a boundary method and the corners of a rectangle are where any of it is
# worst, so sitting a test point exactly on one measures the corner and nothing else.
margin = 0.02
test_pts = vec([SVector(x, y) for x in range(-L + margin, L - margin, length = 25),
                                  y in range(margin, H - margin, length = 25)])
slice_ys = range(0, H, length = 200)
slice_pts = [SVector(0.0, y) for y in slice_ys]      # the benchmark's own slice, x = 0

stress(sol, p) = begin
    fy = field(TractionType(), sol, p, SVector(0.0, 1.0))     # (σ_xy, σ_yy)
    fx = field(TractionType(), sol, p, SVector(1.0, 0.0))     # (σ_xx, σ_xy)
    SVector(fx[1], fy[1], fy[2])                              # (σ_xx, σ_xy, σ_yy)
end
stress_exact(p) = SVector(0.0, 0.0, σyy_exact(p[2]))
stress_error(sol, pts = test_pts) = [norm(stress(sol, p) - stress_exact(p)) / σ_scale for p in pts]

# ---------------------------------------------------------------------------------------------
# 1. CLASSICAL MFS, with the offset and the tolerance the benchmark had to find by hand
# ---------------------------------------------------------------------------------------------
sources_tik = MethodOfFundamentalSolutions.source_positions(bd_tik; relative_source_distance = 1.24)
fsol_tik = solve(Simulation(medium, bd_tik;
    particular_solution = ParticularGravity(),
    solver = TikhonovSolver(tolerance = 2e-7),
    source_positions = sources_tik
))
err_tik = stress_error(fsol_tik)

@info "Tikhonov, offset 1.24 steps and tolerance 2e-7, both hand-tuned" sources = length(sources_tik) rms = round(sqrt(mean(err_tik .^ 2)), sigdigits = 3) worst = round(maximum(err_tik), sigdigits = 3) on_the_x0_slice = round(sqrt(mean(stress_error(fsol_tik, slice_pts) .^ 2)), sigdigits = 3)

# ---------------------------------------------------------------------------------------------
# 2. PHASE 1 (SELECT): candidates everywhere, and let the evidence pick
# ---------------------------------------------------------------------------------------------
# `grid_source_positions` lays a blind regular grid over a box `scale` times the size of the
# block and keeps every grid point OUTSIDE it, at least `clearance` sensor spacings clear of the
# boundary. Sources must be outside because each is singular at its own position and the stress
# must be regular in the block. The grid is deliberately FINE — 41 × 41, so its spacing is 2.6
# boundary steps in x and 1.2 in y — because Phase 1 can only ever choose from what it is
# offered, and the point of the example is to see where the sources want to be rather than where
# a coarse grid allows them to be. Section 2b checks both of the numbers set here.
candidates = grid_source_positions(bd; n = 41, scale = 5.0, clearance = 1.0)

@info "an overcomplete basis, chosen blind" candidates = length(candidates) coefficients = 2 * length(candidates) data_rows = 2 * length(points)

# Phase 1 starts from an EMPTY basis and repeatedly applies the single action — add a
# coefficient at its exact optimal precision, re-estimate one, or delete one — that most
# increases the evidence, stopping when none of them does.
#
# `max_iters = 0` is what makes this Phase 1 AND NOTHING ELSE. A default
# `VariationalBayesianSolver` call runs both phases back to back: it selects, and then falls
# straight into the Algorithm 1 loop on the survivors, so its `elbo_history` is the two phases
# concatenated and its answer is not the selection's. Setting the Phase-2 iteration budget to
# zero stops it after the selection, at the exact per-column evidence maximizers Phase 1 found,
# which is what section 3 then picks up.
solver_select = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 0, elbo_tol = 1e-9)
t_select = @elapsed v_select = solve(Simulation(medium, bd;
    source_positions = candidates, particular_solution = ParticularGravity(), solver = solver_select))

pos_select = copy(v_select.fsol.positions)
err_select = stress_error(v_select)

@info "Phase 1 (select)" kept = length(pos_select) of = length(candidates) seconds = round(t_select, digits = 1) actions = length(v_select.elbo_history) misfit_ratio = round(v_select.misfit_ratio, digits = 3) rms = round(sqrt(mean(err_select .^ 2)), sigdigits = 3) worst = round(maximum(err_select), sigdigits = 3)

# ---------------------------------------------------------------------------------------------
# 2b. THE TWO NUMBERS THE CANDIDATE GRID IS MADE OF
# ---------------------------------------------------------------------------------------------
# Only `scale` and `n` were chosen above, so both deserve a measurement rather than an assertion.
#
# HOW FAR OUT the box reaches is the one that matters. A box only 1.3 times the block keeps the
# deepest candidate 3.4 steps out and spends 59 sources to land SIX TIMES the error; opening it
# to 5 times the block gets a better answer out of 12. A nearly uniform stress is cheap to build
# from far sources and expensive to build from near ones, which have to cancel across the block —
# and with clean data that cost is visible in the error, where a noise realisation would hide it
# under the noise floor.
for sc in (1.3, 2.0, 3.0, 5.0, 8.0)
    cs = grid_source_positions(bd; n = 41, scale = sc, clearance = 1.0)
    v = solve(Simulation(medium, bd; source_positions = cs, particular_solution = ParticularGravity(),
        solver = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 0, elbo_tol = 1e-9)))
    e = stress_error(v)
    @info "candidate box" scale = sc candidates = length(cs) kept = length(v.fsol.positions) deepest_kept_in_steps = round(maximum(depth.(v.fsol.positions)) / step, digits = 1) rms = round(sqrt(mean(e .^ 2)), sigdigits = 3) worst = round(maximum(e), sigdigits = 3)
end

# HOW FINE the grid is barely matters at all, which is the more interesting result: 416, 1600 and
# 3552 candidates all select the SAME TWELVE sources. Phase 1 is not sampling a continuum, it is
# finding a set of isolated evidence maxima, and once the grid resolves them a finer one has
# nothing new to offer. (At 6244 candidates it starts keeping 19: the grid gets fine enough that
# neighbouring candidates are nearly interchangeable, and the selection splits single sources
# across pairs of them without improving the field.)
for gn in (21, 41, 61, 81)
    cs = grid_source_positions(bd; n = gn, scale = 5.0, clearance = 1.0)
    t = @elapsed v = solve(Simulation(medium, bd; source_positions = cs, particular_solution = ParticularGravity(),
        solver = VariationalBayesianSolver(ard_threshold = 1e10, max_iters = 0, elbo_tol = 1e-9)))
    e = stress_error(v)
    @info "candidate grid" n = gn candidates = length(cs) kept = length(v.fsol.positions) seconds = round(t, digits = 2) rms = round(sqrt(mean(e .^ 2)), sigdigits = 3)
end

# ---------------------------------------------------------------------------------------------
# 3. PHASE 2 (MOVE): Algorithm 1 on the survivors
# ---------------------------------------------------------------------------------------------
# Entered exactly as the theory requires — at the positions AND the precisions Phase 1 finished
# with. `select_sources_flag = false` skips a second selection, `prior_variance` hands the
# learned precisions back, and now the source positions are hyperparameters too, moved by
# gradient steps that are accepted only when they increase the bound.
solver_move = VariationalBayesianSolver(
    prior_variance = 1 ./ v_select.prior_precisions,
    select_sources_flag = false,           # Phase 1 is done: do not select again
    ard_prune_flag = true,                 # but Algorithm 1 line 8 still prunes
    optimise_source_positions_flag = true,
    source_position_iters = 5,
    ard_threshold = 1e10, max_iters = 100, elbo_tol = 1e-9
)
t_move = @elapsed v_move = solve(Simulation(medium, bd;
    source_positions = pos_select, particular_solution = ParticularGravity(), solver = solver_move))

pos_move = copy(v_move.fsol.positions)
err_move = stress_error(v_move)

# How much of that is the POSITIONS? The same Phase 2 with the M-step over χ switched off
# re-estimates the precisions and nothing else, so the difference between the two is what moving
# the sources bought.
v_fixed = solve(Simulation(medium, bd; source_positions = pos_select,
    particular_solution = ParticularGravity(),
    solver = VariationalBayesianSolver(
        prior_variance = 1 ./ v_select.prior_precisions,
        select_sources_flag = false, ard_prune_flag = true,
        ard_threshold = 1e10, max_iters = 100, elbo_tol = 1e-9)))
err_fixed = stress_error(v_fixed)

# Phase 2 both moves and prunes, so the two lists are not index-aligned; a source can move at
# most a quarter of its distance to the boundary per accepted step, so the nearest Phase-1
# position is where each survivor came from.
travel = [minimum(norm(b - a) for a in pos_select) / step for b in pos_move]

@info "Phase 2 (move)" sources = length(pos_move) seconds = round(t_move, digits = 1) iterations = length(v_move.elbo_history) misfit_ratio = round(v_move.misfit_ratio, digits = 3) rms = round(sqrt(mean(err_move .^ 2)), sigdigits = 3) worst = round(maximum(err_move), sigdigits = 3) evidence_gained = round(v_move.elbo_history[end] - v_select.elbo_history[end], digits = 3) median_travel_in_steps = round(median(travel), digits = 2) max_travel_in_steps = round(maximum(travel), digits = 2)

# Two things worth reading off this. Most of Phase 2's gain is the EM re-estimation of the
# precisions, not the moving: holding the positions fixed already recovers nearly all of it, and
# the M-step over χ adds the rest. And the EVIDENCE cannot see any of it — the bound ends within
# 0.01 nat of where Phase 1 left it, while the true error falls by a quarter. Phase 1's per-column
# precisions are the exact coordinate maximizers of the same bound, so there is very little left
# for Phase 2 to win on the criterion it optimises; what it improves is the field BETWEEN the
# sensors, which no boundary-data criterion measures.
@info "what Phase 2 actually bought" rms_phase1 = round(sqrt(mean(err_select .^ 2)), sigdigits = 3) rms_precisions_only = round(sqrt(mean(err_fixed .^ 2)), sigdigits = 3) rms_precisions_and_positions = round(sqrt(mean(err_move .^ 2)), sigdigits = 3) worst_phase1 = round(maximum(err_select), sigdigits = 3) worst_precisions_only = round(maximum(err_fixed), sigdigits = 3) worst_precisions_and_positions = round(maximum(err_move), sigdigits = 3)

# ---------------------------------------------------------------------------------------------
# 4. WHERE THE EVIDENCE PUT THE SOURCES
# ---------------------------------------------------------------------------------------------
# Not on an offset curve. Sorted by how much each source actually contributes, the layout reads
# as a physical construction: eight sources on the axis of symmetry, stacked far ABOVE and far
# BELOW the block — the far-field squeeze that puts a slab into uniform vertical compression —
# and two lateral PAIRS at x = ±1.125 and ±2.5 that trim the σ_xx such a squeeze leaves behind.
# Note how far out they all are: at that range a 1/r field is nearly constant across the block,
# which is exactly the state that has to be represented.
amp(sol, j) = norm(sol.fsol.coefficients[2j-1:2j])

order = sortperm([amp(v_move, j) for j in eachindex(pos_move)]; rev = true)
for j in order
    @info "a source, in order of amplitude" at = round.(Tuple(pos_move[j]), digits = 3) depth_in_steps = round(depth(pos_move[j]) / step, digits = 1) amplitude = round(amp(v_move, j), sigdigits = 3)
end
@info "the layout" sources = length(pos_move) depth_range_in_steps = round.((minimum(depth.(pos_move)), maximum(depth.(pos_move))) ./ step, digits = 1) median_depth_in_steps = round(median(depth.(pos_move)) / step, digits = 1) classical_MFS_uses = 1.24

# --- Is it symmetric? The block, the loading and the assumed noise are all symmetric about
#     x = 0, so the ideal layout must be too — but nothing in the algorithm imposes that, and
#     with a noise realisation in the data it came out visibly broken. Measure it: for every
#     source, the distance from its mirror image to the nearest other source.
#
#     Phase 1 scores EXACTLY zero, which is worth pausing on: the selection is a sequence of
#     greedy discrete decisions with no notion of symmetry in it, and it recovers the symmetry of
#     the problem anyway, because with clean data two mirror-image candidates have identical
#     evidence and both get taken. Phase 2's position M-step then loses some of it. That is a
#     limitation of the M-step rather than of the model: it descends on one source at a time, so
#     a symmetric pair is never moved as a pair, and the lateral pair sits in the flattest
#     direction the misfit has — moving one of them barely changes the field, which is exactly
#     why the descent is free to.
asymmetry(pos) = maximum(minimum(norm(SVector(-p[1], p[2]) - q) for q in pos) for p in pos) / step
@info "symmetry about x = 0, in boundary steps" phase1 = round(asymmetry(pos_select), sigdigits = 3) phase2 = round(asymmetry(pos_move), sigdigits = 3)

# ---------------------------------------------------------------------------------------------
# 5. THE ANSWER IS A POSTERIOR: what does its error bar mean when the data are exact?
# ---------------------------------------------------------------------------------------------
# `field_std` propagates the coefficient posterior to the field, so every point of the block
# carries a standard deviation on its stress. That standard deviation is sized for the assumed
# 0.5% sensor noise, and the data do not contain it, so there is no reason in principle for it to
# match the error actually made — and yet it does, to within 25%. The assumed σ is the quantity
# that tells the evidence when adding another source stops being worth it, so it sets the
# accuracy the fit settles at; it also sets the width of the posterior. Both land near a fifth of
# it. The band would only be badly wrong if the two mechanisms came apart — which is what a
# geometry error, an unmodelled load or a boundary the sources cannot represent would do.
z_slice = map(slice_pts) do p
    s = field_std(TractionType(), v_move, p, SVector(0.0, 1.0))[2]
    (abs(stress(v_move, p)[3] - σyy_exact(p[2])) / s, s)
end
@info "posterior calibration on the x = 0 slice" median_z_score = round(median(first.(z_slice)), digits = 2) fraction_inside_2σ = round(mean(first.(z_slice) .< 2), digits = 3) median_posterior_std = round(median(last.(z_slice)) / σ_scale, sigdigits = 3) median_true_error = round(median([abs(stress(v_move, p)[3] - σyy_exact(p[2])) / σ_scale for p in slice_pts]), sigdigits = 3)

# =============================================================================================
# FIGURES
# =============================================================================================
block_x = [-L, L, L, -L, -L]
block_y = [0.0, 0.0, H, H, 0.0]

# --- Figure 1: the total predicted field against the true solution -----------------------------
hx = range(-L, L, length = 121)
hy = range(0, H, length = 61)
Z_true = [σyy_exact(y) for y in hy, x in hx]
Z_var  = [stress(v_move, SVector(x, y))[3] for y in hy, x in hx]
Z_tik  = [stress(fsol_tik, SVector(x, y))[3] for y in hy, x in hx]
E_var  = [norm(stress(v_move, SVector(x, y)) - stress_exact(SVector(x, y))) / σ_scale for y in hy, x in hx]
E_tik  = [norm(stress(fsol_tik, SVector(x, y)) - stress_exact(SVector(x, y))) / σ_scale for y in hy, x in hx]

# what the posterior claims about itself, over the whole stress state
S_var = [begin
        sy = field_std(TractionType(), v_move, SVector(x, y), SVector(0.0, 1.0))
        sx = field_std(TractionType(), v_move, SVector(x, y), SVector(1.0, 0.0))
        norm(SVector(sx[1], sy[1], sy[2])) / σ_scale
    end for y in hy, x in hx]

σlims = (minimum(Z_true), 0.0)
# the errors span two decades between the two methods, so they share a LOG colour scale: on a
# common linear one the variational panel is simply black, which says "smaller" and nothing else
elims = (-1.5, 1.5)
logpct(E) = log10.(clamp.(100 .* E, 10.0^first(elims), 10.0^last(elims)))

function block_map(Z; title, clims, cmap, cbtitle)
    p = heatmap(hx, hy, Z; aspect_ratio = 1, c = cmap, clims = clims,
        title = title, titlefontsize = 9, xlabel = "x", ylabel = "y",
        colorbar_title = cbtitle, colorbar_titlefontsize = 8,
        xlims = (-L - 0.03, L + 0.03), ylims = (-0.03, H + 0.03),
        left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
    plot!(p, block_x, block_y; lc = :black, lw = 1.2, label = "")
    return p
end
err_map(E; title) = block_map(logpct(E); title = title, clims = elims, cmap = :magma,
    cbtitle = "log₁₀(% of ρgH)")

f1 = block_map(Z_true; title = "exact:  σ_yy = ρ g (y - H)", clims = σlims, cmap = :viridis, cbtitle = "σ_yy")
f2 = block_map(Z_var;  title = "variational total field, $(length(pos_move)) sources, nothing tuned", clims = σlims, cmap = :viridis, cbtitle = "σ_yy")
f3 = block_map(Z_tik;  title = "Tikhonov total field, $(length(sources_tik)) sources, offset and tolerance tuned", clims = σlims, cmap = :viridis, cbtitle = "σ_yy")

f4 = err_map(E_var; title = "variational error: rms $(round(100*sqrt(mean(err_move.^2)), digits=2))%, worst $(round(100*maximum(err_move), digits=1))%")
f5 = err_map(E_tik; title = "Tikhonov error: rms $(round(100*sqrt(mean(err_tik.^2)), digits=1))%, worst $(round(100*maximum(err_tik), digits=1))%")
f6 = err_map(S_var; title = "the error the variational solver CLAIMS: posterior std")

savefig(plot(f1, f2, f3, f4, f5, f6; layout = (2, 3), size = (1500, 600)),
    joinpath(FIGDIR, "gravity_square_field.png"))

# --- Figure 2: the x = 0 slice, with the posterior's own error bars ----------------------------
σyy_var = [stress(v_move, p)[3] for p in slice_pts]
σyy_tik = [stress(fsol_tik, p)[3] for p in slice_pts]
s_var   = [field_std(TractionType(), v_move, p, SVector(0.0, 1.0))[2] for p in slice_pts]

ps = plot(slice_ys, [σyy_exact(y) for y in slice_ys]; lc = :black, lw = 2.5, label = "exact  ρ g (y - H)",
    xlabel = "height y", ylabel = "σ_yy", legend = :bottomright, legendfontsize = 8,
    title = "the total predicted field on x = 0", titlefontsize = 10,
    bottom_margin = 6Plots.mm, left_margin = 5Plots.mm)
plot!(ps, slice_ys, σyy_var; ribbon = 2 .* s_var, fc = :crimson, fa = 0.25,
    lc = :crimson, lw = 1.8, label = "variational, Phase 1 + 2  (±2σ posterior)")
plot!(ps, slice_ys, σyy_tik; lc = :steelblue, lw = 1.4, ls = :dash, label = "Tikhonov, hand-tuned")

# the same slice with the exact solution subtracted: at 1% of ρgH the three curves above lie on
# top of one another, and only the residual shows what each method actually did
pd = plot(slice_ys, 100 .* (σyy_var .- [σyy_exact(y) for y in slice_ys]) ./ σ_scale;
    ribbon = 200 .* s_var ./ σ_scale, fc = :crimson, fa = 0.25, lc = :crimson, lw = 1.8,
    label = "variational  (±2σ posterior)", xlabel = "height y", ylabel = "% of ρgH",
    legend = :topleft, legendfontsize = 8, title = "the same slice, minus the exact solution",
    titlefontsize = 10, bottom_margin = 6Plots.mm, left_margin = 5Plots.mm)
plot!(pd, slice_ys, 100 .* (σyy_tik .- [σyy_exact(y) for y in slice_ys]) ./ σ_scale;
    lc = :steelblue, lw = 1.6, ls = :dash, label = "Tikhonov, hand-tuned")
hline!(pd, [0.0]; lc = :black, lw = 1.0, label = "")

# the two components the truth says are zero: a pure error channel
pz = plot(slice_ys, 100 .* [stress(v_move, p)[1] for p in slice_pts] ./ σ_scale;
    lc = :crimson, lw = 1.8, label = "variational σ_xx", xlabel = "height y",
    ylabel = "% of ρgH", legend = :topright, legendfontsize = 8,
    title = "the components that should vanish", titlefontsize = 10,
    bottom_margin = 6Plots.mm, left_margin = 5Plots.mm)
plot!(pz, slice_ys, 100 .* [stress(v_move, p)[2] for p in slice_pts] ./ σ_scale;
    lc = :crimson, lw = 1.2, ls = :dot, label = "variational σ_xy")
plot!(pz, slice_ys, 100 .* [stress(fsol_tik, p)[1] for p in slice_pts] ./ σ_scale;
    lc = :steelblue, lw = 1.4, label = "Tikhonov σ_xx")
plot!(pz, slice_ys, 100 .* [stress(fsol_tik, p)[2] for p in slice_pts] ./ σ_scale;
    lc = :steelblue, lw = 1.0, ls = :dot, label = "Tikhonov σ_xy")
hline!(pz, [0.0]; lc = :black, lw = 1.0, label = "")

savefig(plot(ps, pd, pz; layout = (1, 3), size = (1500, 450)),
    joinpath(FIGDIR, "gravity_square_slice.png"))

# --- Figure 3: the three source layouts --------------------------------------------------------
lims_x = (-2.8, 2.8); lims_y = (-1.3, 1.8)

function source_panel(; title = "")
    p = plot(; aspect_ratio = 1, xlims = lims_x, ylims = lims_y, grid = false,
        title = title, titlefontsize = 10, legend = false, xlabel = "x", ylabel = "y",
        left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
    plot!(p, block_x, block_y; lc = :black, lw = 1.6)
    return p
end

q1 = source_panel(title = "Tikhonov: $(length(sources_tik)) sources, 1.24 steps out")
scatter!(q1, [s[1] for s in sources_tik], [s[2] for s in sources_tik]; mc = :steelblue, ms = 2.2, msw = 0)

q2 = source_panel(title = "Phase 1: $(length(pos_select)) selected from $(length(candidates)) candidates")
scatter!(q2, [c[1] for c in candidates], [c[2] for c in candidates]; mc = :gray85, ms = 1.6, msw = 0)
scatter!(q2, [s[1] for s in pos_select], [s[2] for s in pos_select]; mc = :crimson, ms = 4.0, msw = 0)

# marker area ∝ log amplitude, so the panel shows WHICH sources carry the solution
amps = [amp(v_move, j) for j in eachindex(pos_move)]
msize = 2.5 .+ 6 .* (log10.(amps) .- log10(minimum(amps))) ./ (log10(maximum(amps)) - log10(minimum(amps)))
q3 = source_panel(title = "Phase 2: $(length(pos_move)) sources, sized by amplitude")
for b in pos_move
    a = pos_select[argmin([norm(b - q) for q in pos_select])]
    plot!(q3, [a[1], b[1]], [a[2], b[2]]; lc = :gray55, lw = 0.9)
end
scatter!(q3, [s[1] for s in pos_select], [s[2] for s in pos_select]; mc = :gray75, ms = 2.4, msw = 0)
scatter!(q3, [s[1] for s in pos_move], [s[2] for s in pos_move]; mc = :crimson, ms = msize, msw = 0)

savefig(plot(q1, q2, q3; layout = (1, 3), size = (1400, 330)),
    joinpath(FIGDIR, "gravity_square_sources.png"))

# --- Figure 4: the bound, and the boundary residual --------------------------------------------
# F itself is useless to plot: it climbs from -1.5e5 to O(1) in the first few actions and the
# rest of the run is a flat line. What the two phases are actually doing shows up in how far the
# bound still is from where it ends, on a log scale.
n1 = length(v_select.elbo_history)
F_end = max(v_move.elbo_history[end], maximum(v_select.elbo_history))
gap(h) = max.(F_end .- h, 3e-3)          # floored, so the converged tail does not go to -Inf
pf = plot(1:n1, gap(v_select.elbo_history); lc = :crimson, lw = 1.8, yscale = :log10,
    label = "Phase 1 (one action each)", xlabel = "action / iteration",
    ylabel = "F(best) - F", legend = :topright, legendfontsize = 8,
    title = "how far the bound still has to climb", titlefontsize = 10,
    bottom_margin = 6Plots.mm, left_margin = 8Plots.mm)
plot!(pf, n1 .+ (1:length(v_move.elbo_history)), gap(v_move.elbo_history);
    lc = :steelblue, lw = 1.8, label = "Phase 2 (variational EM)")
vline!(pf, [n1]; lc = :gray50, ls = :dash, lw = 1.0, label = "")
# every iteration at which Phase 2 pruned a source: the model changed there, so F is not
# comparable across it (`baseline_resets` is exactly this list)
vline!(pf, n1 .+ v_move.baseline_resets; lc = :seagreen, ls = :dot, lw = 1.0,
    label = "Phase 2 pruned a source")

resid(sol) = [norm(field(TractionType(), sol, points[i], normals[i]) - tractions[i]) for i in eachindex(points)]
pr = plot(1:length(points), resid(v_move) ./ σ_scale; lc = :crimson, lw = 1.4,
    label = "variational, Phase 1 + 2", xlabel = "boundary sensor", ylabel = "residual / ρgH",
    legend = :topright, legendfontsize = 8, yscale = :log10, ylims = (1e-5, 1e0),
    title = "the traction residual at the sensors", titlefontsize = 10,
    bottom_margin = 6Plots.mm, left_margin = 5Plots.mm)
plot!(pr, 1:length(points), resid(fsol_tik) ./ σ_scale; lc = :steelblue, lw = 1.2, label = "Tikhonov, hand-tuned")
hline!(pr, [sqrt(2) * σ_noise / σ_scale]; lc = :black, ls = :dash, lw = 1.5, label = "the noise the solver was told to assume")

savefig(plot(pf, pr; layout = (1, 2), size = (1250, 470)),
    joinpath(FIGDIR, "gravity_square_evidence.png"))

@info "figures written" dir = FIGDIR
