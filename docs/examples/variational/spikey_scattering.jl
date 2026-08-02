# # Scattering from spikey obstacles: select the sources once, reuse them for any incidence
#
# This example advertises what the variational method (`VariationalBayesianSolver`) buys you
# over classical MFS on a problem where classical MFS is genuinely awkward: scattering from
# THREE strange obstacles covered in spikes, each carrying a Dirichlet (sound-soft)
# condition — the total pressure vanishes on every boundary. Dirichlet obstacles are strong
# scatterers, and the three sit only about half a wavelength apart, so the field visibly
# rattles between them: the multiple scattering shows up as interference structure in the
# gaps between the obstacles. With classical MFS, source placement near each spike is
# make-or-break (see docs/examples/acoustic/teardrop_scattering.jl, where a single missing
# source at a cusp destabilises the whole fit). Here we refuse to think about placement at
# all, and we go one step further — we select the sources ONCE, for many incidences at once:
#
#   * candidate sources are laid down EVERYWHERE the physics allows — a blind regular grid
#     carpeting the inside of every obstacle, with no idea where the spikes are;
#   * sparse Bayesian learning then SELECTS from that grid one source at a time: starting from
#     an empty basis it repeatedly adds the single candidate that most increases the Bayesian
#     evidence, each at its exact optimal prior precision, and stops when no candidate is worth
#     adding — but it is TRAINED on an omni-directional illumination: a ring of incident point
#     sources all around the scene, whose fields are summed (a `Simulation` now accepts a vector
#     of particular solutions). The selected sources therefore have to resolve the whole of
#     every boundary, not just the lit side of a single incidence;
#   * we then FIX those selected positions and re-solve for a SINGLE incident source. Because
#     the solver is Bayesian, that answer is a posterior: every point of the scattered field
#     comes with an uncertainty s(x), mapped in a companion panel as a percentage of the mean
#     field magnitude — it concentrates in the gaps between the obstacles, where the
#     multiply-scattered field is hardest to pin down. It is this single-source posterior that
#     the field picture and the gif animate.
#
# (Sources must lie inside the obstacles: the scattered field must be regular and radiating
# everywhere in the exterior domain, and each MFS source is singular at its own position.
# Inside the obstacles, though, the grid is completely uninformed — that is the point.)
#
# The Dirichlet condition (total pressure = 0 on each boundary) is imposed as noisy data:
# each boundary sensor reports zero total pressure to within a known noise σ. The incident
# field is a point source (or, for training, a ring of them), passed as the particular
# solution, so the solved `FundamentalSolution` evaluates to the total field.
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
# Three spikey star-shaped obstacles, each a polar curve r(θ) around its own centre.
# `tri` is a triangle wave, so the star and the shard have genuinely sharp (corner) spikes;
# the urchin's spikes are smooth but so narrow that they are just as hard to resolve.
# ---------------------------------------------------------------------------------------------
tri(x) = (2 / π) * asin(sin(x))

# The centres put the obstacles only about half a wavelength apart (the wavelength is 1,
# set below): Dirichlet obstacles are strong scatterers, so at this separation the wave
# rattles between them and the multiple scattering is visible in the gaps.
obstacles = [
    (name = "star",   centre = SVector(-1.30, 0.95),
        r = θ -> 0.75 * (1 + 0.35 * tri(5θ))),
    (name = "urchin", centre = SVector(1.45, 0.90),
        r = θ -> 0.55 * (1 + 0.85 * exp(2.5 * (cos(5θ) - 1)))),
    (name = "shard",  centre = SVector(0.05, -1.35),
        r = θ -> 0.70 * (1 + 0.25 * tri(4θ + 0.5) + 0.15 * cos(2θ + 0.8) + 0.10 * tri(7θ))),
]

boundary_point(o, θ) = o.centre + o.r(θ) * SVector(cos(θ), sin(θ))

# outward normal from the finite-difference tangent, oriented away from the obstacle centre
function boundary_normal(o, θ; h = 1e-4)
    t = boundary_point(o, θ + h) - boundary_point(o, θ - h)
    n = SVector(t[2], -t[1]) / norm(SVector(t[2], -t[1]))
    return dot(n, boundary_point(o, θ) - o.centre) < 0 ? -n : n
end

# the obstacles are star-shaped, so "inside" is exact and cheap from the polar radius
inside(o, x) = norm(x - o.centre) < o.r(atan(x[2] - o.centre[2], x[1] - o.centre[1]))
inside_any(x) = any(inside(o, x) for o in obstacles)

# N boundary sensors per obstacle, equally spaced in ARC LENGTH. Uniform θ would be a trap
# here: on the steep spike flanks the curve races through a lot of arc per unit θ, so
# uniform-θ sensors leave exactly the hard parts of the boundary under-sampled. The half-step
# offset keeps every sensor off the triangle-wave corners, where the finite-difference
# normal would be ambiguous.
function arc_length_thetas(o, N; dense = 4000)
    θd = LinRange(0, 2π, dense + 1)
    seg = [norm(boundary_point(o, θd[i + 1]) - boundary_point(o, θd[i])) for i in 1:dense]
    cum = cumsum(seg)
    targets = ((1:N) .- 0.5) .* (cum[end] / N)
    return [θd[searchsortedfirst(cum, s)] for s in targets]
end

N_per = 200
θs_per = [arc_length_thetas(o, N_per) for o in obstacles]
pts     = vcat([[boundary_point(o, θ) for θ in θs] for (o, θs) in zip(obstacles, θs_per)]...)
normals = vcat([[boundary_normal(o, θ) for θ in θs] for (o, θs) in zip(obstacles, θs_per)]...)
centres = [o.centre for o in obstacles]

# ---------------------------------------------------------------------------------------------
# Physics: Dirichlet (sound-soft) obstacles. The MFS sources are SELECTED with ARD from an
# omni-directional training illumination — a ring of incident point sources all around the
# scene — and only then reused for a single incidence. The boundary data are noisy
# measurements of zero total pressure: each sensor is an `MvNormal` over the stacked [Re; Im]
# parts of its field, with the known noise σ as its covariance — this is how the variational
# solver is told the measurement noise.
# ---------------------------------------------------------------------------------------------
ω = 2π                       # wavelength 2π/ω = 1: obstacle diameters span a few wavelengths
medium = Acoustic(2; ω = ω, ρ = 1.0, c = 1.0)

# the trace of an incident field on the boundary sensors, as one complex scalar per sensor.
# `inc` may be a single `PointSource` or a vector of them, in which case their fields are
# summed (see the vector `particular_solution` of the training simulation below).
trace_on_boundary(inc) = [field(TractionType(), medium, inc, p, n)[1] for (p, n) in zip(pts, normals)]

# noisy-zero Dirichlet data (total pressure = 0) with independent noise σ per sensor
dirichlet_data(σ) = [MvNormal(σ .* randn(2), σ^2 * I(2)) for _ in pts]

# ---------------------------------------------------------------------------------------------
# TRAINING illumination: `n_incident` incident point sources equally spaced on a ring that
# encloses every obstacle. Their superposition lights the obstacles from every side at once,
# so the sources ARD keeps have to resolve the whole of every boundary, not just the lit side
# of one incidence.
# ---------------------------------------------------------------------------------------------
n_incident = 20
scene_centre = sum(centres) / length(centres)
incident_radius = 3.5
incident_positions = [
    scene_centre + incident_radius * SVector(cos(θ), sin(θ))
for θ in LinRange(0, 2π, n_incident + 1)[1:n_incident]]
@assert !any(inside_any, incident_positions) "an incident source landed inside an obstacle"
incident_sources = [PointSource([x], [one(ComplexF64)]) for x in incident_positions]

# the combined training trace (the vector sums the ring) and a 3%-of-peak noise level
train_on_boundary = trace_on_boundary(incident_sources)
σ_train = 0.03 * maximum(abs, train_on_boundary)

bd_train = BoundaryData(TractionType();
    boundary_points = pts,
    fields = dirichlet_data(σ_train),
    normals = normals,
    interior_points = centres
)

# ---------------------------------------------------------------------------------------------
# Candidate sources everywhere: one blind regular grid over the whole scene, keeping every
# node that falls inside an obstacle and is not touching the boundary. No shrunken boundary
# curves, no special sources at the spikes — the grid does not know where the spikes are.
# ---------------------------------------------------------------------------------------------
# The spacing and clearance are chosen so that the grid reaches INTO the narrow spike tips:
# a clearance much larger than the spike half-width would leave the tips without any nearby
# candidate, and no solver can keep a source that was never offered. The clearance is
# measured against a dense sampling of the true curves (not just the sensors), and kept
# near the sensor spacing, so no source can hide between two sensors and blow up the field
# there unseen.
spacing = 0.04
clearance = 0.02
curve_pts = vcat([[boundary_point(o, θ) for θ in LinRange(0, 2π, 1501)[1:1500]] for o in obstacles]...)
xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
grid = [SVector(x, y)
    for x in minimum(xs):spacing:maximum(xs), y in minimum(ys):spacing:maximum(ys)]
candidates = filter(vec(grid)) do p
    inside_any(p) && minimum(norm(p - q) for q in curve_pts) > clearance
end

@info "an overcomplete basis: more unknowns than data" candidates = length(candidates) coefficients = 2 * length(candidates) data = 2 * length(pts)

# ---------------------------------------------------------------------------------------------
# Solve #1 — SELECT the sources. Starting from an EMPTY basis, the solver repeatedly adds the
# single candidate coefficient that most increases the evidence, each at its exact optimal prior
# precision αᵢ, and stops when no candidate is worth adding — all under the combined ring
# illumination. The selected sources are the set we keep; this posterior is used only for its
# retained positions.
#
# Selection performs ONE action (add, re-estimate or delete a single coefficient) per iteration,
# so it has its own budget `max_select_iters`, separate from the `max_iters` whole-model sweeps:
# we leave it at its default, which allows four actions per candidate coefficient and stops the
# moment no action increases the evidence (about 1650 actions here). No `prior_variance` either:
# the selection phase derives every precision from the data, so an initial value is ignored.
#
# This is the expensive step of the example — ten seconds or so to sift 2800-odd candidates
# against 1200 data rows. It is also the step you only pay once: that is the whole point of
# selecting the sources for a ring of incidences and then reusing them.
# ---------------------------------------------------------------------------------------------
solver = VariationalBayesianSolver(
    ard_threshold = 1e6,
    max_iters = 300,
    elbo_tol = 1e-9
)

sim_train = Simulation(medium, bd_train;
    source_positions = candidates,
    particular_solution = incident_sources,   # a vector: the fields of the ring are summed
    solver = solver
)
t_prune = @elapsed vsol_train = solve(sim_train)

kept = vsol_train.fsol.positions
@info "the solver selected its sources from the candidate grid" kept = length(kept) of = length(candidates) seconds = round(t_prune, digits = 1) iterations = length(vsol_train.elbo_history) misfit_ratio = round(vsol_train.misfit_ratio, digits = 3)

# ---------------------------------------------------------------------------------------------
# Solve #2 — USE the sources. Keeping the selected positions FIXED (ARD pruning switched off),
# re-solve for ONE incident source. This single-source posterior is what every picture below
# draws from: the total field it predicts, and its uncertainty.
# ---------------------------------------------------------------------------------------------
gif_index = argmax([p[1] - p[2] for p in incident_positions])   # a lower-right incidence
x_source = incident_positions[gif_index]
incident = incident_sources[gif_index]

incident_on_boundary = trace_on_boundary(incident)
σ = 0.03 * maximum(abs, incident_on_boundary)      # 3% of the strongest boundary signal

bd = BoundaryData(TractionType();
    boundary_points = pts,
    fields = dirichlet_data(σ),
    normals = normals,
    interior_points = centres
)

solver_fixed = VariationalBayesianSolver(
    prior_variance = 1.0,
    ard_threshold = 1e6,
    ard_prune_flag = false,     # keep exactly the sources ARD selected above
    max_iters = 300,
    elbo_tol = 1e-9
)

sim = Simulation(medium, bd;
    source_positions = kept,
    particular_solution = incident,
    solver = solver_fixed
)
vsol = solve(sim)

# ---------------------------------------------------------------------------------------------
# Did it work? Check the Dirichlet condition at FRESH boundary points (halfway between the
# sensors — the sensors are fitted by construction), and check that the posterior's own 3σ
# error bars cover the truth (the total field should be 0 there).
# ---------------------------------------------------------------------------------------------
incident_scale = maximum(abs, incident_on_boundary)

# midpoints (in θ) between consecutive sensors of each obstacle
fresh = [
    (boundary_point(o, θm), boundary_normal(o, θm))
    for (o, θs) in zip(obstacles, θs_per)
    for θm in (θs .+ [θs[2:end]; θs[1] + 2π]) ./ 2
]
checks = map(fresh) do (p, n)
    f = field(TractionType(), vsol, p, n)[1]          # total field: should be ≈ 0
    stds = field_std(TractionType(), vsol, p, n)      # std of [Re; Im] of the scattered field
    covered = (abs(real(f)) <= 3 * stds[1] + 1e-9) + (abs(imag(f)) <= 3 * stds[2] + 1e-9)
    (abs(f), covered)
end
residuals = first.(checks)
coverage = sum(last.(checks)) / (2 * length(checks))

@info "Dirichlet condition at fresh boundary points" max_residual_rel = round(maximum(residuals) / incident_scale, digits = 4) mean_residual_rel = round(mean(residuals) / incident_scale, digits = 4) three_sigma_coverage = round(coverage, digits = 3)

# ---------------------------------------------------------------------------------------------
# Evaluate the total field and its uncertainty on a grid. For a complex field the posterior
# stores the covariance of [Re f; Im f]; we collapse it to one number per point,
#     s(x) = √(Var Re f + Var Im f),
# the phase-independent uncertainty of the scattered field at x.
# ---------------------------------------------------------------------------------------------
xlims, ylims = (-3.2, 4.0), (-3.3, 2.8)
nx, ny = 240, 210
xg = LinRange(xlims..., nx)
yg = LinRange(ylims..., ny)

vals = fill(NaN * im, ny, nx)     # total field (incident + scattered)
unc  = fill(NaN, ny, nx)          # s(x)
@time for (iy, y) in enumerate(yg), (ix, x) in enumerate(xg)
    p = SVector(x, y)
    inside_any(p) && continue
    vals[iy, ix] = field(TractionType(), vsol, p)[1]
    stds = field_std(TractionType(), vsol, p)
    unc[iy, ix] = sqrt(stds[1]^2 + stds[2]^2)
end

finite_unc = filter(isfinite, vec(unc))
@info "posterior field uncertainty s(x) over the grid" median = round(median(finite_unc), sigdigits = 3) q95 = round(quantile(finite_unc, 0.95), sigdigits = 3) max = round(maximum(finite_unc), sigdigits = 3)

# ---------------------------------------------------------------------------------------------
# Render the POSTERIOR MEAN of the total field: field value → diverging colour. The
# uncertainty is not mixed into this picture — it gets its own panel below, as a percentage
# of the mean field magnitude. The obstacles are filled dark grey; the retained sources are
# the small dots.
# ---------------------------------------------------------------------------------------------
mags = [abs(vals[iy, ix]) for iy in 1:ny for ix in 1:nx
    if isfinite(abs(vals[iy, ix])) && norm(SVector(xg[ix], yg[iy]) - x_source) > 0.3]
cmax = quantile(mags, 0.99)
cg = cgrad(:balance)
const OBSTACLE_COLOR = RGB{Float64}(0.22, 0.22, 0.25)

function render(t)
    img = Matrix{RGB{Float64}}(undef, ny, nx)
    for iy in 1:ny, ix in 1:nx
        v = vals[iy, ix]
        if !isfinite(abs(v))
            img[iy, ix] = OBSTACLE_COLOR
            continue
        end
        u = clamp((real(v * exp(-im * ω * t)) + cmax) / (2cmax), 0, 1)
        c = get(cg, u)
        img[iy, ix] = RGB{Float64}(c.r, c.g, c.b)
    end
    return img
end

θ_dense = LinRange(0, 2π, 400)
outlines = [([boundary_point(o, θ)[1] for θ in θ_dense], [boundary_point(o, θ)[2] for θ in θ_dense]) for o in obstacles]

function scene(img; title = "", show_sources = true)
    plt = plot(xg, yg, img;
        aspect_ratio = 1, xlims = xlims, ylims = ylims, yflip = false,
        axis = false, grid = false, title = title, titlefontsize = 11)
    for (ox, oy) in outlines
        plot!(plt, ox, oy; lc = :black, lw = 1.2, label = "")
    end
    if show_sources
        scatter!(plt, [s[1] for s in kept], [s[2] for s in kept];
            mc = :crimson, ms = 1.6, msw = 0, label = "")
    end
    scatter!(plt, [x_source[1]], [x_source[2]]; mc = :lime, ms = 6, msw = 1, label = "")
    return plt
end

# --- main advert figure: the posterior mean of the total field, next to the uncertainty
#     map expressed relative to the mean field magnitude ---
plt_field = scene(render(0.0); title = "total field (posterior mean)")

mean_field_mag = mean(mags)
unc_rel = 100 .* unc ./ mean_field_mag
plt_unc = heatmap(xg, yg, unc_rel;
    c = :viridis, aspect_ratio = 1, xlims = xlims, ylims = ylims,
    axis = false, grid = false, colorbar = true,
    title = "uncertainty % of the |mean|. Boundary error 3%", titlefontsize = 11)
for (ox, oy) in outlines
    plot!(plt_unc, Shape(ox, oy); fillcolor = OBSTACLE_COLOR, lc = :black, lw = 1.2, label = "")
end

plt_main = plot(plt_field, plt_unc; layout = (1, 2), size = (1250, 520))
savefig(plt_main, joinpath(FIGDIR, "spikey_field_uncertainty.png"))

# --- the posterior mean of the total field over one period ---
anim = @animate for t in LinRange(0, 2π / ω, 32)[1:31]
    scene(render(t); title = "total field (posterior mean)")
end
gif(anim, joinpath(FIGDIR, "spikey_scattering.gif"), fps = 10)

# --- what ARD did: the blind candidate grid next to the retained sources,
#     sized by the posterior magnitude of their coefficients ---
plt_cand = plot(; aspect_ratio = 1, xlims = xlims, ylims = ylims, axis = false, grid = false,
    title = "$(length(candidates)) candidates · $(n_incident) training incidences", titlefontsize = 11)
for (ox, oy) in outlines
    plot!(plt_cand, ox, oy; lc = :black, lw = 1.2, label = "")
end
scatter!(plt_cand, [s[1] for s in candidates], [s[2] for s in candidates];
    mc = :gray, ms = 1.6, msw = 0, label = "")
scatter!(plt_cand, [s[1] for s in incident_positions], [s[2] for s in incident_positions];
    mc = :orange, ms = 4, msw = 0, label = "")               # the training ring, all around
scatter!(plt_cand, [x_source[1]], [x_source[2]]; mc = :lime, ms = 6, msw = 1, label = "")  # the one animated

amp = abs.(vsol.fsol.coefficients)
plt_kept = plot(; aspect_ratio = 1, xlims = xlims, ylims = ylims, axis = false, grid = false,
    title = "$(length(kept)) selected ($(round(t_prune, digits = 1))s), size = |coefficient|", titlefontsize = 11)
for (ox, oy) in outlines
    plot!(plt_kept, ox, oy; lc = :black, lw = 1.2, label = "")
end
scatter!(plt_kept, [s[1] for s in kept], [s[2] for s in kept];
    mc = :crimson, ms = 1.5 .+ 5 .* amp ./ maximum(amp), msw = 0, label = "")

plt_ard = plot(plt_cand, plt_kept; layout = (1, 2), size = (1250, 520))
savefig(plt_ard, joinpath(FIGDIR, "spikey_ard_pruning.png"))

@info "wrote figures to" FIGDIR
