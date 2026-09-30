# # Point-source scattering from a tear-drop scatterer, and why source placement matters
#
# A sound-soft tear-drop scatterer (total pressure = 0 on its boundary) is illuminated by an
# incident point source. The tear-drop has a *cusp* — a sharp point where the boundary is
# barely resolved — and this is exactly where the Method of Fundamental Solutions struggles
# unless one of the interior sources is placed right next to it.
#
# The scattered field is represented by MFS sources placed INSIDE the obstacle. We solve twice,
# changing only whether one source sits next to the cusp:
#   A. interior sources + one source next to the cusp;
#   B. the same interior sources, cusp left uncovered.
#
# The total field (incident + scattered) is animated over one period, and the boundary residual
# |total field| is compared between A and B — it should blow up at the cusp in case B.
#
# Run from the package root with:  julia --project docs/examples/acoustic/teardrop_scattering.jl

using MethodOfFundamentalSolutions
using LinearAlgebra, Statistics, StaticArrays
using Plots

const FIGDIR = joinpath(@__DIR__, "figures")

# ---------------------------------------------------------------------------------------------
# Tear-drop geometry.  The curve (cos θ, sin θ · sin²(θ/2)) has a cusp at θ = 0, i.e. at (1, 0),
# and a smooth round bulb at θ = π.  Outward normals come from the tangent, oriented away from
# the centroid.
# ---------------------------------------------------------------------------------------------
teardrop(θ; scale = 1.0) = scale .* SVector(cos(θ), sin(θ) * sin(θ / 2)^2)

function teardrop_normal(θ, centroid; scale = 1.0, h = 1e-5)
    t = teardrop(θ + h; scale) - teardrop(θ - h; scale)   # tangent
    n = SVector(t[2], -t[1]) / norm(SVector(t[2], -t[1]))  # rotate 90°
    return dot(n, teardrop(θ; scale) - centroid) < 0 ? -n : n
end

scale    = 1.0
N_bd     = 220
θs       = collect(LinRange(0, 2π, N_bd + 1)[1:N_bd])
pts      = [teardrop(θ; scale) for θ in θs]
centroid = sum(pts) / N_bd
normals  = [teardrop_normal(θ, centroid; scale) for θ in θs]
cusp     = teardrop(0.0; scale)                 # the sharp point, at (scale, 0)

interface = BoundaryShape(boundary_points = pts, normals = normals, interior_points = [centroid])

# ---------------------------------------------------------------------------------------------
# Physics: a sound-soft tear-drop lit by an incident point source to the right, so the wave
# strikes the cusp head-on.  Sound-soft means the total pressure vanishes on the boundary, so
# the scattered field must equal minus the incident field there.  We encode this by giving the
# boundary zero prescribed field and putting the incident wave in the particular solution: the
# solver then matches `scattered = -incident` on the boundary, and the returned solution carries
# the incident field, so `field(fsol, x)` is the total field (incident + scattered).
# ---------------------------------------------------------------------------------------------
ω        = 2π * 1.0
medium   = Acoustic(2; ω = ω, ρ = 1.0, c = 1.0)
x_source = SVector(3.0, 0.6)
incident = PointSource([x_source], [one(ComplexF64)])

bd = BoundaryData(TractionType(); boundary_shape = interface, fields = [[0.0im] for _ in θs])

# ---------------------------------------------------------------------------------------------
# Interior MFS sources.  Shrinking the boundary toward the centroid gives a smooth interior
# source curve — but shrinking pulls the cusp point back to ≈(0.7, 0), so this base leaves the
# cusp tip uncovered.  The single experimental variable is one extra source placed just inside
# the cusp tip.
# ---------------------------------------------------------------------------------------------
sub(N)  = round.(Int, LinRange(1, N_bd, N + 1))[1:N]
shrink  = 0.75
base_in = [centroid + shrink * (pts[i] - centroid) for i in sub(120)]   # smooth interior curve

cusp_source   = SVector(0.96, 0.0)                # just inside the cusp tip
src_in_cusp   = vcat(base_in, [cusp_source])      # A: cusp covered
src_in_nocusp = base_in                           # B: cusp uncovered

@info "closest interior source to the cusp" with_cusp = minimum(norm(s - cusp) for s in src_in_cusp) without = minimum(norm(s - cusp) for s in src_in_nocusp)

# ---------------------------------------------------------------------------------------------
# Solve the scattering problem for a given set of interior sources.
# ---------------------------------------------------------------------------------------------
# plain least squares (λ = 0): without a source at the cusp the fit cannot represent the field
# there and develops spurious oscillations — the whole point of this example
solver = TikhonovSolver(λ = 0.0)

solve_config(src_in) = solve(Simulation(medium, bd;
    source_positions = src_in, particular_solution = incident, solver = solver))

fsol_A = solve_config(src_in_cusp)     # A: cusp covered
fsol_B = solve_config(src_in_nocusp)   # B: cusp uncovered

# ---------------------------------------------------------------------------------------------
# Boundary residual: for a sound-soft scatterer the total field must vanish on the boundary, so
# |field(fsol, p)| at fresh points between the collocation points measures the error.
# ---------------------------------------------------------------------------------------------
θ_fresh = θs .+ (θs[2] - θs[1]) / 2
boundary_residual(fsol) = [abs(field(TractionType(), fsol, teardrop(θ; scale))[1]) for θ in θ_fresh]

res_A = boundary_residual(fsol_A)
res_B = boundary_residual(fsol_B)
near_cusp = @. (θ_fresh < deg2rad(60)) | (θ_fresh > 2π - deg2rad(60))   # cusp error spreads to ≈±60°
@info "boundary residual (should be ≈0)" max_A = maximum(res_A) max_B = maximum(res_B) cusp_A = maximum(res_A[near_cusp]) cusp_B = maximum(res_B[near_cusp])

θ_deg  = θ_fresh .* (180 / π)
presid = plot(θ_deg, res_B; yscale = :log10, label = "B: cusp uncovered", lw = 2, c = :crimson,
    xlabel = "boundary angle θ (degrees)  —  cusp at 0°/360°", ylabel = "|total field| on boundary",
    title = "Boundary residual vs. interior-source placement", legend = :top)
plot!(presid, θ_deg, res_A; label = "A: source at the cusp", lw = 2, c = :seagreen)
vline!(presid, [0, 360]; c = :gray, ls = :dash, label = "cusp")
savefig(presid, joinpath(FIGDIR, "teardrop_boundary_residual.png"))

# ---------------------------------------------------------------------------------------------
# Total field on a grid (incident + scattered); the obstacle interior is masked out.
# ---------------------------------------------------------------------------------------------
xlims, ylims = (-2.1, 3.5), (-2.1, 2.1)
nx, ny = 250, 190
xg = LinRange(xlims..., nx)
yg = LinRange(ylims..., ny)
grid = [SVector(x, y) for y in yg for x in xg]   # x varies fastest

grid_field(fsol) = [in(x, interface) ? NaN * im : field(TractionType(), fsol, x)[1] for x in grid]

vals_A = grid_field(fsol_A)
vals_B = grid_field(fsol_B)

# colour scale from the clean (cusp-covered) solution, so B's spurious oscillations saturate
mags = [abs(v) for (v, x) in zip(vals_A, grid) if isfinite(abs(v)) && norm(x - x_source) > 0.25]
cmax = quantile(mags, 0.99)

outline_x = [p[1] for p in pts]; push!(outline_x, pts[1][1])
outline_y = [p[2] for p in pts]; push!(outline_y, pts[1][2])

function make_gif(vals, src_in, filename, title)
    ts = LinRange(0, 2π / ω, 32)
    anim = @animate for t in ts
        Z = permutedims(reshape(real.(vals .* exp(-im * ω * t)), nx, ny))
        heatmap(xg, yg, Z; clims = (-cmax, cmax), c = :balance, aspect_ratio = 1,
            xlims = xlims, ylims = ylims, colorbar = false, axis = false, title = title)
        plot!(outline_x, outline_y; lc = :black, lw = 1.5, label = "")
        scatter!([s[1] for s in src_in], [s[2] for s in src_in];
            mc = :black, ms = 2.2, msw = 0, label = "")
        scatter!([x_source[1]], [x_source[2]]; mc = :lime, ms = 5, label = "")
    end
    path = joinpath(FIGDIR, filename)
    gif(anim, path, fps = 10)
    return path
end

make_gif(vals_A, src_in_cusp,   "teardrop_with_cusp_source.gif",    "A: interior source at the cusp")
make_gif(vals_B, src_in_nocusp, "teardrop_without_cusp_source.gif", "B: cusp left uncovered")

# ---------------------------------------------------------------------------------------------
# Static zoom on the cusp, at a fixed phase and a sensitive colour scale, so the spurious
# oscillation that appears in B (but not A) right at the cusp is visible side by side.
# ---------------------------------------------------------------------------------------------
function cusp_panel(fsol, src_in, title)
    zx, zy = LinRange(0.4, 1.55, 200), LinRange(-0.62, 0.62, 200)
    Z = [in(SVector(x, y), interface) ? NaN : real(field(TractionType(), fsol, SVector(x, y))[1])
         for y in zy, x in zx]
    p = heatmap(zx, zy, Z; clims = (-0.15, 0.15), c = :balance, aspect_ratio = 1,
        xlims = (0.4, 1.55), ylims = (-0.62, 0.62), colorbar = false, title = title)
    plot!(p, outline_x, outline_y; lc = :black, lw = 2, label = "")
    scatter!(p, [s[1] for s in src_in], [s[2] for s in src_in]; mc = :black, ms = 3, label = "")
    return p
end

cusp_fig = plot(cusp_panel(fsol_B, src_in_nocusp, "B: cusp uncovered — spurious field"),
                cusp_panel(fsol_A, src_in_cusp,   "A: source at the cusp — clean"),
                layout = (1, 2), size = (900, 460))
savefig(cusp_fig, joinpath(FIGDIR, "teardrop_cusp_zoom.png"))

@info "wrote figures to" FIGDIR
