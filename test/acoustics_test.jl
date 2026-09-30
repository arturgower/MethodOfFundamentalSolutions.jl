using LinearAlgebra
using Test
using SpecialFunctions
using MethodOfFundamentalSolutions

@testset "acoustic scattering" begin

N_bd = 100;
r = 1.0;
x0,y0 = 0.0,0.0;
ω=2.0;

ϕ(x, y) =-(im/4)*hankelh1(0, (ω/medium.c) * sqrt((x-x0)^2 + (y-y0)^2))

N_sources=N_bd;
λ = 1e-10;
tolerance = 1e-10;
res = 51;

medium = Acoustic(2; ω = ω,  ρ = 1.0, c = 1.0)

θs = LinRange(0,2pi,N_bd+1)[1:N_bd]
bd_points = [[r*cos(θ), r*sin(θ)] for θ in θs]
normals = [[cos(θ), sin(θ)] for θ in θs]
bd_fields = [[-ϕ(r*cos(θ), r*sin(θ))] for θ in θs]
interior_points = [[0.0, 0.0]]

bd = BoundaryData(TractionType(); 
    boundary_points = bd_points, 
    fields = bd_fields, 
    normals = normals,
    interior_points = interior_points
)

source_pos = source_positions(bd; relative_source_distance = 1.0)

# Solve
solver = TikhonovSolver(λ=λ, tolerance = tolerance)
sim = Simulation(medium,bd, solver=solver, source_positions = source_pos)
fsol = solve(sim)

predict_fields = [field(TractionType(), fsol, bd_points[i], normals[i]) for i in eachindex(bd_points)]
fields = [-ϕ(r*cos(θ),r*sin(θ)) for θ in θs]
f=vcat(bd.fields...)

errors = [abs(fields[i] - predict_fields[i][1]) for i in eachindex(fields)]
@test maximum(errors) < 1e-10

# Splitting one BoundaryData into two disjoint parts must give exactly the same solution:
# the tuple of BoundaryData shares the same sources (columns) and simply stacks its rows,
# so with the same source positions the system matrix and forcing are unchanged.
@testset "split boundary data" begin
    half = N_bd ÷ 2
    idx1, idx2 = 1:half, half+1:N_bd

    bd1 = BoundaryData(TractionType();
        boundary_points = bd_points[idx1], fields = bd_fields[idx1],
        normals = normals[idx1], interior_points = interior_points)
    bd2 = BoundaryData(TractionType();
        boundary_points = bd_points[idx2], fields = bd_fields[idx2],
        normals = normals[idx2], interior_points = interior_points)

    sim_split = Simulation(medium, (bd1, bd2); solver = solver, source_positions = source_pos)
    fsol_split = solve(sim_split)

    # the stacked block system is identical to the single one, so are its coefficients
    @test system_matrix(sim_split) ≈ system_matrix(sim)
    @test fsol_split.coefficients ≈ fsol.coefficients

    predict_split = [field(TractionType(), fsol_split, bd_points[i], normals[i]) for i in eachindex(bd_points)]
    errors_split = [abs(fields[i] - predict_split[i][1]) for i in eachindex(fields)]
    @test maximum(errors_split) < 1e-10
end

end

# ==============================================================================
# Transmission (penetrable) scatterer: on the interface both the field φ (pressure,
# TractionType) and the normal displacement (1/ρ) ∂φ/∂n (DisplacementType) are
# continuous. The incident field is one point source outside the scatterer; the
# scattered field is represented by sources inside the scatterer (exterior medium),
# and the transmitted field by sources outside it (interior medium).
# ==============================================================================
@testset "acoustic transmission" begin

using StaticArrays

@testset "displacement greens is (1/ρ) ∂φ/∂n" begin
    medium = Acoustic(2; ω = 1.3, ρ = 1.7, c = 0.9)
    x = SVector(0.7, -0.4)
    nvec = SVector(cos(0.4), sin(0.4))

    h = 1e-6
    fd = (greens(TractionType(), medium, x + h * nvec)[1] -
          greens(TractionType(), medium, x - h * nvec)[1]) / (2h * medium.ρ)

    @test greens(DisplacementType(), medium, x, nvec)[1] ≈ fd rtol = 1e-8
end

@testset "penetrable circle vs the analytic partial-wave solution" begin
    ω = 2.0
    medium0 = Acoustic(2; ω = ω, ρ = 1.0, c = 1.0)   # exterior
    medium1 = Acoustic(2; ω = ω, ρ = 1.7, c = 0.6)   # interior (the scatterer)
    k0 = ω / medium0.c
    k1 = ω / medium1.c

    a = 1.0                       # radius of the scatterer
    xs = SVector(2.5, 0.5)        # the incident point source, outside the scatterer
    rs, θs = norm(xs), atan(xs[2], xs[1])

    # incident field of the point source and its normal displacement
    φ_inc(p) = greens(TractionType(), medium0, SVector{2,Float64}(p) - xs)[1]
    v_inc(p, n) = greens(DisplacementType(), medium0, SVector{2,Float64}(p) - xs, SVector{2,Float64}(n))[1]

    # --- analytic solution: partial waves with Graf's addition theorem for the source ---
    Jp(n, z) = (besselj(n - 1, z) - besselj(n + 1, z)) / 2
    Hp(n, z) = (hankelh1(n - 1, z) - hankelh1(n + 1, z)) / 2

    ns = -30:30
    Cn = [-(im / 4) * hankelh1(n, k0 * rs) * exp(-im * n * θs) for n in ns]

    AB = map(eachindex(ns)) do i
        n = ns[i]
        M2 = [hankelh1(n, k0 * a)          -besselj(n, k1 * a);
              (k0 / medium0.ρ) * Hp(n, k0 * a)   -(k1 / medium1.ρ) * Jp(n, k1 * a)]
        M2 \ (-Cn[i] .* [besselj(n, k0 * a), (k0 / medium0.ρ) * Jp(n, k0 * a)])
    end
    An = [ab[1] for ab in AB]
    Bn = [ab[2] for ab in AB]

    polar(x) = (norm(x), atan(x[2], x[1]))
    φ_sc_analytic(x) = ((r, θ) = polar(x); sum(An[i] * hankelh1(ns[i], k0 * r) * exp(im * ns[i] * θ) for i in eachindex(ns)))
    φ_int_analytic(x) = ((r, θ) = polar(x); sum(Bn[i] * besselj(ns[i], k1 * r) * exp(im * ns[i] * θ) for i in eachindex(ns)))

    # guard: the Graf expansion of the incident field must match the point source on r = a
    p_check = [a * cos(0.3), a * sin(0.3)]
    φ_inc_series = sum(Cn[i] * besselj(ns[i], k0 * a) * exp(im * ns[i] * 0.3) for i in eachindex(ns))
    @test φ_inc_series ≈ φ_inc(p_check) rtol = 1e-10

    # --- MFS: build the transmission simulation and solve ---
    N_bd = 64
    θ_bd = LinRange(0, 2pi, N_bd + 1)[1:N_bd]
    points = [[a * cos(θ), a * sin(θ)] for θ in θ_bd]
    normals = [[cos(θ), sin(θ)] for θ in θ_bd]
    interface = BoundaryShape(boundary_points = points, normals = normals, interior_points = [[0.0, 0.0]])

    n_src = 32
    θ_src = LinRange(0, 2pi, n_src + 1)[1:n_src]
    src_in = [SVector(0.6 * cos(θ), 0.6 * sin(θ)) for θ in θ_src]    # scattered field
    src_out = [SVector(1.6 * cos(θ), 1.6 * sin(θ)) for θ in θ_src]   # transmitted field

    # the two continuity conditions (φ and (1/ρ)∂φ/∂n) share the same interface; the prescribed
    # boundary fields are zero (pure continuity) — the incident enters via the particular solution
    zero_fields = [[0.0im] for _ in θ_bd]
    conditions = (BoundaryData(TractionType(); boundary_shape = interface, fields = zero_fields),
                  BoundaryData(DisplacementType(); boundary_shape = interface, fields = zero_fields))

    # the incident point source at xs enters as the exterior region's particular solution
    incident = PointSource([xs], [one(ComplexF64)])

    solver = TikhonovSolver(λ = 0.0)   # plain least squares, as in the analytic comparison
    sim_ext = Simulation(medium0, conditions;
        source_positions = src_in, particular_solution = incident, solver = solver)   # scattered field
    sim_int = Simulation(medium1, conditions;
        source_positions = src_out, solver = solver)                                  # transmitted field

    fsol_sc, fsol_int = solve(TransmissionSimulation(sim_ext, sim_int))
    @test fsol_sc.relative_boundary_error < 1e-4

    # the PointSource particular reproduces the incident point-source trace
    @test field(TractionType(), medium0, incident, points[1], normals[1])[1] ≈ φ_inc(points[1]) rtol = 1e-12

    # fsol_sc carries the incident as its particular solution, so field(fsol_sc, ·) is the
    # exterior *total* field (incident + scattered); fsol_int is the transmitted field

    # --- the MFS fields must match the analytic solution away from the boundary ---
    θ_test = LinRange(0, 2pi, 41)[1:40]

    p_ext = [[1.4 * cos(θ), 1.4 * sin(θ)] for θ in θ_test]
    err_sc = maximum(abs(field(TractionType(), fsol_sc, p)[1] - (φ_inc(p) + φ_sc_analytic(p))) for p in p_ext) /
        maximum(abs.(φ_sc_analytic.(p_ext)))
    @test err_sc < 1e-6

    p_int = [[0.5 * cos(θ), 0.5 * sin(θ)] for θ in θ_test]
    err_int = maximum(abs(field(TractionType(), fsol_int, p)[1] - φ_int_analytic(p)) for p in p_int) /
        maximum(abs.(φ_int_analytic.(p_int)))
    @test err_int < 1e-6

    # --- and satisfy both transmission conditions at fresh points on the interface ---
    θ_fresh = θ_test .+ pi / 40   # between the collocation points
    field_scale = maximum(abs(φ_inc(p)) for p in p_int)
    for θ in θ_fresh
        p = [a * cos(θ), a * sin(θ)]
        nvec = [cos(θ), sin(θ)]

        # field(fsol_sc, ·) already includes the incident, so continuity is total_ext = total_int
        jump_φ = field(TractionType(), fsol_sc, p)[1] - field(TractionType(), fsol_int, p)[1]
        jump_v = field(DisplacementType(), fsol_sc, p, nvec)[1] -
            field(DisplacementType(), fsol_int, p, nvec)[1]

        @test abs(jump_φ) < 1e-5 * field_scale
        @test abs(jump_v) < 1e-4 * field_scale * abs(k0)
    end
end

end