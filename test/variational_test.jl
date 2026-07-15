# The evidence lower bound must never decrease, except across iterations where the model
# itself changed (a source was pruned, or the boundary was re-centered).
function elbo_is_monotone(vsol; tol = 1e-6)
    F = vsol.elbo_history
    return all(eachindex(F)) do i
        i == 1 || (i in vsol.baseline_resets) || F[i] >= F[i - 1] - tol * (1 + abs(F[i - 1]))
    end
end

# All the tests below keep the geometry x fixed: update_geometry_flag = false (the default),
# so the algorithm reduces to exact evidence maximization over the prior {αᵢ} and basis χ.
@testset "Variational evidence maximization" begin


# ==============================================================================
# 0. Synthetic data for the Laplace equation, where the exact source positions are
#    known. Each source i has a Gaussian amplitude aᵢ ~ N(μᵢ, vᵢ); the field is linear
#    in the amplitudes, so the boundary field is Gaussian with
#        mean u(x) = Σᵢ μᵢ G(x - χᵢ)   and   variance Var u(x) = Σᵢ vᵢ G(x - χᵢ)².
#    The variance of the measurements is carried by the boundary data alone, as a
#    vector of `MvNormal`s, one per boundary point.
# ==============================================================================
@testset "Laplace equation: circle" begin
    medium = LaplaceMedium{2, Float64}()
    Random.seed!(1234)
    FT = DirichletType()

    # The sources: positions, mean amplitudes, and a variance for each amplitude
    r_source = 1.1
    N_sources = 10
    source_θs = LinRange(0, 2pi, N_sources + 1)[1:N_sources] + rand(N_sources) * 0.1
    source_pos = [[r_source * cos(θ), r_source * sin(θ)] for θ in source_θs]
    
    # want all amps to same strength, as then they can be learned.
    signs = [ θ > pi ? 1.0 : -1.0  for θ in source_θs]
    source_amps = (-1.0) .^ signs
    source_vars = (0.05 .* source_amps) .^ 2   # a 5% standard deviation per source

    # The boundary points
    r = 1.0
    N_points = 40
    points_θs = LinRange(0, 2pi, N_points + 1)[1:N_points]
    points = [[r * cos(θ), r * sin(θ)] for θ in points_θs]

    # 1) evaluate the field of the combined sources: its mean and its variance
    field_mean(x) = sum(source_amps[i] * greens(FT, medium, x - source_pos[i]) for i in eachindex(source_pos))
    field_var(x) = sum(source_vars[i] * greens(FT, medium, x - source_pos[i])^2 for i in eachindex(source_pos))

    bd = BoundaryData(FT;
        boundary_points = points,
        fields = [MvNormal([field_mean(x)], field_var(x) * I(1)) for x in points]
    )

    M = system_matrix([SVector{2, Float64}(p) for p in source_pos], medium, bd)
    Σ_induced = Symmetric(M * Diagonal(source_vars) * M' + 1e-9 * I)
    bd_full = BoundaryData(FT;
        boundary_points = points,
        fields = MvNormal(field_mean.(points), Σ_induced)
    )

    @test cov(bd_full.fields) ≈ Σ_induced

    # using Plots
    # plot(bd)

    # 2) With the true source positions, and the true amplitude variances passed as the
    #    priors, the solver must use exactly the specified variances of the sources and
    #    return exactly the corresponding Gaussian posterior.
    @testset "correct positions and priors: the specified variances are used exactly" begin
        priors = [MvNormal([0.0], source_vars[i] * I(1)) for i in eachindex(source_vars)]
        solver = VariationalBayesianSolver(
            priors = priors,
            learn_prior_flag = false,
            ard_prune_flag = false,
            max_iters = 1
        )
        sim = Simulation(medium, bd_full; solver = solver, source_positions = source_pos)
        vsol = solve(sim)

        # using Plots
        # plot(vsol,bd)

        # the prior variances of the sources are exactly those specified
        @test 1 ./ vsol.prior_precisions ≈ source_vars rtol = 1e-12

        # and the posterior is exactly the Gaussian posterior with those variances
        M = system_matrix([SVector{2, Float64}(p) for p in source_pos], medium, bd_full)
        invΣ = inv(Σ_induced)

        Σ_exact = inv(Symmetric(M' * (invΣ * M) + Diagonal(1 ./ source_vars)))
        μ_exact = Σ_exact * (M' * (invΣ * field_mean.(points)))
        @test vsol.fsol.coefficients ≈ μ_exact rtol = 1e-8
        @test vsol.fsol.coefficients_covariance ≈ Σ_exact rtol = 1e-8
    end

    # 2b) The field covariance induced by independent source amplitudes aᵢ ~ N(μᵢ, vᵢ) is the
    #     FULL matrix Σ = M diag(v) Mᵀ (rank-deficient, so a small jitter makes it a proper
    #     covariance). Fed as the measurement noise, with a weak prior so the data determine
    #     the posterior, the coefficient posterior covariance recovers diag(vᵢ): the true
    #     variances of the sources. (With only the diagonal of Σ, as in test 2, it is close
    #     but not exact; the full covariance makes it exact.)
    @testset "full induced covariance: posterior variance matches the source variances" begin
        M = system_matrix([SVector{2, Float64}(p) for p in source_pos], medium, bd_full)
        Σ_induced = Symmetric(M * Diagonal(source_vars) * M' + 1e-9 * I)
        bd_full = BoundaryData(FT;
            boundary_points = points,
            fields = MvNormal(field_mean.(points), Σ_induced)
        )

        solver = VariationalBayesianSolver(
            prior_variance = 1e8,        # weak prior: the posterior variance is set by the data
            learn_prior_flag = false,
            ard_prune_flag = false,
            max_iters = 1
        )
        sim = Simulation(medium, bd_full; solver = solver, source_positions = source_pos)
        vsol = solve(sim)
        
        @test diag(vsol.fsol.coefficients_covariance) ≈ source_vars rtol = 1e-5

        # posterior variance is diagonal
        @test norm(vsol.fsol.coefficients_covariance - diagm(diag(vsol.fsol.coefficients_covariance))) / norm(vsol.fsol.coefficients_covariance) < 1e-5
        
        # if we learn the prior from the data, the posterior variance approximately recover the true source variances
        solver = VariationalBayesianSolver(
            prior_variance = 1e0,        # weak prior: the posterior variance is set by the data
            learn_prior_flag = true,
            ard_prune_flag = false,
            max_iters = 100
        )
        sim = Simulation(medium, bd_full; solver = solver, source_positions = source_pos)
        vsol = solve(sim)
        
        @test diag(vsol.fsol.coefficients_covariance) ≈ source_vars rtol = 1e-2
        @test norm(vsol.fsol.coefficients_covariance - diagm(diag(vsol.fsol.coefficients_covariance))) / norm(vsol.fsol.coefficients_covariance) < 1e-5

    end

    # 3) Misspecified priors: evidence maximization re-learns the prior, so the
    #    posterior mean still recovers the amplitudes exactly (the data are the exact
    #    mean field of the sources).
    @testset "misspecified priors still recover the exact solution" begin
        priors_wrong = [MvNormal([0.0], 1e3 * source_vars[i] * I(1)) for i in eachindex(source_vars)]
        solver = VariationalBayesianSolver(
            priors = priors_wrong,
            ard_prune_flag = false,
            max_iters = 200,
            elbo_tol = 1e-10
        )
        sim = Simulation(medium, bd_full; solver = solver, source_positions = source_pos)
        vsol = solve(sim)

        # recovery well below the 5% standard deviation of the amplitudes
        @test vsol.fsol.coefficients ≈ source_amps rtol = 2e-2
        pred = [field(FT, vsol, x)[1] for x in points]
        @test pred ≈ field_mean.(points) rtol = 2e-2
        @test elbo_is_monotone(vsol)
    end

    # 4) A large number of candidate sources: ARD prunes the superfluous ones and the
    #    source position update moves the survivors, recovering the true field.
    #    (docs/examples/variational/laplace_source_learning.jl draws the frames of this
    #    run: the sources moving, being pruned, and the field they generate.)
    @testset "ARD with source position updates learns the sources" begin
        N_dense = 60
        dense_θs = LinRange(0, 2pi, N_dense + 1)[1:N_dense]
        sources0 = [[1.35 * cos(θ), 1.35 * sin(θ)] for θ in dense_θs]

        solver = VariationalBayesianSolver(
            prior_variance = 1.0,
            optimise_source_positions_flag = true,
            ard_threshold = 1e6,
            max_iters = 150,
            elbo_tol = 1e-9
        )
        sim = Simulation(medium, bd; solver = solver, source_positions = sources0)
        vsol = solve(sim)

        # ARD pruned most of the candidates, down to near the number of true sources
        @test length(vsol.fsol.positions) < length(sources0) / 3
        @test elbo_is_monotone(vsol)

        # the retained sources reproduce the field, on fresh boundary points and in the
        # interior, to within the 5% noise level the boundary data specified
        θ_test = LinRange(0, 2pi, 101)[1:100]
        for r_test in (r, 0.9r)
            p_test = [[r_test * cos(θ), r_test * sin(θ)] for θ in θ_test]
            scale = maximum(abs.(field_mean.(p_test)))
            errs = [abs(field(FT, vsol, p)[1] - field_mean(p)) for p in p_test] ./ scale
            @test mean(errs) < 0.05
        end
    end

end

# ==============================================================================
# 1. The analytic gradient of the expected misfit over the source positions χ
#    must match finite differences.
# ==============================================================================
@testset "source position gradient consistency" begin
    Random.seed!(123)
    medium = Elastostatic(2; ρ = 1.0, cp = 2.0, cs = 1.0)

    n_bd = 8
    θs = LinRange(0, 2pi, n_bd + 1)[1:n_bd]
    points = [[1.3cos(θ), 1.3sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]
    bd = BoundaryData(TractionType();
        boundary_points = points,
        fields = [randn(2) for _ in θs],
        normals = normals,
        interior_points = [[0.0, 0.0]]
    )

    n_src = 5
    θsrc = LinRange(0, 2pi, n_src + 1)[1:n_src]
    chi = vcat([[2.1cos(θ), 2.1sin(θ)] for θ in θsrc]...)

    FD = 2
    K = FD * n_src
    N = FD * n_bd
    μ = randn(K)
    A = randn(K, K)
    Σpost = A' * A + I
    w = rand(N) .+ 0.5
    g = randn(N)

    R(c) = MFS._chi_misfit(c, medium, bd, false, μ, Σpost, w, g, 2)

    G = zeros(length(chi))
    MFS._chi_misfit_gradient!(G, chi, medium, bd, μ, Σpost, w, g, 2, FD)

    h = 1e-6
    G_fd = map(eachindex(chi)) do i
        cp = copy(chi); cm = copy(chi)
        cp[i] += h; cm[i] -= h
        (R(cp) - R(cm)) / (2h)
    end

    @test isapprox(G, G_fd; rtol = 1e-5)
end

# ==============================================================================
# 2. With the learning switched off, the variational solver must reproduce the
#    exact fixed-boundary Gaussian posterior of the existing Bayesian solver.
# ==============================================================================
@testset "reduces to the exact fixed-boundary posterior" begin
    Random.seed!(11)
    medium = Elastostatic(2; ρ = 1.0, cp = 3.0, cs = 2.0)

    n = 10; L = 0.5; H = 0.5
    x = [
        LinRange(-L, L, n + 2)[2:end - 1]; zeros(n) .+ L;
        LinRange(L, -L, n + 2)[2:end - 1]; zeros(n) .- L
    ]
    y = [
        zeros(n); LinRange(0, H, n + 2)[2:end - 1];
        zeros(n) .+ H; LinRange(H, 0, n + 2)[2:end - 1]
    ]
    points = [[x[i], y[i]] for i in eachindex(x)]
    normals = [
        [[0.0, -1.0] for _ = 1:n]; [[1.0, 0.0] for _ = 1:n];
        [[0.0, 1.0] for _ = 1:n]; [[-1.0, 0.0] for _ = 1:n]
    ]

    weight = 2L * H * medium.ρ * 9.81
    traction = [
        [[0.0, weight / (2L)] for _ = 1:n]; [[0.0, 0.0] for _ = 1:n];
        [[0.0, 0.0] for _ = 1:n]; [[0.0, 0.0] for _ = 1:n]
    ]
    flat_traction = vcat(traction...)

    σ_noise = 0.01 * weight / (2L)
    g_noisy = flat_traction .+ σ_noise .* randn(length(flat_traction))
    field_distribution = MvNormal(g_noisy, σ_noise^2 * I(length(g_noisy)))

    bd = BoundaryData(TractionType();
        boundary_points = points,
        fields = field_distribution,
        normals = normals
    )

    ns = 10
    sources = Vector{Float64}[]
    for x_val in LinRange(-2.0, 2.0, ns + 1)[1:end - 1]; push!(sources, [x_val, -1.75]); end
    for y_val in LinRange(-1.75, 2.25, ns + 1)[1:end - 1]; push!(sources, [2.0, y_val]); end
    for x_val in LinRange(2.0, -2.0, ns + 1)[1:end - 1]; push!(sources, [x_val, 2.25]); end
    for y_val in LinRange(2.25, -1.75, ns + 1)[1:end - 1]; push!(sources, [-2.0, y_val]); end

    prior = MvNormal(zeros(2 * length(sources)), 100.0^2 * I(2 * length(sources)))

    # exact posterior from the existing machinery; the boundary points are deterministic,
    # so its geometric covariance Cx vanishes and it is the exact Gaussian posterior
    bsim = Simulation(medium, bd;
        solver = BayesianSolver(prior),
        source_positions = sources,
        particular_solution = ParticularGravity(height = H)
    )
    μ_exact, Σ_exact = compute_coefficient_posterior(bsim, vcat(sources...))

    vsim = Simulation(medium, bd;
        solver = VariationalBayesianSolver(prior;
            learn_prior_flag = false,
            ard_prune_flag = false,
            max_iters = 1
        ),
        source_positions = sources,
        particular_solution = ParticularGravity(height = H)
    )
    vsol = solve(vsim)

    @test isapprox(vsol.fsol.coefficients, μ_exact; rtol = 1e-8)
    @test isapprox(vsol.fsol.coefficients_covariance, Σ_exact; rtol = 1e-8)
end

# ==============================================================================
# 3. Acoustic scattering (the example of acoustics_test.jl): complex Helmholtz
#    data on the unit circle, with more sources than measurements and ARD.
# ==============================================================================
@testset "acoustic scattering" begin
    Random.seed!(1234)
    ω = 2.0
    medium = Acoustic(2; ω = ω, ρ = 1.0, c = 1.0)
    x0, y0 = 0.0, 0.0
    ϕ(x, y) = -(im / 4) * hankelh1(0, (ω / medium.c) * sqrt((x - x0)^2 + (y - y0)^2))

    r = 1.0
    N_bd = 40
    θs = LinRange(0, 2pi, N_bd + 1)[1:N_bd]
    points = [[r * cos(θ), r * sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

    # known noise on each of the real and imaginary parts, carried by the boundary
    # data: each point is an MvNormal over the stacked [Re; Im] parts of its field
    σ = 1e-3
    fields = map(points) do p
        f = -ϕ(p[1], p[2]) + σ * (randn() + im * randn())
        MvNormal([real(f), imag(f)], σ^2 * I(2))
    end

    bd = BoundaryData(TractionType();
        boundary_points = points,
        fields = fields,
        normals = normals,
        interior_points = [[0.0, 0.0]]
    )

    # fresh boundary angles (not the measurement points) to test the prediction
    θ_test = LinRange(0, 2pi, 101)[1:100]
    field_scale = maximum(abs(ϕ(r * cos(θ), r * sin(θ))) for θ in θ_test)

    relative_prediction_error(vsol) = maximum(θ_test) do θ
        p = [r * cos(θ), r * sin(θ)]
        pred = field(TractionType(), vsol, p, [cos(θ), sin(θ)])[1]
        abs(pred - (-ϕ(p[1], p[2])))
    end / field_scale

    solver = VariationalBayesianSolver(
        prior_variance = 1.0,
        ard_threshold = 1e6,
        max_iters = 150,
        elbo_tol = 1e-9
    )

    @testset "built-in initial source positions" begin
        # take the sources from a denser sampling of the boundary, so that there are
        # more sources (and coefficients) than measurements
        N_dense = 80
        θd = LinRange(0, 2pi, N_dense + 1)[1:N_dense]
        bd_dense = BoundaryData(TractionType();
            boundary_points = [[r * cos(θ), r * sin(θ)] for θ in θd],
            normals = [[cos(θ), sin(θ)] for θ in θd],
            interior_points = [[0.0, 0.0]]
        )
        sources = source_positions(bd_dense; relative_source_distance = 2.0)
        @test length(sources) > N_bd   # more sources than measurements

        sim = Simulation(medium, bd; solver = solver, source_positions = sources)
        vsol = solve(sim)

        @test relative_prediction_error(vsol) < 0.05
        @test length(vsol.fsol.positions) < length(sources)   # ARD pruned sources
        @test elbo_is_monotone(vsol)
        @test 0.1 < vsol.misfit_ratio < 5.0

        # the posterior covariance is propagated for complex problems: field_covariance
        # returns the 2FD × 2FD covariance of the stacked [Re; Im] field, and its 3σ error
        # bars cover the true field at fresh boundary points
        @test size(field_covariance(TractionType(), vsol, [0.5, 0.0])) == (2, 2)
        covered = sum(θ_test) do θ
            p = [r * cos(θ), r * sin(θ)]
            pred = field(TractionType(), vsol, p, [cos(θ), sin(θ)])[1]
            stds = field_std(TractionType(), vsol, p, [cos(θ), sin(θ)])
            err = pred - (-ϕ(p[1], p[2]))
            (abs(real(err)) <= 3 * stds[1] + 1e-6) + (abs(imag(err)) <= 3 * stds[2] + 1e-6)
        end
        @test covered / (2 * length(θ_test)) >= 0.9
    end

    @testset "sources everywhere and ARD" begin
        sources = grid_source_positions(bd; n = 15, scale = 2.5, clearance = 1.0)
        @test length(sources) > N_bd   # more sources than measurements

        sim = Simulation(medium, bd; solver = solver, source_positions = sources)
        vsol = solve(sim)

        @test relative_prediction_error(vsol) < 0.05
        @test length(vsol.fsol.positions) < length(sources)
        @test elbo_is_monotone(vsol)
        @test 0.1 < vsol.misfit_ratio < 5.0
    end

    @testset "learning the source positions" begin
        N_small = 20
        θsmall = LinRange(0, 2pi, N_small + 1)[1:N_small]
        points_small = [[r * cos(θ), r * sin(θ)] for θ in θsmall]
        fields_small = map(points_small) do p
            f = -ϕ(p[1], p[2]) + σ * (randn() + im * randn())
            MvNormal([real(f), imag(f)], σ^2 * I(2))
        end
        bd_small = BoundaryData(TractionType();
            boundary_points = points_small,
            fields = fields_small,
            normals = [[cos(θ), sin(θ)] for θ in θsmall],
            interior_points = [[0.0, 0.0]]
        )

        n_src = 12
        θsrc = LinRange(0, 2pi, n_src + 1)[1:n_src]
        sources = [[2.5cos(θ), 2.5sin(θ)] for θ in θsrc]

        solver_chi = VariationalBayesianSolver(
            prior_variance = 1.0,
            optimise_source_positions_flag = true,
            source_position_iters = 3,
            ard_prune_flag = false,
            max_iters = 30,
            elbo_tol = 1e-9
        )
        sim = Simulation(medium, bd_small; solver = solver_chi, source_positions = sources)
        vsol = solve(sim)

        @test relative_prediction_error(vsol) < 0.1
        @test elbo_is_monotone(vsol)
    end
end

# ==============================================================================
# 4. Elastostatics on a circular domain (the Airy benchmark of benchmarks_test.jl),
#    with sources placed everywhere and ARD.
# ==============================================================================
@testset "Airy stress circular benchmark" begin
    Random.seed!(42)
    medium = Elastostatic(2; ρ = 1.0, cp = 2.0, cs = 1.0)

    # Airy stress solution, independent of r
    σrr(r, θ) = -2 * cos(2θ)
    σrθ(r, θ) = 2 * sin(2θ)

    R_bd = 1.3
    n_bd = 40
    θs = LinRange(0, 2pi, n_bd + 1)[1:n_bd]
    points = [[R_bd * cos(θ), R_bd * sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

    traction_true = [
        radial_to_cartesian_transform([R_bd, θ]) * [σrr(R_bd, θ), σrθ(R_bd, θ)]
    for θ in θs]

    σ_noise = 0.01 * maximum(norm.(traction_true))
    fields = [MvNormal(t .+ σ_noise .* randn(2), σ_noise^2 * I(2)) for t in traction_true]

    bd = BoundaryData(TractionType();
        boundary_points = points,
        fields = fields,
        normals = normals,
        interior_points = [[0.0, 0.0]]
    )

    sources = grid_source_positions(bd; n = 14, scale = 2.2, clearance = 1.0)
    @test 2 * length(sources) > 2 * n_bd   # more coefficients than measurements

    solver = VariationalBayesianSolver(
        prior_variance = 10.0^2,
        ard_threshold = 1e8,
        max_iters = 100,
        elbo_tol = 1e-9
    )
    sim = Simulation(medium, bd; solver = solver, source_positions = sources)
    # plot(sim; source_positions = true, boundary_points = true)
    vsol = solve(sim)

    # predict the traction on an interior circle r = 1 and compare with the analytic solution
    r_test = 1.0
    errors = Float64[]
    covered = 0
    total = 0
    for θ in θs
        p = [r_test * cos(θ), r_test * sin(θ)]
        nvec = [cos(θ), sin(θ)]
        truth = radial_to_cartesian_transform([r_test, θ]) * [σrr(r_test, θ), σrθ(r_test, θ)]
        pred = field(TractionType(), vsol, p, nvec)
        stds = field_std(TractionType(), vsol, p, nvec)

        push!(errors, norm(pred - truth))
        for c in 1:2
            total += 1
            covered += abs(pred[c] - truth[c]) <= 3 * stds[c] + 1e-6 ? 1 : 0
        end
    end

    @test maximum(errors) / mean(norm.(traction_true)) < 0.05
    @test covered / total >= 0.9   # the 3σ error bars must cover the analytic solution
    @test length(vsol.fsol.positions) < length(sources)
    @test elbo_is_monotone(vsol)
    @test 0.1 < vsol.misfit_ratio < 5.0
end

# ==============================================================================
# 5. Rectangle under gravity (the benchmark of benchmarks_test.jl), with the
#    built-in source initializer and a particular solution.
# ==============================================================================
@testset "gravity rectangle benchmark" begin
    Random.seed!(7)
    medium = Elastostatic(2; ρ = 1.0, cp = 3.0, cs = 2.0)

    n = 15; L = 0.5; H = 0.5
    rectangle_boundary(m) = begin
        xs = [
            LinRange(-L, L, m + 2)[2:end - 1]; zeros(m) .+ L;
            LinRange(L, -L, m + 2)[2:end - 1]; zeros(m) .- L
        ]
        ys = [
            zeros(m); LinRange(0, H, m + 2)[2:end - 1];
            zeros(m) .+ H; LinRange(H, 0, m + 2)[2:end - 1]
        ]
        pts = [[xs[i], ys[i]] for i in eachindex(xs)]
        nrm = [
            [[0.0, -1.0] for _ = 1:m]; [[1.0, 0.0] for _ = 1:m];
            [[0.0, 1.0] for _ = 1:m]; [[-1.0, 0.0] for _ = 1:m]
        ]
        pts, nrm
    end

    points, normals = rectangle_boundary(n)

    weight = 2L * H * medium.ρ * 9.81
    traction_true = [
        [[0.0, weight / (2L)] for _ = 1:n]; [[0.0, 0.0] for _ = 1:n];
        [[0.0, 0.0] for _ = 1:n]; [[0.0, 0.0] for _ = 1:n]
    ]

    σ_noise = 0.01 * weight / (2L)
    fields = [MvNormal(t .+ σ_noise .* randn(2), σ_noise^2 * I(2)) for t in traction_true]

    bd = BoundaryData(TractionType();
        boundary_points = points,
        fields = fields,
        normals = normals
    )

    # built-in initializer, applied to a denser sampling of the same boundary so that
    # there are more sources (and coefficients) than measurements
    points_dense, normals_dense = rectangle_boundary(25)
    bd_dense = BoundaryData(TractionType();
        boundary_points = points_dense,
        normals = normals_dense
    )
    sources = source_positions(bd_dense; relative_source_distance = 1.24)
    @test length(sources) > length(points)   # more sources than measurements

    ys_test = range(0, H, length = 50)
    f_true = [[0.0, medium.ρ * 9.81 * y - weight / (2L)] for y in ys_test]

    function slice_error_and_coverage(vsol)
        fs = [field(TractionType(), vsol, [0.0, y], [0.0, 1.0]) for y in ys_test]
        stds = [field_std(TractionType(), vsol, [0.0, y], [0.0, 1.0]) for y in ys_test]
        err = norm(norm.(fs .- f_true)) / norm(norm.(f_true))
        total = 2 * length(ys_test)
        covered = sum(eachindex(ys_test)) do i
            count(abs.(fs[i] .- f_true[i]) .<= 3 .* stds[i] .+ 1e-6)
        end
        return err, covered / total
    end

    @testset "fixed source positions with ARD" begin
        solver = VariationalBayesianSolver(
            prior_variance = 100.0^2,
            ard_threshold = 1e8,
            max_iters = 100,
            elbo_tol = 1e-9
        )
        sim = Simulation(medium, bd;
            particular_solution = ParticularGravity(height = H),
            solver = solver,
            source_positions = sources
        )
        vsol = solve(sim)

        err, coverage = slice_error_and_coverage(vsol)
        @test err < 0.05
        @test coverage >= 0.9
        @test length(vsol.fsol.positions) < length(sources)
        @test elbo_is_monotone(vsol)
        @test 0.1 < vsol.misfit_ratio < 5.0
    end

    @testset "learning the source positions" begin
        solver = VariationalBayesianSolver(
            prior_variance = 100.0^2,
            ard_threshold = 1e8,
            optimise_source_positions_flag = true,
            source_position_iters = 2,
            max_iters = 40,
            elbo_tol = 1e-9
        )
        sim = Simulation(medium, bd;
            particular_solution = ParticularGravity(height = H),
            solver = solver,
            source_positions = sources
        )
        vsol = solve(sim)

        err, coverage = slice_error_and_coverage(vsol)
        @test err < 0.05
        @test elbo_is_monotone(vsol)
    end
end

end
