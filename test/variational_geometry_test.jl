using LinearAlgebra
using Test
using Statistics
using Random
using Distributions
using StaticArrays
using MethodOfFundamentalSolutions

# Tests for learning the boundary x in 2D with the VariationalBayesianSolver
# (update_geometry_flag = true, E-step II of Algorithm 1). The measured data are generated
# at the TRUE sensor positions, while the solver is given perturbed nominal positions with
# a diagonal Gaussian prior; it must infer the perturbation δx from the data.
#
# The tests deliberately cover both regimes:
#
# When it works:
#  - perturbations along the local field gradient (the data sees them through ∇u · δx), and
#    especially perturbations that are rough from sensor to sensor, which the smooth
#    fundamental-solution field cannot absorb;
#  - the posterior variance is a faithful identifiability diagnostic: it collapses in the
#    informed directions and stays at the prior in the uninformed ones.
#
# When it does not work:
#  - perturbations orthogonal to the field gradient carry no signal, and the estimate stays
#    at the nominal boundary (the prior);
#  - smooth low-order perturbations are partly absorbed by re-fitting the coefficients, so
#    they are only partially recovered;
#  - perturbations larger than the length scale of the field violate the linearization
#    eq. (linearization_M): the boundary is not recovered, while the mean-field posterior
#    still reports small uncertainty — a wrong and overconfident answer.
#
# These tests motivated two revisions to Algorithm 1 as implemented here: the prior of δx
# stays anchored at the original nominal boundary across re-centerings (otherwise the
# unidentified directions random-walk away with each re-linearization), and re-centering
# only counts as a model change when the step is significant relative to the prior scale
# (otherwise the bound can never be monitored and convergence is never declared).

# The bound is guaranteed to be monotone only between re-centerings of the boundary
# (Algorithm 1 re-linearizes the model there), so those iterations are excluded.
function geometry_elbo_is_monotone(vsol; tol = 1e-6)
    F = vsol.elbo_history
    return all(eachindex(F)) do i
        i == 1 || (i in vsol.baseline_resets) || (i in vsol.recenter_iterations) ||
            F[i] >= F[i - 1] - tol * (1 + abs(F[i - 1]))
    end
end

@testset "Learning the boundary" begin

    n_bd = 40
    θs = LinRange(0, 2pi, n_bd + 1)[1:n_bd]
    x_true_circle = [SVector(cos(θ), sin(θ)) for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

    n_src = 16
    θsrc = LinRange(0, 2pi, n_src + 1)[1:n_src]
    sources = [[2.0cos(θ), 2.0sin(θ)] for θ in θsrc]

    σ_noise = 1e-4

    # Solve the Laplace Dirichlet problem with data g measured at the true points but
    # nominal (perturbed) sensor positions, learning the boundary.
    function solve_laplace_geometry(x_nom, g, σx; max_iters = 100)
        medium = LaplaceMedium{2, Float64}()
        bd = BoundaryData(DirichletType();
            boundary_points = MvNormal(vcat(Vector.(x_nom)...), σx^2 * I(2n_bd)),
            fields = [MvNormal([gi], σ_noise^2 * I(1)) for gi in g],
            normals = normals,
            interior_points = [[0.0, 0.0]]
        )
        solver = VariationalBayesianSolver(
            prior_variance = 10.0^2,
            update_geometry_flag = true,
            ard_prune_flag = false,
            max_iters = max_iters,
            elbo_tol = 1e-10
        )
        sim = Simulation(medium, bd; solver = solver, source_positions = sources)
        return solve(sim)
    end

    coordinate_error(xs, x_true, c) = norm([xs[i][c] - x_true[i][c] for i in eachindex(x_true)])
    # mean posterior standard deviation of the boundary coordinate c = 1 (x) or 2 (y)
    boundary_std(vsol, c) = mean(sqrt.(diag(cov(vsol.boundary_shape.boundary_points))[c:2:end]))

    # ==========================================================================
    # Works: rough sensor shifts along the field gradient. The true field is
    # u = x, whose gradient (1,0) makes x-shifts visible in the data, and an
    # alternating ±σx shift cannot be absorbed by the smooth field.
    # ==========================================================================
    @testset "recovers rough perturbations along the field gradient" begin
        Random.seed!(101)
        σx = 0.05
        ε = [σx * (isodd(i) ? 1.0 : -1.0) for i in 1:n_bd]
        x_nom = [x_true_circle[i] + SVector(ε[i], 0.0) for i in 1:n_bd]
        g = [p[1] + σ_noise * randn() for p in x_true_circle]   # u = x at the true points

        vsol = solve_laplace_geometry(x_nom, g, σx)
        x_est = mean_points(vsol.boundary_shape)

        err_nom = coordinate_error(x_nom, x_true_circle, 1)
        err_est = coordinate_error(x_est, x_true_circle, 1)

        @test err_est < 0.25 * err_nom                      # strong recovery
        @test coordinate_error(x_est, x_true_circle, 2) < σx  # no drift in the blind y-direction
        @test boundary_std(vsol, 1) < 0.1 * σx              # x is reported as well determined
        @test boundary_std(vsol, 2) > 0.9 * σx              # y is reported as prior-limited
        @test geometry_elbo_is_monotone(vsol)
        @test 0.1 < vsol.misfit_ratio < 5.0
    end

    # ==========================================================================
    # Partial: smooth random shifts along the gradient. The smooth part of the
    # perturbation can be traded against a re-fitted smooth field, so only the
    # rough part is identifiable: the error is reduced but NOT to the noise level.
    # ==========================================================================
    @testset "smooth perturbations are only partly recovered" begin
        Random.seed!(102)
        σx = 0.05
        ε = σx .* randn(n_bd)
        x_nom = [x_true_circle[i] + SVector(ε[i], 0.0) for i in 1:n_bd]
        g = [p[1] + σ_noise * randn() for p in x_true_circle]

        vsol = solve_laplace_geometry(x_nom, g, σx)
        x_est = mean_points(vsol.boundary_shape)

        err_nom = coordinate_error(x_nom, x_true_circle, 1)
        err_est = coordinate_error(x_est, x_true_circle, 1)

        @test err_est < 0.6 * err_nom     # the boundary does improve ...
        @test err_est > 0.1 * err_nom     # ... but the smooth residual is absorbed by the field
        @test geometry_elbo_is_monotone(vsol)
    end

    # ==========================================================================
    # Does not work: shifts orthogonal to the field gradient. With u = x the data
    # do not depend on the y-coordinates at all, so the estimate must stay at the
    # nominal boundary — and the posterior variance must say so, staying at the prior.
    # ==========================================================================
    @testset "perturbations orthogonal to the gradient are unidentifiable" begin
        Random.seed!(103)
        σx = 0.05
        ε = σx .* randn(n_bd)
        x_nom = [x_true_circle[i] + SVector(0.0, ε[i]) for i in 1:n_bd]
        g = [p[1] + σ_noise * randn() for p in x_true_circle]

        vsol = solve_laplace_geometry(x_nom, g, σx)
        x_est = mean_points(vsol.boundary_shape)

        err_nom = coordinate_error(x_nom, x_true_circle, 2)
        err_est = coordinate_error(x_est, x_true_circle, 2)

        @test abs(err_est - err_nom) < 0.1 * err_nom   # no (false) recovery, no drift
        @test boundary_std(vsol, 2) > 0.9 * σx         # the algorithm reports: no information in y
        @test geometry_elbo_is_monotone(vsol)
    end

    # ==========================================================================
    # Does not work: perturbations beyond the linearization. The field of a source
    # at distance 0.35 from the boundary varies on that length scale, so shifts of
    # size 0.25 violate the first-order expansion eq. (linearization_M): the
    # boundary is not recovered, and — the dangerous part — the mean-field
    # posterior still reports a variance well below the prior, i.e. the answer is
    # wrong and overconfident. With small shifts the same field is recovered.
    # ==========================================================================
    @testset "fails beyond the linearization regime" begin
        Random.seed!(104)
        s_near = SVector(1.35, 0.0)
        u_near(p) = -1 / (2π) * log(norm(p - s_near))
        g = [u_near(p) + σ_noise * randn() for p in x_true_circle]

        total_error(xs) = norm(vcat([Vector(xs[i] - x_true_circle[i]) for i in 1:n_bd]...))

        # large radial shifts, comparable to the distance to the near source
        σx_big = 0.25
        ε = σx_big .* randn(n_bd)
        x_nom = [x_true_circle[i] + ε[i] * SVector(cos(θs[i]), sin(θs[i])) for i in 1:n_bd]

        vsol = solve_laplace_geometry(x_nom, g, σx_big)

        @test total_error(mean_points(vsol.boundary_shape)) > 0.7 * total_error(x_nom)  # not recovered
        # overconfidence: the reported uncertainty shrinks although the estimate is wrong
        @test boundary_std(vsol, 1) < 0.5 * σx_big
        @test geometry_elbo_is_monotone(vsol)

        # the same nonlinear field with small shifts is (partly) recovered
        σx_small = 0.02
        ε_small = σx_small .* randn(n_bd)
        x_nom_small = [x_true_circle[i] + ε_small[i] * SVector(cos(θs[i]), sin(θs[i])) for i in 1:n_bd]

        vsol_small = solve_laplace_geometry(x_nom_small, g, σx_small)

        @test total_error(mean_points(vsol_small.boundary_shape)) < 0.95 * total_error(x_nom_small)
        @test geometry_elbo_is_monotone(vsol_small)
    end

    # ==========================================================================
    # Elastostatics: traction data from the Airy solution σrθ = 6 r² sin(2θ),
    # whose r-dependence makes radial sensor shifts observable (except near the
    # zeros of sin(2θ)): the radial error is reduced, the tangential one is not.
    # ==========================================================================
    @testset "elastostatic traction data" begin
        Random.seed!(105)
        medium = Elastostatic(2; ρ = 1.0, cp = 2.0, cs = 1.0)
        traction_at(r, θ) = radial_to_cartesian_transform([r, θ]) * [0.0, 6 * r^2 * sin(2θ)]

        σx = 0.03
        σ_noise_e = 6e-3
        ε = σx .* randn(n_bd)
        x_true_e = [SVector((1 + ε[i]) * cos(θs[i]), (1 + ε[i]) * sin(θs[i])) for i in 1:n_bd]
        x_nom = x_true_circle   # the nominal boundary is the unperturbed unit circle
        g = [traction_at(1 + ε[i], θs[i]) .+ σ_noise_e .* randn(2) for i in 1:n_bd]

        bd = BoundaryData(TractionType();
            boundary_points = MvNormal(vcat(Vector.(x_nom)...), σx^2 * I(2n_bd)),
            fields = [MvNormal(gi, σ_noise_e^2 * I(2)) for gi in g],
            normals = normals,
            interior_points = [[0.0, 0.0]]
        )
        solver = VariationalBayesianSolver(
            prior_variance = 10.0^2,
            update_geometry_flag = true,
            ard_prune_flag = false,
            max_iters = 100,
            elbo_tol = 1e-10
        )
        sim = Simulation(medium, bd; solver = solver, source_positions = sources)
        vsol = solve(sim)
        x_est = mean_points(vsol.boundary_shape)

        radial_error(xs) = norm([dot(xs[i] - x_true_e[i], SVector(cos(θs[i]), sin(θs[i]))) for i in 1:n_bd])

        @test radial_error(x_est) < 0.8 * radial_error(x_nom)
        @test geometry_elbo_is_monotone(vsol)
        @test 0.1 < vsol.misfit_ratio < 5.0
    end

end
