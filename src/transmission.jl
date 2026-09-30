# Transmission (penetrable-scatterer) problems.
#
# A closed interface separates two regions, each modelled by its own [`Simulation`](@ref):
# `sim1` (e.g. the exterior, carrying the incident field and the unknown scattered field) and
# `sim2` (e.g. the interior, carrying the unknown transmitted field). The two share the same
# interface geometry and the same boundary conditions (the tuple of `FieldType`s of their
# `boundary_data`), but each has its own medium and its own MFS `source_positions`.
#
# Because `boundary_data` may itself be a tuple of `BoundaryData`, `system_matrix(sim)` already
# returns the full stacked block-column for one region across every condition, so the coupled
# system is simply
#     [ system_matrix(sim1)  -system_matrix(sim2) ] * [x1; x2] = boundary_forcing(sim1) - boundary_forcing(sim2),
# i.e. the residual `system_matrix(sim)*x - boundary_forcing(sim)` (the region's fundamental
# field minus its known/particular data) is matched across the interface. Any known excitation
# — an incident wave — enters through the region's boundary `fields` or its
# `particular_solution` (see [`PointSource`](@ref)), exactly as for a single-domain solve.

"""
    PointSource(positions, amplitudes) <: ParticularSolution

A known field radiated by point sources at `positions` with the given `amplitudes`, in
whatever medium it is evaluated against. Its trace is `∑ᵢ amplitudes[i] * greens(F, medium,
x - positions[i], n)`, so it can be used as the `particular_solution` of a [`Simulation`](@ref)
to represent an incident field (e.g. the incident wave of a transmission problem).
"""
struct PointSource{Dim, T} <: ParticularSolution
    positions::Vector{SVector{Dim, Float64}}
    amplitudes::Vector{T}
end

function PointSource(positions, amplitudes)
    Dim = length(first(positions))
    return PointSource([SVector{Dim, Float64}(p) for p in positions], collect(amplitudes))
end

function field(field_type::FieldType, medium::PhysicalMedium, ps::PointSource, x::AbstractVector, outward_normal::AbstractVector)
    x = SVector(x...)
    n = SVector(outward_normal...) ./ norm(outward_normal)
    Gs = [greens(field_type, medium, x - p, n) for p in ps.positions]
    return hcat(Gs...) * ps.amplitudes[:]
end

"""
    TransmissionSimulation(sim1, sim2; solver = sim1.solver)

An MFS transmission problem coupling two [`Simulation`](@ref)s across their shared interface:
`sim1` and `sim2` must use the same boundary conditions (the `FieldType`s of their
`boundary_data`) on the same interface points, but each carries its own medium and MFS
`source_positions`. Calling [`solve`](@ref) enforces continuity of the region residuals across
the interface and returns the two [`FundamentalSolution`](@ref)s `(fsol1, fsol2)`.

The two simulations are expected to share the same `solver`; the coupled solve uses `solver`
(the shared one by default).
"""
struct TransmissionSimulation{S <: AbstractSolver, SIM1 <: Simulation, SIM2 <: Simulation}
    sim1::SIM1
    sim2::SIM2
    solver::S
end

TransmissionSimulation(sim1::Simulation, sim2::Simulation; solver::AbstractSolver = sim1.solver) =
    TransmissionSimulation(sim1, sim2, solver)

"""
    solve(tsim::TransmissionSimulation)

Solve the coupled transmission problem, returning the two [`FundamentalSolution`](@ref)s
`(fsol1, fsol2)`. The coupled system
```
[ system_matrix(sim1)  -system_matrix(sim2) ] * [x1; x2] = boundary_forcing(sim1) - boundary_forcing(sim2)
```
is assembled here (solver-independent); the actual solve of `M * x = forcing` is delegated to
[`transmission_coefficients`](@ref), which dispatches on `tsim.solver`. The returned solutions
therefore follow the solver: a deterministic [`TikhonovSolver`](@ref) gives plain coefficients,
while a probabilistic solver returns a posterior covariance carried by each `FundamentalSolution`.
"""
function solve(tsim::TransmissionSimulation)

    sim1, sim2 = tsim.sim1, tsim.sim2

    M1 = system_matrix(sim1)
    M2 = system_matrix(sim2)
    M = hcat(M1, -M2)

    forcing = boundary_forcing(sim1) - boundary_forcing(sim2)

    coes, covariance, relative_error = transmission_coefficients(tsim.solver, M, forcing)

    # the first columns are sim1's unknowns, the rest are sim2's
    n1 = size(M1, 2)
    cov1, cov2 = split_covariance(covariance, n1)

    fsol1 = FundamentalSolution(sim1.medium;
        positions = sim1.source_positions,
        coefficients = coes[1:n1],
        coefficients_covariance = cov1,
        particular_solution = sim1.particular_solution,
        relative_boundary_error = relative_error
    )
    fsol2 = FundamentalSolution(sim2.medium;
        positions = sim2.source_positions,
        coefficients = coes[n1+1:end],
        coefficients_covariance = cov2,
        particular_solution = sim2.particular_solution,
        relative_boundary_error = relative_error
    )

    return fsol1, fsol2
end

"""
    transmission_coefficients(solver, M, forcing) -> (coefficients, covariance, relative_error)

Solve the assembled coupled system `M * coefficients = forcing` with `solver`. This is the one
method a solver must provide to work with [`solve(::TransmissionSimulation)`](@ref); the
returned `covariance` is `0.0 * I` for a deterministic solve and the joint posterior covariance
for a probabilistic one.
"""
function transmission_coefficients(solver::TikhonovSolver, M::AbstractMatrix, forcing::AbstractVector)
    coes, relative_error, condM = tikhonov_solve(M, forcing, solver)
    @info "Solved the transmission system with condition number $(condM), and with a relative error of the boundary data of $(relative_error), using the tolerance $(solver.tolerance)"
    return coes, 0.0 * I, relative_error
end

# split a joint coefficient covariance into the sim1 / sim2 diagonal blocks. Only the
# deterministic (UniformScaling) case is needed today; a probabilistic transmission solver
# should add the method matching its covariance convention (real K×K vs complex 2K×2K).
split_covariance(cov::UniformScaling, n1::Int) = (cov, cov)
