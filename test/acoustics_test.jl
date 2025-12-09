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

apply(x) = real(x)
#apply(x) = imag(x)

medium = Acoustic(2; ω = ω,  ρ = 1.0, c = 1.0)

θs = LinRange(0,2pi,N_bd+1)[1:N_bd]
bd_points = [[r*cos(θ), r*sin(θ)] for θ in θs]
normals = [[cos(θ), sin(θ)] for θ in θs]
bd_fields = [[-ϕ(r*cos(θ), r*sin(θ))] for θ in θs]
interior_points = [[0.0, 0.0]]

bd = BoundaryData(TractionType(); 
    boundary_points = bd_points, 
    fields = bd_fields, 
    outward_normals = normals,
    interior_points = interior_points
)


source_pos=source_positions(bd; relative_source_distance = 1.0)

# Solve
solver = TikhonovSolver(λ=λ, tolerance = tolerance)
sim=Simulation(medium,bd, solver=solver, source_positions = source_pos)
fsol = solve(sim)

predict_fields = [field(TractionType(), fsol, bd_points[i], normals[i]) for i in eachindex(bd_points)]
fields = [-ϕ(r*cos(θ),r*sin(θ)) for θ in θs]
f=vcat(bd.fields...)

errors = [abs(fields[i] - predict_fields[i][1]) for i in eachindex(fields)]
@test maximum(errors) < 1e-10

end