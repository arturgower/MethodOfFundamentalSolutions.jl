using MethodOfFundamentalSolutions
using LinearAlgebra
using Test
using Statistics
using SpecialFunctions
#using MultipleScattering

N_bd = 100;
r = 1.0;
x0,y0 = 0.0,0.0;
ω=2.0;
#ϕ(x, y) = exp(im *( ω / medium.c)  *x )
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
maximum(errors)
    
using Plots
using MultipleScattering
    
bottomleft = [-1.;-1.]
topright = [1.;1.]
region = Box([bottomleft, topright])
    
x_vec, inds = points_in_shape(bd; res = res)
xs = x_vec[inds]

fs = [
    field(TractionType(), fsol, x, x / norm(x)) 
for x in xs];
    
field_mat = [[0.0+0.0im] for x in x_vec]
field_mat[inds] = [[fs[i][1]] for i in eachindex(fs)];
field_inc=[[ϕ(x[1],x[2])] for x in x_vec]
#field_predict = FieldResult(x_vec, [apply.(field_mat[i]+field_inc[i]) for i in eachindex(field_mat)]);
field_scat = FieldResult(x_vec, [apply.(field_mat[i]) for i in eachindex(field_mat)]);

fs = map(xs) do x
    r, θ = cartesian_to_radial_coordinates(x)
    [ϕ(r*cos(θ),r*sin(θ))]
end

maxc = maximum(field_scat.field)[1]
minc = minimum(field_scat.field)[1]
#plot(field_predict, field_apply = first, clim=(-1.0,1.0));
#plot(field_predict, clim=(minc,maxc));
plot(field_scat)
scatter!([x0],[y0], markersize=5, markercolor=:green)
#plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)

ϕ_exact(x, y) = (im/4)*(hankelh1(0, (ω/medium.c) * r)/besselj0((ω/medium.c) * r)) * besselj0((ω/medium.c) * sqrt((x)^2 + (y)^2))
field_exact=[[ϕ_exact(x[1],x[2])] for x in x_vec]
field_error = FieldResult(x_vec, [apply.(field_exact[i]-field_mat[i]) for i in eachindex(field_exact)])

plot(field_error, field_apply = first)


