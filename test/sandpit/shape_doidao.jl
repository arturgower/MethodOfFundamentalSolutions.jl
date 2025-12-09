using MethodOfFundamentalSolutions
using LinearAlgebra
using Test
using Statistics
using SpecialFunctions

N_bd = 100;
r = 1.0;
ω=6.0;
ϕ(x, y) = exp(im *( ω / medium.c)  *x )

#rsource=0.69;
N_sources=N_bd;

λ = 1e-8;
#λ = 0.0;
tolerance = 1e-11;
#tolerance = 0.0;

res = 101;

apply(x) = real(x)
#apply(x) = imag(x)

radius(θ)=1.5+0.3*sin(3*θ)+0.15*sin(5θ+0.6*sin(2θ))

medium = Acoustic(2; ω = ω, ρ = 1.0, c = 1.0)
θs = LinRange(0,2pi,N_bd+1)[1:N_bd]

bd_points = [[radius(θ)*cos(θ), radius(θ)*sin(θ)] for θ in θs]
normals = [[cos(θ), sin(θ)] for θ in θs]
bd_fields = [[-ϕ(radius(θ)*cos(θ), radius(θ)*sin(θ))] for θ in θs]
interior_points = [[0.0, 0.0]]
    
bd = BoundaryData(TractionType(); 
        boundary_points = bd_points, 
        fields = bd_fields, 
        outward_normals = normals,
        interior_points = interior_points
    )

bdplot=BoundaryData(TractionType(); 
        boundary_points = bd_points, 
        fields = real.(bd_fields), 
        outward_normals = normals,
        interior_points = interior_points
    )    
#Nsources=length(bds[n].boundary_points)
source_pos=source_positions(bd; relative_source_distance = 1.0)
rsource2(θ) = radius(θ)+ 0.2
θsource=LinRange(0,2pi,N_sources+1)[1:(N_sources)]
source_pos=[ [rsource2(θ)*cos(θ), rsource2(θ)*sin(θ)] for θ in θsource ]
# Solve
#solver = TikhonovSolver(tolerance = tolerance)
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
    
bottomleft = [-5.;-5.]
topright = [5.;5.]
region = Box([bottomleft, topright])
    
x_vec, inds = points_in_shape(region; res = res)
xs = x_vec[inds]

fs = [
    field(TractionType(), fsol, x, x / norm(x)) 
for x in xs];

field_mat = [[0.0+0.0im] for x in x_vec]
field_mat[inds] = [[fs[i][1]] for i in eachindex(fs)];
field_inc=[[ϕ(x[1],x[2])] for x in x_vec]
field_predict = FieldResult(x_vec, [real.(field_mat[i]+field_inc[i]) for i in eachindex(field_mat)]);
field_scat = FieldResult(x_vec, [apply.(field_mat[i]) for i in eachindex(field_mat)]);



maxc = maximum(apply, field_predict.field)[1]
minc = minimum(apply, field_predict.field)[1]
plot(field_predict, clim=(-1.0,1.0))
plot(field_predict, clim=(minc,maxc))
#scatter!(source_pos)
#maxc = maximum(field_scat.field)[1]
plot(field_scat, clim=(minc,maxc));
plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)
#plot!(fsol,color="gold")
# plot(bdplot)
#plot!(bd)
# err0r = [field_scat.field[i] .- apply.(result.field[:][i]) for i in eachindex(field_predict.field)]
# error_field=FieldResult(x_vec,err0r)
# plot(error_field, field_apply = norm, clim=(0.0,1.0))
# plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)