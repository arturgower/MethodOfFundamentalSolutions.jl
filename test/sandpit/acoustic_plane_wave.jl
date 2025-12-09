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

res = 51;

apply(x) = real(x)
#apply(x) = imag(x)

medium = Acoustic(2; ρ = 1.0, c = 1.0)
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
    
#Nsources=length(bds[n].boundary_points)
source_pos=source_positions(bd; relative_source_distance = 1.0)
rsource = r-(source_pos[1][1]-r)
θsource=LinRange(0,2pi,N_sources+1)[1:(N_sources)]
source_pos=[ [rsource*cos(θ), rsource*sin(θ)] for θ in θsource ]
# Solve
#solver = TikhonovSolver(tolerance = tolerance)
solver = TikhonovSolver(λ=λ, tolerance = tolerance)

sim=Simulation(medium,bd, solver=solver, source_positions = source_pos;ω=ω)

fsol = solve(sim)

predict_fields = [field(TractionType(), fsol, bd_points[i], normals[i]; ω=ω) for i in eachindex(bd_points)]
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
    field(TractionType(), fsol, x, x / norm(x); ω=ω) 
for x in xs];
   


field_mat = [[0.0+0.0im] for x in x_vec]
field_mat[inds] = [[fs[i][1]] for i in eachindex(fs)];
field_inc=[[ϕ(x[1],x[2])] for x in x_vec]
#field_predict = FieldResult(x_vec, [real.(field_mat[i]+field_inc[i]) for i in eachindex(field_mat)]);
field_scat = FieldResult(x_vec, [apply.(field_mat[i]) for i in eachindex(field_mat)]);

fs = map(xs) do x
    r, θ = cartesian_to_radial_coordinates(x)
    [ϕ(r*cos(θ),r*sin(θ))]
end



particle_shape = Circle(r)
particle = Particle(Acoustic(2; ρ = 0.0, c = 0.0), particle_shape)
source =  plane_source(medium; direction = [1.0,0.])
simulation = FrequencySimulation([particle], source)
result = run(simulation, region, ω, basis_order=7; only_scattered_waves = true, res=res)
#plot(result, ω; seriestype = :heatmap, field_apply = real, clim=(-1.0,1.0), title = "");
plot(result, ω; seriestype = :heatmap, field_apply = apply, title = "");
plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)

maxc = maximum(apply, result.field)[1]
minc = minimum(apply, result.field)[1]
#plot(field_predict, clim=(-1.0,1.0));
#plot(field_predict, clim=(minc,maxc));
#scatter!(source_pos)
#maxc = maximum(field_scat.field)[1]
plot(field_scat, clim=(minc,maxc));
plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)

err0r = [field_scat.field[i] .- apply.(result.field[:][i]) for i in eachindex(field_predict.field)]
error_field=FieldResult(x_vec,err0r)
plot(error_field, field_apply = norm, clim=(0.0,1.0))
plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)