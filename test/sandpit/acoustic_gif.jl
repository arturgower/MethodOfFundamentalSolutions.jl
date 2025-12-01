using MethodOfFundamentalSolutions
using LinearAlgebra
using Test
using Statistics
using SpecialFunctions
#using MultipleScattering

N_bd = 100;
r = 1.0;
x0,y0 = 0.6,0.5;
ω=12.0;
#ϕ(x, y) = exp(im *( ω / medium.c)  *x )
ϕ(x, y) =hankelh1(0, (ω/medium.c) * sqrt((x-x0)^2 + (y-y0)^2))

rsource=1.8;
N_sources=100;

λ = 1e-8;
tolerance = 1e-10;

res = 50;

medium = Acoustic(2; ρ = 1.0, c = 1.0)
θs = LinRange(0,2pi,N_bd+1)[1:N_bd]

bd_points = [[r*cos(θ), r*sin(θ)] for θ in θs]
normals = [[cos(θ), sin(θ)] for θ in θs]
bd_fields = [[-ϕ(r*cos(θ), r*sin(θ))] for θ in θs]
interior_points = [[0.0, 0.0]]
    
bd = BoundaryData(DirichletType(); 
        boundary_points = bd_points, 
        fields = bd_fields, 
        outward_normals = normals,
        interior_points = interior_points
    )
    
#Nsources=length(bds[n].boundary_points)
θsource=LinRange(0,2pi,N_sources+1)[1:(N_sources)]
#source_pos=[ [rsource*cos(θ), rsource*sin(θ)] for θ in θsource ]
source_pos=source_positions(bd; relative_source_distance = 1.0)
# Solve
solver = TikhonovSolver(λ = λ, tolerance = tolerance)

sim=Simulation(medium,bd, solver=solver, source_positions = source_pos;ω=ω)

fsol = solve(sim)

predict_fields = [field(DirichletType(), fsol, bd_points[i], normals[i]; ω=ω) for i in eachindex(bd_points)]
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
    field(DirichletType(), fsol, x, x / norm(x)) 
for x in xs];
    
field_mat = [[0.0+0.0im] for x in x_vec]
field_mat[inds] = [[fs[i][1]] for i in eachindex(fs)];
field_inc=[[ϕ(x[1],x[2])] for x in x_vec]
field_predict = FieldResult(x_vec, [real.(field_mat[i]+field_inc[i]) for i in eachindex(field_mat)]);
field_scat = FieldResult(x_vec, [abs.(field_mat[i]) for i in eachindex(field_mat)]);

fs = map(xs) do x
    r, θ = cartesian_to_radial_coordinates(x)
    [ϕ(r*cos(θ),r*sin(θ))]
end

maxc = maximum(field_predict.field)[1]
minc = minimum(field_predict.field)[1]
#plot(field_predict, field_apply = first, clim=(-1.0,1.0));
plot(field_predict, field_apply = first, clim=(minc,maxc));
#plot(field_scat, field_apply = first)
scatter!([x0],[y0], markersize=5, markercolor=:green)
#plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)

ts = LinRange(0.,2pi/ω,30)
fields = [f[1] for f in field_predict.field]
field_complex = complex.(getindex.(field_predict.field, 1))

freq_result=FrequencySimulationResult(field_complex, field_predict.x, [ω])
anim = @animate for t in ts
    plot(freq_result,ω; seriestype = :heatmap,
        phase_time=t, clim=(minc,maxc) , c=:balance
    )
    scatter!([x0],[y0], markersize=5, markercolor=:green, label = "")
    plot!(colorbar=false, title="",axis=false, xguide ="",yguide ="")
end

gif(anim,"gap-diffraction.gif", fps = 7)





particle_shape = Circle(r)
particle = Particle(medium, particle_shape)
#source =  plane_source(medium; direction = [1.0,0.])
source =  point_source(Acoustic(2;ρ = 0.00001, c = 0.00001), [x0,y0], 1.0)
simulation = FrequencySimulation([particle], source)
result = run(simulation, region, ω, basis_order=5; res=res)
plot(result, ω; seriestype = :heatmap, field_apply = real, clim=(-1.0,1.0), title = "");
plot(result, ω; seriestype = :heatmap, field_apply = real, title = "");
plot!(Circle(r),fill = true, fillcolor = :gray, linecolor = :black)

err0r = [abs.(field_predict.field[i] .- result.field[:][i]) for i in eachindex(field_predict.field)]
error_field=FieldResult(x_vec,err0r)
   

#plot(error_field, field_apply = norm, clim=(0.0,1.0))
    x_vec
#plot!(bd)
norm.(x_vec)
hankelh1.(0,(ω/medium.c)*norm.(x_vec))

plot(x_vec,abs.(hankelh1.(0,(ω/medium.c)*norm.(x_vec))))