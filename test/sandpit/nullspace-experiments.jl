using MethodOfFundamentalSolutions
using MultipleScattering
using StaticArrays:SVector
using LinearAlgebra
using Statistics
using BlockArrays
medium = Elastostatic(2; ρ = 1.0, cp = 2.0, cs = 1.0)
    θs_arr = [LinRange(0,2pi,n)[1:(n-1)] for n in (10,20,40,80,140)];
    i=5
    θs = θs_arr[i]
    r = 1.3
    r = 1.0
    points = [[r*cos(θ), r*sin(θ)] for θ in θs]
    normals = [[cos(θ), sin(θ)] for θ in θs]

       # Circumferential stress
    # From the Airy stress function we have the solution inside a circular domain which does not depend on the radius r: 
    σrr(r, θ) = -2 * cos(2θ)
    σθθ(r, θ) = 2 * cos(2θ)
    σrθ(r, θ) = 2 * sin(2θ)
    bds = map(θs_arr) do θs
        points = [[r*cos(θ), r*sin(θ)] for θ in θs]
        normals = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]

        # WRRRROOONG needs the basis vectors
        fields = [radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)] for θ in θs]
        BoundaryData(TractionType(); 
            boundary_points = points, 
            fields = fields, 
            normals = normals,
            interior_points = interior_points
        )
    end
    n=i
    #Nsources=length(bds[n].boundary_points)
    rsource=1.0051989236591015
    N_sources=1081
    θsource=LinRange(0,2pi,N_sources)[1:(N_sources-1)]
    source_pos=[ [rsource*cos(θ), rsource*sin(θ)] for θ in θsource ]

    # Solve
    solver = TikhonovSolver(tolerance = 1e-8)
    fsols = map(bds) do bd
        fsol = solve(medium, bd; 
            solver=solver, 
         #   source_positions = source_positions(bd; relative_source_distance = 1.0) 
            source_positions = source_pos 
        )
    end

    
    
    fsol=fsols[n]
    #f
    predict_fields = [field(TractionType(), fsol, points[i], normals[i]) for i in eachindex(points)]
    fields = [radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)] for θ in θs]

    Msystem=system_matrix(fsol,bds[n])
    typeof(Msystem)
    Msystem=Matrix(Msystem)
    
    F=svd(Msystem)
    U=F.U
    S=F.S
    Vt=F.Vt
    
    S
    
    f=vcat(bds[n].fields...)
    
    errors = [norm(fields[i] - predict_fields[i]) for i in eachindex(fields)]
    errors = errors ./ mean([norm(f) for f in fields])
    
    maximum(errors)
    bd=bds[n]
    using Plots
    x_vec, inds = points_in_shape(bd; res = 15)
    x_vec, inds = points_in_shape(bd; res = 25)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.
    xs = x_vec[inds]

    fs = [
        field(TractionType(), fsol, x, x / norm(x)) 
    for x in xs];
    
    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_predict = FieldResult(x_vec, field_mat[:]);

     # fs = [field(wave,x,fieldtype) for x in xs];
    fs = map(xs) do x
        r, θ = cartesian_to_radial_coordinates(x)
       radial_to_cartesian_transform(SVector(r,θ))*[σrr(r,θ), σrθ(r,θ)]
    end

    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_true = FieldResult(x_vec, field_mat[:]);

    using Plots
    plot(field_predict, field_apply = first, c = :inferno)
    plot(field_true, field_apply = first, c = :inferno)
    
    plot(field_predict - field_true, field_apply = norm, c = :inferno)
    plot!(fsol, xlims = (-1.6,1.6), ylims = (-1.6,1.6))
    plot!(bd)
    
    M = system_matrix(fsol, bd)
    cond(M)
    svdM = svd(M)
    svdM.S[end-6:end]

    norm(M * svdM.Vt[end,:] - svdM.S[end] .* svdM.U[:,end]) / norm(M * svdM.Vt[end,:])

    fsol_null = FundamentalSolution(fsol.medium; 
        positions = fsol.positions, 
        coefficients = svdM.Vt[end,:]
    )

    fs = [
        field(TractionType(), fsol_null, x, x / norm(x)) 
    for x in xs];
    
    field_mat = [[0.0, 0.0] for x in x_vec]
    field_mat[inds] = fs;
    field_null = FieldResult(x_vec, field_mat[:]);

    plot(field_null, field_apply = norm, c = :inferno)