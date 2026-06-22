@testset "boundarydata" begin

    # Use points sampled on the boundary with different resolutions
    n = 15
    θs = LinRange(0,2pi,n)[1:(n-1)];

    r = 1.3
    ε = 0.0

    # Create boundary data with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]

        outward_ns2 = compute_outward_normals(points, interior_points)
        errors = [norm(outward_ns[i] - outward_ns2[i]) for i in eachindex(outward_ns)]

    @test maximum(errors) < 1e-10

    ε = 0.1
    # Create point cloud with errors
        points = [[r*cos(θ) + (rand() -0.5) * ε, r*sin(θ)+ (rand() -0.5) * ε] for θ in θs]
        outward_ns = [[cos(θ), sin(θ)] for θ in θs]
        interior_points = [[0.0, 0.0]]
        
        outward_ns2 = compute_outward_normals(points, interior_points)
        errors = [norm(outward_ns[i] - outward_ns2[i]) for i in eachindex(outward_ns)]

    @test maximum(errors) < 0.1
   
    # Let us check whether the interior is correctly defined.
    points = [[r*cos(θ), r*sin(θ)] for θ in θs]
    cloud = BoundaryData(DisplacementType();
        boundary_points = points,
        interior_points = interior_points
    )

    x_vec, inds = points_in_shape(cloud; res = 40)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.
    xs = x_vec[inds]

    x1s = [x[1] for x in x_vec] |> unique;
    x2s = [x[2] for x in x_vec] |> unique;
    dx = abs(x1s[2] - x1s[1])
    dy = abs(x2s[2] - x2s[1])

    # zs = zeros(x_vec |> length) 
    # zs[inds] .= 1.0
    # zs = reshape(zs, (length(x1s), length(x2s))) |> transpose;
    # heatmap(x1s |> sort,x2s |> sort,zs)
    # plot!(cloud)

    @test all([norm(x) < r + norm([dx,dy]) for x in xs]) 

    # rough estimate of volume
    volume = dx*dy * length(xs)
    exact_volume = pi * r^2

    @test abs(volume - exact_volume) < dx

    fields = [
        [0.0, 0.0]
    for x in x_vec]

    res = FieldResult(x_vec, fields)

    # using Plots
    # plot(res, clims = (-1,1))
    # plot!(cloud)
end