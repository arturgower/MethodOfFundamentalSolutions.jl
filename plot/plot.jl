# recipe for plotting a PointCloud. Assumes 2D boundary_points 
@recipe function plot(cloud::BoundaryData; fields = false, outward_normals = true, interior_points = true)
    # Set default attributes
    markershape --> :circle
    legend --> false
    aspect_ratio --> 1.0

    bs = cloud.boundary_points
    ns = cloud.outward_normals
    is = cloud.interior_points
    fs = cloud.fields

    # used to determine the length to plot the vectors
    lengthscale = mean([norm(b - mean(bs)) for b in bs])

    bx = [b[1] for b in bs]
    by = [b[2] for b in bs]
    ix = [b[1] for b in is]
    iy = [b[2] for b in is]

    fx = [b[1] for b in fs]
    fy = [b[2] for b in fs]

    f_lengthscale = mean([norm(f) for f in fs]) > 0 ? mean([norm(f) for f in fs]) : 1.0
    fx = fx .* (lengthscale/4) ./ (sqrt(2.0) * f_lengthscale)
    fy = fy .* (lengthscale/4) ./ (sqrt(2.0) * f_lengthscale)

    nx = [b[1] for b in ns]
    ny = [b[2] for b in ns]

    nx = nx .* (lengthscale/4)
    ny = ny .* (lengthscale/4)

    # fallback: repeat or zero-length normals
    if length(nx) != length(bx)
        nx = zeros(length(bx))
        ny = zeros(length(by))
    end

    @series begin
        label --> "Boundary points"
        seriestype --> :scatter
        markersize --> 4.0
        (bx, by)
    end

    if interior_points
    @series begin
        label --> "Interior points"
        seriestype --> :scatter
        markersize --> 6.0
        (ix, iy)
    end
    end

    if outward_normals
        @series begin
            label --> ""
            color --> :black
            seriestype --> :quiver
            quiver --> (nx, ny)
            marker = :none           # remove marker symbol at arrow tips
            markersize --> 0.0          # fallback to ensure no marker is drawn
            (bx, by)
        end
    end

    if fields
        @series begin
            label --> ""
            color --> :green
            seriestype --> :quiver
            quiver --> (fx, fy)
            marker = :none           # remove marker symbol at arrow tips
            markersize --> 0.0          # fallback to ensure no marker is drawn
            (bx, by)
        end
    end
end

@recipe function plot(sim::Simulation; source_positions = true)
    # Set default attributes
    legend --> false

    # Extract source and evaluation points
    sources = sim.source_positions

    # Split into x and y coordinates
    sx = [s[1] for s in sources]
    sy = [s[2] for s in sources]

    # Plot source points
    if source_positions
        @series begin
            label --> "Source points"
            seriestype --> :scatter
            markersize --> 4.0
            color --> :red
            (sx, sy)
        end
    end

    # Plot boundary data
    @series sim.boundary_data
end

@recipe function plot(fsol::FundamentalSolution)
    # Set default attributes
    legend --> false

    # Extract source and evaluation points
    sources = fsol.positions

    # Split into x and y coordinates
    sx = [s[1] for s in sources]
    sy = [s[2] for s in sources]

    # Plot source points
    @series begin
        label --> "Source points"
        seriestype --> :scatter
        markersize --> 4.0
        color --> :red
        (sx, sy)
    end
end



# Plot the result in space (across all x) for a specific angular frequency
@recipe function plot(res::FieldResult;
        region_shape = :empty, field_apply := first)

    x = [x[1] for x in res.x]
    y = [x[end] for x in res.x] # y will actually be z for 3D...

    seriestype --> :heatmap
    seriescolor --> :balance
    aspect_ratio --> 1.0

    st = get(plotattributes, :seriestype, :surface)
    
    if st == :heatmap
        # We could check here to see if x and y have the right structure
        x = unique(x)
        y = unique(y)

        n_x = length(x)
        n_y = length(y)

        fill --> true

        if region_shape != :empty
            bounds = bounding_box(region_shape)

            # If user has not set xlims and ylims, set them to the rectangle
            xlims --> (bottomleft(bounds)[1], topright(bounds)[1])
            ylims --> (bottomleft(bounds)[2], topright(bounds)[2])
        else
            xlims --> (minimum(x), maximum(x))
            ylims --> (minimum(y), maximum(y))
        end

        x, y, field_apply.(transpose(reshape(field(res),n_x,n_y)))

    else
        (x, y, field_apply.(field(res)))
    end

end

function plot_reconstructed_fields(
    points_flat, 
    mean_field, 
    variance_field, 
    boundary_positions_flat;          # Semicolon marks the start of keyword arguments
    source_positions_flat = nothing   # Defaults to nothing if omitted
)
    
    # 1. Extract X and Y coordinates using your reshaping function
    pos_matrix = flat_to_pos_matrix(points_flat)
    x_pts, y_pts = pos_matrix[1, :], pos_matrix[2, :]

    boundary_positions_matrix = flat_to_pos_matrix(boundary_positions_flat)
    x_boundary, y_boundary = boundary_positions_matrix[1, :], boundary_positions_matrix[2, :]

    # 2. Plot Left: Reconstructed Mean Field
    p1 = scatter(x_pts, y_pts, 
                 zcolor=mean_field, 
                 colorbar=true,
                 cmap=:balance,              
                 markersize=5,
                 markerstrokewidth=0,        
                 label="",
                 aspect_ratio=:equal,
                 xlabel="X", ylabel="Y",
                 title="Mean Field Intensity",
                 legend=:outertop
                 
                 )

    # Overlay the boundary line on Left plot
    plot!(p1, x_boundary, y_boundary, color=:black, linestyle=:dash, linewidth=1.5, label="r = 1")

    # 3. Plot Right: Reconstructed Variance Field (Uncertainty)
    p2 = scatter(x_pts, y_pts, 
                 zcolor=variance_field, 
                 colorbar=true,
                 cmap=:inferno,              
                 markersize=5,
                 markerstrokewidth=0, 
                 label="",
                 aspect_ratio=:equal,
                 xlabel="X", ylabel="Y",
                 title="Field Variance (Uncertainty)",
                 legend=:outertop)

    # Overlay the boundary line on Right plot
    plot!(p2, x_boundary, y_boundary, color=:black, linestyle=:dash, linewidth=1.5, label="r = 1")

    # 4. Optional: Overlay Sources only if they are provided
    if source_positions_flat !== nothing && !isempty(source_positions_flat)
        source_positions_matrix = flat_to_pos_matrix(source_positions_flat)
        x_sources, y_sources = source_positions_matrix[1, :], source_positions_matrix[2, :]
        
        # Explicitly target p1
        scatter!(p1, x_sources, y_sources, 
                 label="Sources", 
                 color=:red, 
                 markersize=6, 
                 marker=:circle) 
        
        # Explicitly target p2
        scatter!(p2, x_sources, y_sources, 
                 label="Sources", 
                 color=:red, 
                 markersize=6, 
                 marker=:circle) 
    end

    # 5. Combine plots side-by-side
    full_layout = plot(p1, p2, layout=(1, 2), size=(1150, 550), margin=5Plots.mm)

    return full_layout
end

function plot_reconstructed_fields_mean(
    points_flat, 
    mean_field, 
    boundary_positions_flat;         # Semicolon marks the start of keyword arguments
    source_positions_flat = nothing  # Defaults to nothing if omitted
)
    
    # 1. Extract X and Y coordinates using your reshaping function
    pos_matrix = flat_to_pos_matrix(points_flat)
    x_pts, y_pts = pos_matrix[1, :], pos_matrix[2, :]

    boundary_positions_matrix = flat_to_pos_matrix(boundary_positions_flat)
    x_boundary, y_boundary = boundary_positions_matrix[1, :], boundary_positions_matrix[2, :]

    # 2. Plot: Reconstructed Mean Field
    p1 = scatter(x_pts, y_pts, 
                 zcolor=mean_field, 
                 colorbar=true,
                 cmap=:balance,              
                 markersize=5,
                 markerstrokewidth=0,        
                 label="",
                 aspect_ratio=:equal,
                 xlabel="X", ylabel="Y",
                 title="Mean Field Intensity"
       )

    # Overlay the unit circle boundary line
    plot!(p1, x_boundary, y_boundary, color=:black, linestyle=:dash, linewidth=1.5, label="r = 1")

    # 3. Optional: Overlay Sources if provided
    if source_positions_flat !== nothing && !isempty(source_positions_flat)
        source_positions_matrix = flat_to_pos_matrix(source_positions_flat)
        x_sources, y_sources = source_positions_matrix[1, :], source_positions_matrix[2, :]
        
        scatter!(p1, x_sources, y_sources, 
                 label="Sources", 
                 color=:red, 
                 markersize=6, 
                 marker=:circle) 
    end

    full_layout = plot(p1, layout=(1, 1), size=(450, 450), margin=5Plots.mm)
    return full_layout
end

function plot_reconstructed_fields_variance(
    points_flat,  
    variance_field, 
    boundary_positions_flat;         # Semicolon marks the start of keyword arguments
    source_positions_flat = nothing  # Defaults to nothing if omitted
)
    
    # 1. Extract X and Y coordinates using your reshaping function
    pos_matrix = flat_to_pos_matrix(points_flat)
    x_pts, y_pts = pos_matrix[1, :], pos_matrix[2, :]

    boundary_positions_matrix = flat_to_pos_matrix(boundary_positions_flat)
    x_boundary, y_boundary = boundary_positions_matrix[1, :], boundary_positions_matrix[2, :]

    # 2. Plot: Reconstructed Variance Field
    p1 = scatter(x_pts, y_pts, 
                 zcolor=variance_field, 
                 colorbar=true,
                 cmap=:inferno,              
                 markersize=5,
                 markerstrokewidth=0,        
                 label="",
                 aspect_ratio=:equal,
                 xlabel="X", ylabel="Y",
                 title="Variance Field Intensity"
                )

    # Overlay the unit circle boundary line
    plot!(p1, x_boundary, y_boundary, color=:black, linestyle=:dash, linewidth=1.5, label="r = 1")

    # 3. Optional: Overlay Sources if provided
    if source_positions_flat !== nothing && !isempty(source_positions_flat)
        source_positions_matrix = flat_to_pos_matrix(source_positions_flat)
        x_sources, y_sources = source_positions_matrix[1, :], source_positions_matrix[2, :]
        
        scatter!(p1, x_sources, y_sources, 
                 label="Sources", 
                 color=:red, 
                 markersize=6, 
                 marker=:circle) 
    end

    full_layout = plot(p1, layout=(1, 1), size=(450, 450), margin=5Plots.mm)
    return full_layout
end