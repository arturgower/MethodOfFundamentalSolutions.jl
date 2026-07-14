# recipe for plotting a BoundaryData. Assumes 2D boundary points. The boundary points and
# the fields may each be a plain vector of vectors, a single `MvNormal` over the flattened
# values, or a vector of `MvNormal`s (one per point); the mean is plotted, and when the
# boundary points carry a covariance it is drawn as position error bars.
@recipe function plot(cloud::BoundaryData{F, Dim}; fields = false, outward_normals = true, interior_points = true) where {F, Dim}
    # Set default attributes
    markershape --> :circle
    legend --> false
    aspect_ratio --> 1.0

    bs = struct_points(cloud.boundary_points, Dim)
    ns = cloud.outward_normals
    is = cloud.interior_points

    # used to determine the length to plot the vectors
    lengthscale = mean([norm(b - mean(bs)) for b in bs])

    bx = [b[1] for b in bs]
    by = [b[2] for b in bs]
    ix = [b[1] for b in is]
    iy = [b[2] for b in is]

    # position uncertainty: the per-point standard deviations from the boundary covariance
    Σx = cov(cloud.boundary_points)
    point_std = Σx isa UniformScaling ? nothing : sqrt.(max.(diag(Σx), 0.0))

    @series begin
        label --> "Boundary points"
        seriestype --> :scatter
        markersize --> 4.0
        if point_std !== nothing
            xerror --> [point_std[(i - 1) * Dim + 1] for i in eachindex(bs)]
            yerror --> [point_std[(i - 1) * Dim + 2] for i in eachindex(bs)]
        end
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
        nx = [b[1] for b in ns] .* (lengthscale / 4)
        ny = [b[2] for b in ns] .* (lengthscale / 4)
        # fallback: repeat or zero-length normals
        if length(nx) != length(bx)
            nx = zeros(length(bx))
            ny = zeros(length(by))
        end
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

    # quiver of the mean field, for vector-valued fields only (a scalar field has no
    # direction to draw)
    if fields
        fmean = flat_fields(cloud.fields)
        FD = length(fmean) ÷ length(bs)
        if FD >= 2
            fs = [fmean[(i - 1) * FD .+ (1:2)] for i in eachindex(bs)]
            f_lengthscale = mean(norm.(fs)) > 0 ? mean(norm.(fs)) : 1.0
            fx = [f[1] for f in fs] .* (lengthscale / 4) ./ (sqrt(2.0) * f_lengthscale)
            fy = [f[2] for f in fs] .* (lengthscale / 4) ./ (sqrt(2.0) * f_lengthscale)
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
end

@recipe function plot(sim::Simulation; source_positions = true, boundary_points = true)
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
    if boundary_points
        @series sim.boundary_data
    end
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

# Sample and plot the field of a fundamental solution over the interior of `bd`, then
# overlay the source points. The heatmap is produced by delegating to the FieldResult
# recipe below; see [`predict_field`](@ref) for the keyword arguments.
@recipe function plot(fsol::FundamentalSolution, bd::BoundaryData;
        xres = 30, yres = 30, normal_vec = nothing, field_transform = identity, sources = true)

    aspect_ratio --> 1.0
    legend --> false

    fr = predict_field(fsol, bd;
        xres = xres, yres = yres, normal_vec = normal_vec, field_transform = field_transform)

    @series begin
        seriestype := :heatmap
        region_shape := bd
        fr
    end

    if sources
        @series begin
            label --> "Source points"
            seriestype := :scatter
            markersize --> 4.0
            markercolor --> :red
            ([p[1] for p in fsol.positions], [p[2] for p in fsol.positions])
        end
    end
end

# Per-source prior standard deviation 1/√αᵢ learned by the variational solver, one value per
# source (averaged over the field/real-imaginary degrees of freedom of each source). A large
# value marks a relevant source, a small one a source ARD has switched off.
function _source_prior_std(vsol::VariationalSolution)
    pos = vsol.fsol.positions
    n_src = length(pos)
    α = vsol.prior_precisions
    per = max(1, length(α) ÷ n_src)
    src_std = [sqrt(mean(1 ./ α[(j - 1) * per + 1:j * per])) for j in 1:n_src]
    return [p[1] for p in pos], [p[2] for p in pos], src_std
end

# Scatter of the learned source positions, sized and coloured by their prior standard
# deviation 1/√αᵢ (the ARD relevance): the location and the variance of the sources.
@recipe function plot(vsol::VariationalSolution)
    aspect_ratio --> 1.0
    colorbar_title --> "prior std  1/√α"

    sx, sy, src_std = _source_prior_std(vsol)
    ms = 4.0 .+ 8.0 .* src_std ./ max(maximum(src_std), eps())

    @series begin
        seriestype := :scatter
        label --> "sources (size, colour ∝ prior std)"
        marker_z := src_std
        markersize := ms
        seriescolor --> :viridis
        (sx, sy)
    end
end

# The variational field over the interior of `bd`, with the learned sources overlaid and
# sized/coloured by their prior standard deviation.
@recipe function plot(vsol::VariationalSolution, bd::BoundaryData;
        xres = 30, yres = 30, normal_vec = nothing, field_transform = identity)

    aspect_ratio --> 1.0
    legend --> false

    fr = predict_field(vsol.fsol, bd;
        xres = xres, yres = yres, normal_vec = normal_vec, field_transform = field_transform)

    @series begin
        seriestype := :heatmap
        region_shape := bd
        fr
    end

    sx, sy, src_std = _source_prior_std(vsol)
    ms = 4.0 .+ 8.0 .* src_std ./ max(maximum(src_std), eps())
    @series begin
        seriestype := :scatter
        label --> "sources (size ∝ prior std)"
        markersize := ms
        markercolor --> :red
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
