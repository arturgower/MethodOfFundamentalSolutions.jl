function interior_points_along_coordinate(points; offset_percent = 0.1,
    coordinate_number = 2, segments = 10)

    if offset_percent < 0.0 || offset_percent > 1.0
        throw(ArgumentError("offset_percent must be between 0.0 and 1.0"))
    end
    
    cs = [xy[coordinate_number] for xy in points]

    minc = minimum(cs)
    maxc = maximum(cs)

    edges = collect(range(minc, stop=maxc, length=segments+1))

    # collect indices for each y-bin
    bin_indices = [Int[] for _ in 1:segments]
    
    for (idx, y) in enumerate(cs)
        b = clamp(searchsortedlast(edges, y), 1, segments)
        push!(bin_indices[b], idx)
    end

    bin_indices = filter(v -> !isempty(v), bin_indices)

    interior_points = map(bin_indices) do inds
        p = mean(points[inds])
        
        inner_low = minc + offset_percent * (maxc - minc)
        inner_high = maxc - offset_percent * (maxc - minc)
        pc = p[coordinate_number]
        t = clamp((pc - minc) / (maxc - minc), 0.0, 1.0)

        p = @set p[coordinate_number] = inner_low + t * (inner_high - inner_low)
        p
    end

    return interior_points
end


# function predict_field(fsol::FundamentalSolution, bd::BoundaryData; normal_vec = [0.0,1.0], yres = 35, xres = 15)

#     x_vec, inds = points_in_shape(bd; yres = yres, xres = xres)
#     xs = x_vec[inds]

#     fs = [
#         field(TractionType(), fsol, x, normal_vec) 
#     for x in xs];

#     field_mat = [[0.0, 0.0] for x in x_vec]
#     field_mat[inds] = fs;
#     field_predict = FieldResult(x_vec, field_mat[:]);

#     return field_predict
# end

# f1 = predict_field(fsol1, bd1; normal_vec = [0.0,1.0], yres = yres, xres = xres);
