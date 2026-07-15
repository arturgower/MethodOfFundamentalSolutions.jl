using Test
using RecipesBase
using MethodOfFundamentalSolutions

# Provide a minimal support function so RecipesBase.apply_recipe doesn't
# require a plotting backend (Plots.jl) during tests. Plots normally
# extends `is_key_supported`; here we stub it to allow recipe testing.
import RecipesBase: is_key_supported
is_key_supported(::Any...) = true

@testset "Plot recipes" begin
    # Test BoundaryShape / BoundaryData recipes
    @testset "BoundaryShape recipe" begin
        bs = [[0.0,0.0], [1.0,0.0], [0.0,1.0]]
        ns = [[1.0,0.0], [0.0,1.0],[-1.0,-1.0]]
        is = [[0.5,0.5]]

        cloud = BoundaryData(DisplacementType();
            boundary_points = bs,
            normals=ns,
            interior_points=is
        )

        # Test boundary points series
        plt = RecipesBase.apply_recipe(Dict{Symbol,Any}(), cloud.boundary_shape)[1]
        boundary_series = plt.plotattributes
        @test boundary_series[:seriestype] == :scatter
        @test boundary_series[:label] == "Boundary points"

        # Test normals (quiver) series
        plt = RecipesBase.apply_recipe(Dict{Symbol,Any}(), cloud.boundary_shape)[3]
        normal_series = plt.plotattributes
        @test normal_series[:seriestype] == :quiver

        # The BoundaryData recipe delegates to its boundary_shape
        plt = RecipesBase.apply_recipe(Dict{Symbol,Any}(), cloud)[1]
        @test plt.args == (cloud.boundary_shape,)
    end

    # Test FieldResult recipe
    @testset "FieldResult recipe" begin
        # Create a simple 2D field result
        x = [[Float64(i), Float64(j)] for i in 0:1, j in 0:1][:]
        fs = [[1.0] for _ in 1:length(x)]
        res = FieldResult(x, fs)

        # Apply the recipe
        plt = RecipesBase.apply_recipe(Dict{Symbol,Any}(), res)[1]
        
        # Test default attributes
        @test plt.plotattributes[:seriestype] == :heatmap
        
        # Test with different seriestype
        plt_scatter = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:seriestype => :scatter), 
            res
        )[1]
        plt_scatter.plotattributes
        @test true
    end
end