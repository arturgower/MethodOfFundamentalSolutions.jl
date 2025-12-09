module MethodOfFundamentalSolutions

using Accessors
using LinearAlgebra
using BlockArrays: mortar
using StaticArrays: SVector
using MultipleScattering
using SpecialFunctions: hankelh1,besselj

import MultipleScattering: PhysicalMedium, ScalarMedium, spatial_dimension, field_dimension, Shape, Box, bounding_box, points_in_shape, cartesian_to_radial_coordinates, radial_to_cartesian_transform, cartesian_to_radial_transform, field
export cartesian_to_radial_coordinates, radial_to_cartesian_transform, cartesian_to_radial_transform

import Statistics: mean

# for ploting recipes
using RecipesBase


export BoundaryData  # types
export outward_normals, points_in_shape
include("boundarydata.jl")

export interior_points_along_coordinate
include("utils.jl")

export Simulation, ParticularSolution, TikhonovSolver # types
export greens, source_positions
export solve, system_matrix
include("solve.jl")

export FieldResult, FundamentalSolution
export field # types
include("results.jl")

export DisplacementType, TractionType, Elastostatic, Acoustic
export ParticularGravity # types
include("physics/elastic.jl")
include("physics/acoustic.jl")
include("../plot/plot.jl")

end