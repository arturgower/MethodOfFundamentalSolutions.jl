module MethodOfFundamentalSolutions

using Accessors
using LinearAlgebra
using BlockArrays: mortar
using StaticArrays: SVector, SMatrix
using MultipleScattering
using SpecialFunctions: hankelh1,besselj
using Distributions
using Optim


import MultipleScattering: PhysicalMedium, ScalarMedium, spatial_dimension, field_dimension, Shape, Box, bounding_box, points_in_shape, cartesian_to_radial_coordinates, radial_to_cartesian_transform, cartesian_to_radial_transform, field
export cartesian_to_radial_coordinates, radial_to_cartesian_transform, cartesian_to_radial_transform

import Statistics: mean

# for ploting recipes
using RecipesBase


export BoundaryData  # types
export compute_outward_normals, points_in_shape
include("boundarydata.jl")

export interior_points_along_coordinate, flat_to_pos_matrix
include("utils.jl")

export log_marginal_likelihood, optimise_hyperparameters, compute_coefficient_posterior, reconstruct_full_field, construct_prior, geometric_covariance # Bayesian functions
export ProbabilityDistribution, GaussianDistribution, BayesianSolver # Prior types
include("bayesian.jl")

export Simulation, ParticularSolution, TikhonovSolver, NoParticularSolution # types
export greens, source_positions, greens_gradient
export solve, system_matrix, system_matrix_gradient
include("solve.jl")

export FieldResult, FundamentalSolution
export field, field_covariance, field_std # types
include("results.jl")

export DisplacementType ,TractionType, Elastostatic, Acoustic, LaplaceMedium, DirichletType, NeumannType, laplace_M, laplace_grad_M # types and functions
export ParticularGravity # types
include("physics/elastic.jl")
include("physics/acoustic.jl")
include("physics/laplace.jl")
include("../plot/plot.jl")


end