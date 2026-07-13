using MethodOfFundamentalSolutions
using Test, Statistics, LinearAlgebra, Distributions, StaticArrays
using SpecialFunctions
using Random

const MFS = MethodOfFundamentalSolutions

include("boundarydata_test.jl")
include("benchmarks_test.jl")
include("acoustics_test.jl")
include("plot_test.jl")
include("bayesian_test.jl")
include("variational_test.jl")
include("variational_geometry_test.jl")

