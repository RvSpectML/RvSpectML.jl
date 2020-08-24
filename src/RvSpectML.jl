"""
[]**RvSpectML.jl**](https://github.com/eford/RvSpectML.jl) is a Julia package for measuring radial velocities (RVs) from stellar spectroscopic timeseries via machine learning (ML).

Author: Eric Ford and collaborators
Created: August 2020
Contact: https://github.com/eford/RvSpectML.jl
"""

module RvSpectML

using LinearAlgebra, Statistics
using DataFrames, Query
using Distributions, Interpolations, MultivariateStats, PDMats
#using Optim, Plots

include("types/types.jl")
include("util/util.jl")
export calc_doppler_factor
include("instruments/instruments.jl")
include("alg/alg.jl")

end
