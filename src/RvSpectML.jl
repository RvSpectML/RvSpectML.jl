"""
[]**RvSpectML.jl**](https://github.com/eford/RvSpectML.jl) is a Julia package for measuring radial velocities (RVs) from stellar spectroscopic timeseries via machine learning (ML).

Author: Eric Ford and collaborators
Created: August 2020
Contact: https://github.com/eford/RvSpectML.jl
"""

module RvSpectML

# Packages we know we'll use in many places
using LinearAlgebra, Statistics
using DataFrames, Query
# Packages thare are being used and likely can be shared
using Distributions, Interpolations, MultivariateStats, PDMats
# Packages we might use soon
#using Optim, Stheno
# Let's try to keep Plots out the main package to reduce compilation times when not plotting.
#using  Plots

# types.jl is responsible for exporting its own types and functions
include("types/types.jl")

include("util/util.jl")
 export calc_doppler_factor, apply_doppler_boost!
 #export predict_line_width
 export λ_vac_to_air, λ_air_to_vac
 #export find_cols..., find_orders..., findall_line,...
 export make_chunk_list, make_orders_into_chunks, filter_bad_chunks, make_grid_for_chunk
 export merge_lines
 export calc_normalization, normalize_spectrum!

# alg.jl  is responsible for exporting its own types, modules & functions
include("alg/alg.jl")

# instruments.jl & the ihstruments it contains are responsible for exporting their own functions & modules
include("instruments/instruments.jl")

# util/plots.jl  is responsible for exporting its own functions
include("util/plots.jl")
export plot_spectrum_chunks
export plot_basis_vectors, plot_basis_scores

end
