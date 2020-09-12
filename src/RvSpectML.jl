"""
Author: Eric Ford and collaborators (see header of each file)
Created: August 2020
Contact: https://github.com/eford/RvSpectML.jl
"""

"""
[**RvSpectML.jl**](https://github.com/eford/RvSpectML.jl) is a Julia package for measuring radial velocities (RVs) from stellar spectroscopic timeseries via machine learning (ML).

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
export absorption_line
export speed_of_light_mps
#export searchsortednearest

include("util/spectra.jl")
export calc_normalization, normalize_spectrum!, get_λ_range
include("util/chunks.jl")
export make_chunk_list, make_orders_into_chunks, filter_bad_chunks, make_grid_for_chunk
#export find_cols..., find_orders..., findall_line,...

# alg.jl  is responsible for exporting its own types, modules & functions
include("alg/alg.jl")

# instruments.jl & the instruments it contains are responsible for exporting their own functions & modules
include("instruments/instruments.jl")

# util/files.jl
include("util/files.jl")
export make_manifest, code_to_include_param_jl

include("util/pipeline.jl")
using .Pipeline
export PipelinePlan
export make_plot, save_plot, save_data, need_to, need_to!, dont_need_to!, has_cache, reset_all_needs!

# util/plots.jl  is responsible for exporting its own functions
include("util/plots.jl")
export plot_spectrum_chunks
export plot_basis_vectors, plot_basis_scores
export add_time_gap_lines

end
