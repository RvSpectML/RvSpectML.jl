"""
Author: Eric Ford and collaborators (see header of each file)
Created: August 2020
Contact: https://github.com/eford/RvSpectML.jl
"""

"""
[**RvSpectML.jl**](https://github.com/eford/RvSpectML.jl) is a Julia package for measuring radial velocities (RVs) from stellar spectroscopic timeseries via machine learning (ML).

"""
module RvSpectML

# Packages from RvSpectML ecosystem
using RvSpectMLBase
export speed_of_light_mps
export calc_normalization, normalize_spectrum!, get_λ_range
export make_chunk_list, make_orders_into_chunks
export make_chunk_list_timeseries, make_order_list_timeseries
export make_grid_for_chunk, filter_bad_chunks

using EchelleInstruments
export EchelleInstruments

using EchelleCCFs
export EchelleCCFs
CCFs = EchelleCCFs

using Scalpels
export Scalpels

#using Experimental
#using RvSpectMLPlots

# Packages we know we'll use in many places
using LinearAlgebra, Statistics
using DataFrames, Query
using Dates
using AstroLib  # For solar position
# Packages that are being used and likely can be shared
# using Distributions, Interpolations, MultivariateStats, PDMats
# Packages we might use soon
#using Optim
#using Stheno, TemporalGPs
# Let's try to keep Plots out the main package to reduce compilation times when not plotting.
#using  Plots

include("pipeline/pipeline.jl")
using .Pipeline
export Pipeline
export extract_orders, prepare_line_list, ccf_total, ccf_orders, calc_rvs_from_ccf_total

#=
export PipelinePlan
export make_plot,  save_plot, save_data, need_to, has_cache, read_cache, set_cache!  # Query pipeline
export need_to!, dont_need_to!, reset_all_needs!             # Write to pipeline
export make_plot!, dont_make_plot!, make_all_plots!,  make_no_plots!
=#

include("interp/interp.jl")
#=
export TemporalGPInterpolation
export construct_gp_posterior, gp_marginal, predict_gp
export predict_mean, predict_deriv, predict_deriv2, predict_mean_and_deriv, predict_mean_and_derivs
=#


include("line_finder/line_finder.jl")
export LineFinderPlan

include("line_shapes/get_line_shapes.jl")

# alg.jl  will be responsible for exporting its own types, modules & functions
#include("alg/alg.jl")
# For now treating parts individually
include("alg/project_flux_common_wavelengths.jl")
include("alg/make_template_spectrum.jl")
#export interp_chunk_to_shifted_grid_gp_temporal
include("alg/combine_obs.jl")
export bin_times_consecutive, bin_rvs_consecutive, bin_spectra_consecutive
export bin_times_nightly, bin_rvs_nightly, bin_times_and_rvs_nightly, bin_spectra_nightly, rms_rvs_within_night
export bin_times_max_Δt, bin_rvs_max_Δt,  bin_times_and_rvs_max_Δt, bin_spectra_max_Δt

include("alg/dcpca.jl")
#export DCPCA
using .DCPCA
export clean_rvs_dcpca, calc_sigma_pca_scores

#include("alg/ppcca.jl")
#include("alg/rvs_from_gp_pairs.jl")

include("util/spectra.jl")
export calc_depth_and_expected_rv_precission, calc_formal_rv_precission
include("util/sun.jl")
export calc_solar_alt

end
