"""
Author: Eric Ford
Adapted from: https://github.com/eford/RvSpectraKitLearn.jl/blob/master/src/deriv_spectra_gp.jl

GPs using Stheno.jl, TemporalGPs., KernelFunctions.jl, etc.
"""

"""
Module for interpolating via Gaussian Process Regression based on Stheno and TemporalGPs packages.
"""
module TemporalGPInterpolation
using LinearAlgebra
using PDMats
using StaticArrays
using Stheno, TemporalGPs
using Dates
using RvSpectMLBase

#export make_kernel_data, make_kernel_obs_pred,
export gp_marginal
export predict_mean, predict_deriv, predict_deriv2, predict_mean_and_deriv, predict_mean_and_derivs

import Stheno: AbstractGP
import Distributions: AbstractMvNormal

using Distributions
const Fx_PosteriorType = Distribution{Multivariate,Continuous}

include("low_level.jl")
export  construct_gp_prior, construct_gp_posterior, log_pdf_gp_posterior

include("convenience.jl")
export interp_chunk_to_shifted_grid_gp_temporal, interp_chunk_to_shifted_grid_gp_temporal!
export interp_chunk_to_grid_gp_temporal, interp_chunk_to_grid_gp_temporal!


end # module
