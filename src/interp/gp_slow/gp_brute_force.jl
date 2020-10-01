"""
Author: Eric Ford
Adapted from: https://github.com/eford/RvSpectraKitLearn.jl/blob/master/src/deriv_spectra_gp.jl

Uses brute-force, but could update GaussianProcesses.jl, Stheno.jl, TemporalGPs., KernelFunctions.jl, etc.
Note:  If make faster version using pacakge(s) that take advantage of specific kernels,
       please keep this around since it gets the job done with has minimal dependancies.
"""

"""
Module for interpolating via Gaussian Process Regression.
"""
module GPInterpolation
using PDMats
using LinearAlgebra
using RvSpectMLBase

export make_kernel_data, make_kernel_obs_pred, log_pdf_gp_posterior
#export predict_mean, predict_deriv, predict_deriv2, predict_mean_and_deriv, predict_mean_and_derivs

include("low_level.jl")
include("convenience.jl")
export interp_chunk_to_shifted_grid_gp_brute_force, interp_chunk_to_shifted_grid_gp_brute_force!
export interp_chunk_to_grid_gp_brute_force, interp_chunk_to_grid_gp_brute_force!


end
