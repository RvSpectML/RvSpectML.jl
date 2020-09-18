""" Delegates loading of code with modules, types, functions and parameters for various algorithms """

include("interp/interp.jl")
export pack_chunk_list_timeseries_to_matrix
export pack_shifted_chunk_list_timeseries_to_matrix

include("line_finder/line_finder.jl")
export LineFinder, LineFinderPlan
export find_lines_in_chunk, find_lines_in_chunklist, find_lines_in_chunklist_timeseries

include("combine_obs.jl")
export make_template_spectra
export bin_spectra_consecutive, bin_spectra_max_Δt, bin_spectra_nightly
export bin_times_consecutive, bin_rvs_consecutive
export bin_times_max_Δt, bin_rvs_max_Δt, bin_times_and_rvs_max_Δt
export bin_times_nightly, bin_rvs_nightly, bin_times_and_rvs_nightly
export rms_rvs_within_night


include("ccf/ccf.jl")
export CCF
export calc_ccf_Δv_grid
export calc_ccf_chunk, calc_ccf_chunklist, calc_ccf_chunklist_timeseries

include("calc_rv/ccf_std.jl")
using .RVFromCCF
export measure_rv_from_ccf, measure_rvs_from_ccf
export AbstractMeasureRvFromCCF, MeasureRvFromCCFCentroid, MeasureRvFromCCFQuadratic, MeasureRvFromCCFGaussian, MeasureRvFromCCFBestFit

include("calc_rv/project_flux.jl")
export RVFromCCF
export calc_dfluxdlnlambda, calc_d2fluxdlnlambda2
export calc_mean_spectrum, calc_mean_dfluxdlnlambda, calc_mean_d2fluxdlnlambda2
export calc_rvs_from_taylor_expansion, calc_chunk_rvs_from_taylor_expansion

include("scalpels.jl")
export Scalpels
export clean_rvs_scalpels, scalpels_rms_vs_num_basis

include("dcpca.jl")
export DCPCA
export doppler_constrained_pca
export compute_spectra_perp_doppler_shift

#include("ppca.jl")       # TODO:  Need to fold into this package, adjust types, and add test/example
