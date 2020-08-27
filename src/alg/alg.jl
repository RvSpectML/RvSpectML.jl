include("interp.jl")
export make_interpolator_linear_flux, make_interpolator_linear_var
export make_interpolator_gp, interp_to_grid
export interp_chunk_to_grid_linear!, pack_chunk_list_timeseries_to_matrix
export calc_mean_spectrum, calc_deriv, calc_mean_deriv
export calc_rvs_from_taylor_expansion, calc_chunk_rvs_from_taylor_expansion
export compute_spectra_perp_doppler_shift

include("gp.jl")     # TODO:  Need to fold into this package, adjust types, and add test/example


#incldue("ccf.jl")    # TODO:  Need to fold into this package, adjust types, and add test/example
#include("ppca.jl")   # TODO:  Need to fold into this package, adjust types, and add test/example
