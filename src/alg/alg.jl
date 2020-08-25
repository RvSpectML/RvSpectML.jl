include("interp.jl")
export make_interpolator_linear_flux, make_interpolator_linear_var
export make_interpolator_gp, interp_to_grid, pack_chunks_into_matrix
export calc_deriv, calc_mean_spectrum, calc_mean_deriv

#incldue("ccf.jl")    # TODO:  Need to fold into this package, adjust types, and add test/example
#include("ppca.jl")   # TODO:  Need to fold into this package, adjust types, and add test/example
