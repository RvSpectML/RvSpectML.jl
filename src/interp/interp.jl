"""
Delegates loading of code for different interpolation aglorithms (in their own modules).
Also wraps these modules, so they're called easily by `pack_chunk_list_timeseries_to_matrix`.

Author: Eric Ford
Created: August 2020
"""

# First include modules for different interpolation algoritihms
# File called module will include whatever packages it depends on, and also
# provide at least one function with a common interface (currently named) interp_chunk_to_grid_ALGNAME!.

include("common.jl")
export numerical_deriv

include("linear/linear.jl")
import .LinearInterpolation
export interp_chunk_to_shifted_grid_linear, interp_chunk_to_shifted_grid_linear!
export interp_chunk_to_grid_linear, interp_chunk_to_grid_linear!

include("sinc/sinc.jl")
using .SincInterpolation
export interp_chunk_to_shifted_grid_sinc, interp_chunk_to_shifted_grid_sinc!
export interp_chunk_to_grid_sinc, interp_chunk_to_grid_sinc!

#=
include("gp_slow/gp_brute_force.jl")             # Basic, brute-force GPs with no dependancies, hasn't been updated recently
import .GPInterpolation
export interp_chunk_to_shifted_grid_gp_brute_force, interp_chunk_to_shifted_grid_gp_brute_force!
export interp_chunk_to_grid_gp_brute_force, interp_chunk_to_grid_gp_brute_force!
=#

include("gp_fast/temporalgps.jl")    # Dependenacies on Stheno, Temporal GPs and Static Arrays
import .TemporalGPInterpolation
export interp_chunk_to_shifted_grid_gp_temporal, interp_chunk_to_shifted_grid_gp_temporal!
export interp_chunk_to_grid_gp_temporal, interp_chunk_to_grid_gp_temporal!
