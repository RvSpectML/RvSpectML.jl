"""
   Module for performing Sinc interpolation
From
https://github.com/christiangil/GP-Research/blob/master/julia/src/interpolation_functions.jl
https://github.com/christiangil/GP-Research/blob/master/julia/src/general_functions.jl

Authors: Joe Ninan (original)
         Christian Gilbertson (Converted to Julia and optimized)
         Eric Ford (Further adapted/optimizations)
"""
module SincInterpolation

using Base.Math
using DSP
using Interpolations
using RvSpectMLBase
#import ..RvSpectML: searchsortednearest

export spectra_interpolate

include("low_level.jl")

include("convenience.jl")
export interp_chunk_to_shifted_grid_sinc, interp_chunk_to_shifted_grid_sinc!
export interp_chunk_to_grid_sinc, interp_chunk_to_grid_sinc!

end # module
