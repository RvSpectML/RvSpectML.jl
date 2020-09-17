"""
    Code to compute CCF
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactored and optimized by Eric Ford
"""

""" Module for computing CCFs """
module CCF

import SpecialFunctions: erf

import ..RvSpectML
import ..RvSpectML: AbstractChuckOfSpectrum, AbstractChunkList, AbstractChunkListTimeseries
import ..RvSpectML: AbstractInstrument
import ..RvSpectML: num_chunks

# physical constants
const c_ms = 2.99782458e8    # m/s
#const c_kms = 2.99782458e5   # km/s

include("mask_shapes/mask_shapes.jl")
export AbstractCCFMaskShape, TopHatCCFMask

include("line_list.jl")
export AbstractLineList, BasicLineList

# default parameters
const default_v_center = 0.0    # m/s
const default_v_step = 250.0    # m/s
const default_v_max = 15.0e3    # m/s
const default_v_width = 410.0   # m/s
const default_v_range_no_mask_change = default_v_max # m/s

include("plan.jl")
export AbstractCCFPlan, BasicCCFPlan
export calc_ccf_v_grid

include("calc_ccf.jl")
export ccf_1D, ccf_1D!

#using ThreadTools
include("convenience.jl")
export calc_ccf_chunk, calc_ccf_chunklist, calc_ccf_chunklist_timeseries

# Will soon comment out to reduce precompilation time
include("calc_ccf_old.jl")
include("convenience_old.jl")
export ccf_1D_old, ccf_1D_old!

end # module
