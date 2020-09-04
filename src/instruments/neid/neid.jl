"""
   Delegates loading functions & traits for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

module NEID
using ..RvSpectML
import ..RvSpectML: AbstractInstrument, AbstractInstrument1D, AbstractInstrument2D
using DataFrames, Query, FITSIO

#type NEID <: AbstractInstrument end
""" Trait for 1D Extracted spectra from NEID """
struct NEID1D <: AbstractInstrument1D end

""" Trait for 2D Extracted spectra from NEID """
struct NEID2D <: AbstractInstrument2D end

# Trait for any spectra from NEID (could improve by using SimpleTraits)
const AnyNEID = Union{NEID1D,NEID2D}
export NEID, NEID1D, NEID2D, AnyNEID

include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default
export default_ccf_mask_v_width

include("io.jl")
export make_manifest
# export make_manifest
export read_metadata, read_data, read_solar_data
# read_header not exported to avoid conflict with FITSIO.read_header
export read_drift_corrections!, read_barycentric_corrections!

end
