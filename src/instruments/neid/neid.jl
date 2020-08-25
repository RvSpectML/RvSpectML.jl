"""
   Functions & traits for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

module NEID
using RvSpectML
using DataFrames, FITSIO

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

include("io.jl")
export make_manifest, read_metadata, read_header, read_data, read_solar_data
export read_drift_corrections!, read_barycentric_corrections!

end
