"""
   Delegates loading functions & traits for the HARPS-N spectrograph
Author: Eric Ford
Created: August 2020
"""
module HARPSN
using RvSpectML
using DataFrames, FITSIO

#type EXPRES <: AbstractInstrument end
struct HARPSN1D <: AbstractInstrument1D end
struct HARPSN2D <: AbstractInstrument2D end
const AnyHARPSN = Union{HARPSN1D,HARPSN2D}
export HARPSN1D, HARPSN2D, AnyHARPSN

#include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default

#include("io.jl")
export make_manifest, read_data, read_solar_data
# read_header not exported to avoid conflict with FITSIO.read_header

end
