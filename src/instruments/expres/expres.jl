"""
   Delegates loading functions & traits for the EXPRES spectrograph
   http://exoplanets.astro.yale.edu/expresBlog/
   https://ui.adsabs.harvard.edu/abs/2016SPIE.9908E..6TJ/abstract
Author: Eric Ford and collaborators
Created: August 2020
"""
module EXPRES
using RvSpectML
using DataFrames, FITSIO
using Dates  # If need to use datetime2julian() to get jd.  Need to check about getting BJD.

#type EXPRES <: AbstractInstrument end
struct EXPRES1D <: AbstractInstrument1D end
struct EXPRES2D <: AbstractInstrument2D end
const AnyEXPRES = Union{EXPRES1D,EXPRES2D}
export EXPRES1D, EXPRES2D, AnyEXPRES

#include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default

#include("io.jl")
export make_manifest, read_data, read_solar_data
# read_header not exported to avoid conflict with FITSIO.read_header

end
