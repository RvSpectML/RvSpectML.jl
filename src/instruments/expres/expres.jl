module EXPRES
using RvSpectML
using DataFrames, FITSIO

#type EXPRES <: AbstractInstrument end
struct EXPRES1D <: AbstractInstrument1D end
struct EXPRES2D <: AbstractInstrument2D end
AnyEXPRES = Union{EXPRES1D,EXPRES2D}
export EXPRES1D, EXPRES2D, AnyEXPRES

#include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default

#include("io.jl")
export make_manifest, read_header, read_header, read_data, read_solar_data

end
