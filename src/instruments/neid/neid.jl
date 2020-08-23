module NEID
using RvSpectML
using DataFrames, FITSIO

#type NEID <: AbstractInstrument end
struct NEID1D <: AbstractInstrument1D end
struct NEID2D <: AbstractInstrument2D end
AnyNEID = Union{NEID1D,NEID2D}
export NEID, NEID1D, NEID2D, AnyNEID

include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default

include("io.jl")
export make_manifest, read_header, read_header, read_data, read_solar_data

end
