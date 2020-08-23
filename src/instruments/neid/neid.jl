module NEID
using RvSpectML
using DataFrames, FITSIO

#type NEID <: AbstractInstrument end
struct NEID1D <: AbstractInstrument1D end
struct NEID2D <: AbstractInstrument2D end
AnyNEID = Union{NEID1D,NEID2D}

min_order(::NEID2D) = 1
max_order(::NEID2D) = 90
min_pixel_in_order(::NEID2D) = 1
max_pixel_in_order(::NEID2D) = 9216

min_pixel(::NEID1D) = 1
max_pixel(::NEID1D) = 90*9216 # TODO: Update once know size of NEID's 1d extracted spectra

orders_to_use_default(::NEID2D) = 1:60
min_col_default(::NEID2D) = 451
max_col_default(::NEID2D) = 9216 - (min_col_default(NEID2D())-1)

metadata_symbols_default(::AnyNEID) = Symbol[:bjd, :target, :ssbz]
metadata_strings_default(::AnyNEID) = String["OBSJD", "SKY-OBJ", "SSBZ000"]

export NEID, NEID1D, NEID2D, AnyNEID
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default

include("io.jl")
export make_manifest, read_header, read_header, read_data, read_solar_data

end
