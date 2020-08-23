abstract type AbstractInstrument end

abstract type AbstractInstrument1D <: AbstractInstrument end
abstract type AbstractInstrument2D <: AbstractInstrument end

#=
For each instrument/data file type, need to create a sub-type of either AbstractInstrument2D or AbstractInstrument1D and
to implement following functions for each instrument:
metadata_symbols_default(::NEID) = Symbol[:bjd, :target, :ssbz]
metadata_strings_default(::NEID) = String["OBSJD", "SKY-OBJ", "SSBZ000"]
and
  If 2d:
    min_order(::NEID2D) = 1
    max_order(::NEID2D) = 90
    min_pixel_in_order(::NEID2D) = 1
    max_pixel_in_order(::NEID2D) = 9216
    orders_to_use_default(::NEID2D) = neid_orders_all
    min_col_default(::NEID2D) = 451
    max_col_default(::NEID2D) = 9216 #- min_col_default(NEID2D())
  or if 1d:
    min_pixel(::NEID1D) = 1
    max_pixel(::NEID1D) = 90*9216  # TODO: Update once know size of NEID's 1d extracted spectra
    min_col_default(::NEID1D) = 451
    max_col_default(::NEID1D) = 9216 - min_col_default(NEID1D())
=#

export min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_col_default, max_col_default
export orders_all, pixels_all, max_pixels_in_spectra       # generic implementations avaliable
export metadata_symbols_default, metadata_strings_default  # need to specialize

include("common.jl")
export read_manifest, read_metadata_from_fits

include("neid/neid.jl")
export NEID1D, NEID2D, AnyNEID

#include("expres/expres.jl")
#include("neid/harps-n.jl")
