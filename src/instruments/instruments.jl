""" Subdirectories of src/instruments include provide functions specialized for each instrument,
    typically the file I/O and pre-processing, so data ends up in a common format.
    src/instruments/common.jl provides routines that can be shared by instruments. """


include("common.jl")
export read_manifest, read_header, read_metadata_from_fits
export read_mask_espresso, read_mask_vald

include("neid/neid.jl")
export NEID, NEID1D, NEID2D, AnyNEID

include("expres/expres.jl")
export EXPRES, EXPRES1D, EXPRES2D, AnyEXPRES

include("harps-n/harps-n.jl")
export HARPSN, HARPSN1D, HARPSN2D, AnyHARPSN

export min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_col_default, max_col_default
export orders_all, pixels_all, max_pixels_in_spectra       # generic implementations avaliable
export metadata_symbols_default, metadata_strings_default  # need to specialize

#=
For each instrument/data file type, need to create a sub-type of either AbstractInstrument2D or AbstractInstrument1D and
to implement the functions in neid/traits.jl for each instrument trait.
=#
