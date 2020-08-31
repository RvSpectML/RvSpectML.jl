"""
   Traits for a theoretical spectrograph
Author: Eric Ford
Created: August 2020
"""

""" Delegates loading of code specifying types essential to the package.  """

import ..RvSpectML: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel
min_order(::TheoreticalInstrument2D) = 1
max_order(::TheoreticalInstrument2D) = 128
min_pixel_in_order(inst::TheoreticalInstrument2D) = 1
max_pixel_in_order(inst::TheoreticalInstrument2D) = 8192

min_pixel(inst::TheoreticalInstrument1D) = 1
max_pixel(inst::TheoreticalInstrument1D) = 128*8192

import ..RvSpectML: orders_to_use_default, min_col_default, max_col_default
orders_to_use_default(inst::TheoreticalInstrument2D) = min_order(inst):max_order(inst)
min_col_default(::TheoreticalInstrument2D, ord::Integer) = 1
max_col_default(::TheoreticalInstrument2D, ord::Integer) = max_pixel_in_order(inst)

#import ..RvSpectML: metadata_symbols_default, metadata_strings_default
#metadata_symbols_default(::AnyD) = Symbol[:bjd, :target, :ssbz]
#metadata_strings_default(::AnyD) = String["OBSJD", "SKY-OBJ", "SSBZ000"]

import ..RvSpectML: default_ccf_v_width
default_ccf_v_width(::AnyTheoreticalInstrument) = 620.953  # from NEID
