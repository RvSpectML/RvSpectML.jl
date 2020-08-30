"""
   Delegates loading functions & traits for the EXPRES spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

min_order(::EXPRES2D) = 1
max_order(::EXPRES2D) = 86
min_pixel_in_order(::EXPRES2D) = 1
max_pixel_in_order(::EXPRES2D) = 7920

min_pixel(::EXPRES1D) = 1
max_pixel(::EXPRES1D) = orders_to_use_default(EXPRES2D())*(max_col_default(EXPRES2D())-min_col_default(EXPRES2D()+1)

orders_to_use_default(::EXPRES2D) = 42:74
min_col_default(::EXPRES2D) = 770
max_col_default(::EXPRES2D) = 6650

metadata_symbols_default(::AnyEXPRES) = Symbol[:midpoint] #, :target, :ssbz]
metadata_strings_default(::AnyEXPRES) = String["MIDPOINT"] #, "SKY-OBJ", "SSBZ000"]
