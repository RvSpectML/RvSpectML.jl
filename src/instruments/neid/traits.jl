import RvSpectML: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel
#import RvSpectML.min_order, RvSpectML.max_order, RvSpectML.min_pixel_in_order, RvSpectML.max_pixel_in_order, RvSpectML.min_pixel, RvSpectML.max_pixel
min_order(::NEID2D) = 1
max_order(::NEID2D) = 90
min_pixel_in_order(inst::NEID2D) = 1
max_pixel_in_order(inst::NEID2D) = 9216

min_pixel(::NEID1D) = 1
max_pixel(::NEID1D) = 90*9216 # TODO: Update once know size of NEID's 1d extracted spectra

import RvSpectML: orders_to_use_default, min_col_default, max_col_default
#orders_to_use_default(inst::NEID2D) = 1:52   # Avoiding redder orders due to tellurics
orders_to_use_default(inst::NEID2D) = 1:71   # Avoiding 72 because of NaNs in solar data
min_col_default(::NEID2D) = 451              # Avoiding smaller columns due to NaNs
#min_col_default(::NEID2D) = 2000              # Avoiding smaller columns due to lower flux and distortions
max_col_default(::NEID2D) = 9216 - (min_col_default(NEID2D())-1)   # Avoiding larger columns for symmetry

import RvSpectML: metadata_symbols_default, metadata_strings_default
metadata_symbols_default(::AnyNEID) = Symbol[:bjd, :target, :ssbz]
metadata_strings_default(::AnyNEID) = String["OBSJD", "SKY-OBJ", "SSBZ000"]
