"""
   Delegates loading functions & traits for the EXPRES spectrograph
   http://exoplanets.astro.yale.edu/expresBlog/
   https://ui.adsabs.harvard.edu/abs/2016SPIE.9908E..6TJ/abstract
Author: Eric Ford and collaborators
Created: August 2020
"""

import ..RvSpectML: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel
min_order(::EXPRES2D) = 1
max_order(::EXPRES2D) = 86
min_pixel_in_order(::EXPRES2D) = 1
max_pixel_in_order(::EXPRES2D) = 7920

min_pixel(::EXPRES1D) = 1
#max_pixel(::EXPRES1D) = orders_to_use_default(EXPRES2D())*(max_col_default(EXPRES2D())-min_col_default(EXPRES2D()+1))

import ..RvSpectML: orders_to_use_default, min_col_default, max_col_default

# Values hard coded based on HD 101501.  Someone will need to make automatic or ask if these will stay fixed if these are to be used for anything other than exploring.
const minmax_col_nonnan = [2149:5973, 1920:6121, 1787:6333, 1634:6532, 1485:6744, 1327:6934, 1179:7096, 1061:7242, 956:7282, 853:7292, 777:7291, 711:7292, 678:7292, 639:7292, 581:7292, 547:7292, 518:7296, 456:7297, 382:7296, 365:7293, 356:7289, 348:7289, 336:7292, 318:7295, 295:7297, 272:7296, 258:7296, 247:7297, 238:7295, 229:7294, 222:7292, 216:7291, 211:7288, 206:7285, 202:7282, 199:7279, 197:7278, 193:7275, 189:7275, 183:7273, 180:7273, 173:7273, 167:7274, 162:7274, 157:7273, 154:7272, 150:7271, 146:7269, 143:7268, 140:7267, 137:7266, 133:7265, 130:7264, 127:7262, 124:7261, 121:7260, 118:7258, 115:7258, 112:7256, 110:7256, 107:7254, 104:7253, 101:7252, 99:7250, 96:7248, 94:7247, 91:7245, 89:7243, 87:7241, 85:7239, 83:7237, 81:7234, 79:7230, 78:7227, 76:7223, 75:7218, 75:7211, 75:7204, 76:7194, 79:7181, 84:7167, 93:7142, 109:7097, 157:6352, 500:5119, 3961:3961 ]
min_col_nonnan(::EXPRES2D, ord::Integer) = first(minmax_col_nonnan[ord])
max_col_nonnan(::EXPRES2D, ord::Integer) = last(minmax_col_nonnan[ord])

const minmax_col_excalibur_avail = vcat(fill(1:1,31),
    [ 2658:5073, 1962:5830, 1650:6030, 1550:6364, 1322:6405, 1156:6476, 1063:6848, 843:6927, 779:7053, 672:7136, 469:7172, 487:7177, 399:7148, 331:7203, 337:7191, 257:7196, 264:7201, 226:7207, 221:7212, 228:7218, 189:7206, 184:7212, 168:7200, 163:7206, 158:7212, 153:7200, 148:7206, 143:7194, 150:7200, 145:7207, 140:7194, 135:7201, 130:7188, 138:7195, 133:7182, 128:7189, 123:7176, 132:7183, 127:7170, 122:7177, 131:7164, 126:7150, 121:7113, 145:7098, 155:5401 ],
    fill(1:1, max_order(EXPRES2D())-77+1) )

min_col_excalibur(::EXPRES2D, ord::Integer) = first(minmax_col_excalibur_avail[ord])
max_col_excalibur(::EXPRES2D, ord::Integer) = last(minmax_col_excalibur_avail[ord])

orders_to_use_default(::EXPRES2D) = 43:75    # Based on methods for inital RV described at http://exoplanets.astro.yale.edu/science/activity.php
min_col_default(inst::EXPRES2D, ord::Integer) = 770
max_col_default(inst::EXPRES2D, ord::Integer) = 6650

import ..RvSpectML: metadata_symbols_default, metadata_strings_default
metadata_symbols_default(::AnyEXPRES) = Symbol[:midpoint, :target, :exposure_time, :airmass, :moondist, :sundist]
metadata_strings_default(::AnyEXPRES) = String["MIDPOINT", "OBJECT", "AEXPTIME", "AIRMASS", "MOONDIST", "SUNDIST"]

#metadata_hdu2_symbols_default(::AnyEXPRES) = Symbol[:S_indicator, :Hα_indicator, :Hα_width, :ccf_width, :σ_ccf_width, :bis_indicator]
#metadata_hdu2_strings_default(::AnyEXPRES) = String["S-VALUE", "HALPHA", "HWIDTH", "CCFFWHM", "CCFFWHME", "BIS"]
# Left out S and BIS since their type is to be standardized in next update of input files
metadata_hdu2_symbols_default(::AnyEXPRES) = Symbol[:Hα_indicator, :Hα_width, :ccf_width, :σ_ccf_width]
metadata_hdu2_strings_default(::AnyEXPRES) = String["HALPHA", "HWIDTH", "CCFFWHM", "CCFFWHME"]

import ..RvSpectML: default_ccf_v_width
default_ccf_v_width(::AnyEXPRES) = 448.0   # TODO: Update value for EXPRES
