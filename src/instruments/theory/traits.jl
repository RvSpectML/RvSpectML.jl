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

import ..RvSpectML: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyTheoreticalInstrument) = 500.0  #

import ..RvSpectML: get_inst_module
get_inst_module(::AnyTheoreticalInstrument) = TheoreticalInstrument

import ..RvSpectML: get_λ_range
function get_λ_range(data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyTheoreticalInstrument, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = extrema(data.λ)
   return (min=λmin, max=λmax)
end

import ..RvSpectML: filter_line_list, find_worst_telluric_in_each_chunk

export filter_line_list, find_worst_telluric_in_each_chunk
export get_inst_module
