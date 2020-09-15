"""
   Traits for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

""" Delegates loading of code specifying types essential to the package.  """

import ..RvSpectML: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel
#import RvSpectML.min_order, RvSpectML.max_order, RvSpectML.min_pixel_in_order, RvSpectML.max_pixel_in_order, RvSpectML.min_pixel, RvSpectML.max_pixel
min_order(::NEID2D) = 1
max_order(::NEID2D) = 90
min_pixel_in_order(inst::NEID2D) = 1
max_pixel_in_order(inst::NEID2D) = 9216

min_pixel(::NEID1D) = 1
max_pixel(::NEID1D) = 90*9216 # TODO: Update once know size of NEID's 1d extracted spectra

import ..RvSpectML: orders_to_use_default, min_col_default, max_col_default
#orders_to_use_default(inst::NEID2D) = 1:52   # Avoiding redder orders due to tellurics
orders_to_use_default(inst::NEID2D) = 1:71   # Avoiding 72 because of NaNs in solar data
min_col_default(::NEID2D, ord::Integer) = 451              # Avoiding smaller columns due to NaNs
#min_col_default(::NEID2D) = 2000              # Avoiding smaller columns due to lower flux and distortions
max_col_default(::NEID2D, ord::Integer) = 9216 - (min_col_default(NEID2D(),ord)-1)   # Avoiding larger columns for symmetry

import ..RvSpectML: metadata_symbols_default, metadata_strings_default
metadata_symbols_default(::AnyNEID) = Symbol[:bjd, :target, :ssbz]
metadata_strings_default(::AnyNEID) = String["OBSJD", "SKY-OBJ", "SSBZ000"]

import ..RvSpectML: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyNEID) = 620.953

import ..RvSpectML: get_inst_module
get_inst_module(::AnyNEID) = NEID

import ..RvSpectML: get_λ_range
function get_λ_range(data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyNEID, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = extrema(data.λ)
   return (min=λmin, max=λmax)
end

default_λmin = 3950.0  # Based on HD solar data from PSU, should generalize
default_λmax = 9500.0  #

import ..RvSpectML: filter_line_list, find_worst_telluric_in_each_chunk
function filter_line_list(df::DataFrame, inst::IT ; λmin::Real = default_λmin, λmax::Real = default_λmax ) where { # data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:NEID.AnyNEID } #, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   df |> @filter(λmin <= _.lambda <= λmax) |>
    #    @filter( _.lambda < 6000.0 ) |>                       # Avoid tellurics at redder wavelengths
    #    @filter( _.lambda >6157 || _.lambda < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame
end
