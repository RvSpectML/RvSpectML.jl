#=
const min_orders_in_spectra = 1
const max_orders_in_spectra = 128
const min_pixels_in_order = 1
const max_pixels_in_order = 9128
const min_pixels_in_spectra = min_orders_in_spectra*min_pixels_in_order
const max_pixels_in_spectra = max_orders_in_spectra*max_pixels_in_order
min_usable_pixels_in_order = 128
=#

""" Abstract type for any Spectrum (or region of spectrum) """
abstract type AbstractSpectra end
""" Abstract type for any 1-d Spectrum (or region of spectrum) """
abstract type AbstractSpectra1D <: AbstractSpectra end
""" Abstract type for any 2-d Spectrum """
abstract type AbstractSpectra2D <: AbstractSpectra end

""" Basic struct for Spectra1D (or region of specturm)
    Instruments can specialize their own if additional data is avaliable. """
struct Spectra1DBasic{T1<:Real,T2<:Real,T3<:Real,#=T4<:Real,=# AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1} } <: AbstractSpectra1D
    λ::AA1
    flux::AA2
    var::AA3
    #doppler_factor::T4    # Move to metadata
    metadata::Dict{Symbol,Any}
end

""" Basic struct for Spectra2D (or region of specturm)
    Instruments can specialize their own if additional data is avaliable. """
struct Spectra2DBasic{T1<:Real,T2<:Real,T3<:Real,#=T4<:Real,=# AA1<:AbstractArray{T1,2},AA2<:AbstractArray{T2,2},AA3<:AbstractArray{T3,2}} <: AbstractSpectra2D
    λ::AA1
    flux::AA2
    var::AA3
    #doppler_factor::T4    # Move to metadata
    metadata::Dict{Symbol,Any}
end


function Spectra1DBasic(λ::A1, flux::A2, var::A3; #= doppler_factor::T4=one(eltype(λ)), =# metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() ) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1}} #, T4<:Real }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert min_pixels_in_spectra <= length(λ) <= max_pixels_in_spectra
    Spectra1DBasic{eltype(λ),eltype(flux),eltype(var),#=typeof(doppler_factor),=# typeof(λ),typeof(flux),typeof(var)}(λ,flux,var,#=doppler_factor,=# metadata)
end


function Spectra2DBasic(λ::A1, flux::A2, var::A3; #= doppler_factor::T4=one(eltype(λ)), =# metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() ) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}} #, T4<:Real }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert min_pixels_in_spectra <= size(λ,1) <= max_pixels_in_spectra
    @assert min_orders_in_spectra <= size(λ,2) <= max_orders_in_spectra
    Spectra2DBasic{eltype(λ),eltype(flux),eltype(var),#=typeof(doppler_factor),=# typeof(λ),typeof(flux),typeof(var)}(λ,flux,var,#=doppler_factor,=# metadata)
end
