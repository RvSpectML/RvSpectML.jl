"""
Code declaring Spectra1DBasic, Spectra1DBasic, SpectralTimeSeriesCommonWavelengths and their abstract versions.

Author: Eric Ford
Created: August 2020
"""

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

MetadataT = Dict{Symbol,Any}

""" Basic struct for Spectra1D (or region of specturm)
    Instruments can specialize their own if additional data is avaliable. """
struct Spectra1DBasic{T1<:Real,T2<:Real,T3<:Real,
    AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1},
    InstT<:AbstractInstrument
    } <: AbstractSpectra1D
    λ::AA1
    flux::AA2
    var::AA3
    inst::InstT
    metadata::MetadataT #Dict{Symbol,Any}
end

""" Basic struct for Spectra2D (or region of specturm)
    Instruments can specialize their own if additional data is avaliable. """
struct Spectra2DBasic{T1<:Real,T2<:Real,T3<:Real,
     AA1<:AbstractArray{T1,2},AA2<:AbstractArray{T2,2},AA3<:AbstractArray{T3,2},
     InstT<:AbstractInstrument
    } <: AbstractSpectra2D
    λ::AA1
    flux::AA2
    var::AA3
    inst::InstT
    metadata::MetadataT # Dict{Symbol,Any}
end


function Spectra1DBasic(λ::A1, flux::A2, var::A3, inst::InstT;
        metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() ) where {
          T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1},
          InstT<:AbstractInstrument }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= length(λ) <= max_pixel(inst)-min_pixel(inst)+1
    Spectra1DBasic{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var),typeof(inst)}(λ,flux,var,inst,metadata)
end

function Spectra2DBasic(λ::A1, flux::A2, var::A3, inst::InstT;
     metadata::MetadataT = MetadataT() ) where {
     T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
     InstT<:AbstractInstrument  }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= size(λ,1) <= max_pixel_in_order(inst)-min_pixel_in_order(inst)+1
    @assert 1 <= size(λ,2) <= max_order(inst)-min_order(inst)+1
    Spectra2DBasic{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var),typeof(inst)}(λ,flux,var,inst,metadata)
end


""" Abstract type for a time series of spectra that share a common wavelength grid. """
abstract type AbstractSpectralTimeSeriesCommonWavelengths <: AbstractSpectra1D   end

""" Time series of spectra that share a common wavelength grid. """
struct SpectralTimeSeriesCommonWavelengths{T1<:Real,T2<:Real,T3<:Real,AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,2},AA3<:AbstractArray{T3,2},
            AA4<:AbstractArray{UnitRange{Int64},1}, InstT<:AbstractInstrument } <: AbstractSpectralTimeSeriesCommonWavelengths
    λ::AA1
    flux::AA2
    var::AA3
    chunk_map::AA4
    inst::InstT
    metadata::MetadataT
end

function SpectralTimeSeriesCommonWavelengths(λ::A1, flux::A2, var::A3, chunk_map::A4, inst::InstT;
        metadata::MetadataT = MetadataT() ) where {
          T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
          A4<:AbstractArray{UnitRange{Int64},1}, InstT<:AbstractInstrument }
          println("len(λ) = ", length(λ))
          println("size(flux) = ", size(flux))
    @assert length(λ) == size(flux,1)
    @assert length(λ) == size(var,1)
    @assert 1 <= length(λ)
    SpectralTimeSeriesCommonWavelengths{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var),typeof(chunk_map),typeof(inst)}(λ,flux,var,chunk_map,inst,metadata)
end

""" Extract the metadata from a time series of spectra and return it as an array. """
function make_vec_metadata_from_spectral_timeseries(spec_arr::AA) where { AS<:AbstractSpectra, AA<:AbstractArray{AS,1} }
    map(s->s.metadata,spec_arr)
end
