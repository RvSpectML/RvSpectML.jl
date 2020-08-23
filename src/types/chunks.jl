""" Abstract type for any ChunkOfSpectrum """
abstract type AbstractChuckOfSpectrum end

""" ChuckOfSpectrum for views into Spectra1DBasic or Spectra2DBasic """
struct ChuckOfSpectrum{T1<:Real,T2<:Real,T3<:Real,AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1}} <: AbstractChuckOfSpectrum
    λ::AA1
    flux::AA2
    var::AA3
end

function ChuckOfSpectrum{T1,T2,T3}(λ::A1, flux::A2, var::A3) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1} }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    min_pixels_in_chunk = 6
    max_pixels_in_chunk = 10000
    @assert min_pixels_in_chunk <= length(λ) <= max_pixels_in_chunk
    ChuckOfSpectrum{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var)}(λ,flux,var)
end

function ChuckOfSpectrum(λ::A1, flux::A2, var::A3, order::Integer, pixels::AUR) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= order <= size(λ,2)
    @assert 1 <= first(pixels) < last(pixels) <= size(λ,1)
    ChuckOfSpectrum{T1,T2,T3}(view(λ,pixels,order),view(flux,pixels,order),view(var,pixels,order))
end


function ChuckOfSpectrum(λ::A1, flux::A2, var::A3, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange, I1<:Integer }
    ChuckOfSpectrum(λ,flux,var,loc.order,loc.pixels)
end

function ChuckOfSpectrum(spectra::AS, order::Integer, pixels::AUR) where { AS<:AbstractSpectra, AUR<:AbstractUnitRange }
    ChuckOfSpectrum(spectra.λ,spectra.flux,spectra.var,order,pixels)
end

function ChuckOfSpectrum(spectra::AS, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  AS<:AbstractSpectra, AUR<:AbstractUnitRange, I1<:Integer }
    ChuckOfSpectrum(spectra.λ,spectra.flux,spectra.var,loc.order,loc.pixels)
end

abstract type AbstractChunkList end

""" Struct containing an array of ChuckOfSpectrum """
mutable struct ChunkList{CT<:AbstractChuckOfSpectrum, AT<:AbstractArray{CT,1} } <: AbstractChunkList
      data::AT
end

#=
function ChunkList(in::AT) where {CT<:AbstractChuckOfSpectrum, AT<:AbstractArray{CT,1} }
    ChunkList{CT,AT}(in)
end
=#

abstract type AbstractChunkListTimeseries end

""" Mtching lists of times and array of ChunkLists """
mutable struct ChunkListTimeseries{CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1}, TT<:Real, AT<:AbstractArray{TT,1} } <: AbstractChunkListTimeseries
    times::AT
    chunk_list::ACLT
end

#=
function ChunkListTimeseries(t::AT, cl::ACLT) where {CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1}, TT<:Real, AT<:AbstractArray{TT,1} }
    ChunkListTimeseries{CLT,ACLT,TT,AT}(t,cl)
end
=#

import Base.length
length(cl::CLT) where {CLT<:AbstractChunkList} = length(cl.data)
length(clts::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(clts.chunk_list)

""" Number of chunks in first chunk_list in chunk_list_timeseries. """
num_chunks(chunk_list_timeseries::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(first(chunk_list_timeseries.chunk_list))
