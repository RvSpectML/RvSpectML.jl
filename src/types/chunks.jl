"""
Code declaring ChunkOfSpectrum, ChunkList, ChunkListTimeseries, and their abstract versions.

Author: Eric Ford
Created: August 2020
"""

""" Abstract type for any ChunkOfSpectrum """
abstract type AbstractChunkOfSpectrum end

""" ChunkOfSpectrum for views into Spectra1DBasic or Spectra2DBasic """
struct ChunkOfSpectrum{T1<:Real,T2<:Real,T3<:Real,AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1}} <: AbstractChunkOfSpectrum
    λ::AA1
    flux::AA2
    var::AA3
end

function ChunkOfSpectrum{T1,T2,T3}(λ::A1, flux::A2, var::A3) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1} }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    min_pixels_in_chunk = 3
    max_pixels_in_chunk = 10000
    @assert min_pixels_in_chunk <= length(λ) <= max_pixels_in_chunk
    ChunkOfSpectrum{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var)}(λ,flux,var)
end

function ChunkOfSpectrum(λ::A1, flux::A2, var::A3, order::Integer, pixels::AUR) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= order <= size(λ,2)
    @assert 1 <= first(pixels) < last(pixels) <= size(λ,1)
    ChunkOfSpectrum{T1,T2,T3}(view(λ,pixels,order),view(flux,pixels,order),view(var,pixels,order))
end

function ChunkOfSpectrum(λ::A1, flux::A2, var::A3, pixels::AUR) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1}, AUR<:AbstractUnitRange }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= first(pixels) < last(pixels) <= size(λ,1)
    ChunkOfSpectrum{T1,T2,T3}(view(λ,pixels),view(flux,pixels),view(var,pixels))
end


function ChunkOfSpectrum(λ::A1, flux::A2, var::A3, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange, I1<:Integer }
    ChunkOfSpectrum(λ,flux,var,loc.order,loc.pixels)
end

function ChunkOfSpectrum(spectra::AS, order::Integer, pixels::AUR) where { AS<:AbstractSpectra2D, AUR<:AbstractUnitRange }
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,order,pixels)
end

function ChunkOfSpectrum(spectra::AS, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  AS<:AbstractSpectra2D, AUR<:AbstractUnitRange, I1<:Integer }
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,loc.order,loc.pixels)
end

function ChunkOfSpectrum(spectra::AS, pixels::AUR) where { AS<:AbstractSpectra1D, AUR<:AbstractUnitRange }
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,pixels)
end


get_order_index(chunk::AC) where { AC<:AbstractChunkOfSpectrum }  = chunk.flux.indices[2]
get_pixels_range(chunk::AC) where { AC<:AbstractChunkOfSpectrum } = chunk.flux.indices[1]



#AbstractSpectra2D

abstract type AbstractChunkList end

""" Struct containing an array of ChunkOfSpectrum """
struct ChunkList{CT<:AbstractChunkOfSpectrum, AT<:AbstractArray{CT,1} } <: AbstractChunkList
      data::AT

      #=
      function ChunkList{CT,AT}(in::AT) where {CT<:AbstractChunkOfSpectrum, AT<:AbstractArray{CT,1} }
        @assert length(data)>=1
        ChunkList{CT,AT}(in)
    end
    =#
end

import Base.getindex, Base.setindex!, Base.append!

""" Allow direct access to data, an AbstractArray of ChunkOfSpectrum's via [] operator """
function getindex(x::CLT,idx) where {CLT<:AbstractChunkList}
    x.data[idx]
end

""" Allow direct access to data, an AbstractArray of ChunkOfSpectrum's via [] operator """
function setindex!(x::CLT, idx, y::CLT) where {CLT<:AbstractChunkList}
    x.data[idx] .= y.data[idx]
end

function append!(x::CLT, y::CLT) where {CLT<:AbstractChunkList}
    append!(x.data,y.data)
end


abstract type AbstractChunkListTimeseries end

""" Mtching lists of times and array of ChunkLists """
struct ChunkListTimeseries{ TT<:Real, AT<:AbstractArray{TT,1}, CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1}, InstT<:AbstractInstrument, AVMT<:AbstractVector{MetadataT} } <: AbstractChunkListTimeseries
    times::AT
    chunk_list::ACLT
    inst::InstT
    metadata::AVMT
end

function ChunkListTimeseries(times::AT, clt::ACLT; inst::InstT = Generic1D(),
            metadata::AbstractVector{MetadataT} = fill(MetadataT,length(times)) ) where {
                TT<:Real, AT<:AbstractArray{TT,1}, CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1},  InstT<:AbstractInstrument, AVMT<:AbstractVector{MetadataT} }
    ChunkListTimeseries{eltype(times),typeof(times),eltype(clt),typeof(clt),typeof(inst),typeof(metadata)}(times,clt,inst,metadata)
end

#=
function ChunkListTimeseries(t::AT, cl::ACLT) where {CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1}, TT<:Real, AT<:AbstractArray{TT,1} }
    ChunkListTimeseries{CLT,ACLT,TT,AT}(t,cl)
end
=#

import Base.length
""" Return number of chunks in ChunkList """
length(cl::CLT) where {CLT<:AbstractChunkList} = length(cl.data)

""" Return number of times/ChunkLists in a ChunkListTimeseries """
length(clts::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(clts.chunk_list)
""" Return number of times/ChunkLists in a ChunkListTimeseries """
num_times(clts::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(clts.chunk_list)

""" Number of chunks in first chunk_list in chunk_list_timeseries. """
num_chunks(chunk_list_timeseries::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(first(chunk_list_timeseries.chunk_list))
