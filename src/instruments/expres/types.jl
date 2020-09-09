""" Struct for 2D spectra specialized for EXPRES and its additional data products.
    Drafted, but not yet folded into examples, as still using BasicSpectra2D instead.
    For now extra info is showing up in metadata.
    """
struct Spectra2DEXPRES{T1<:Real,T2<:Real,T3<:Real,
     AA1<:AbstractArray{T1,2},AA2<:AbstractArray{T2,2},AA3<:AbstractArray{T3,2},
     InstT<:AbstractInstrument
    } <: AbstractSpectra2D
    λ::AA1
    flux::AA2
    var::AA3
    blaze::AA2
    continuum::AA2
    tellurics::AA2
    inst::InstT
    metadata::Dict{Symbol,Any} # Dict{Symbol,Any}
end


function Spectra2DEXPRES(λ::A1, flux::A2, var::A3, inst::InstT; blaze::A2, continuum::A2, tellurics::A2,
     metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() ) where {
     T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
     InstT<:AbstractInstrument  }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= size(λ,1) <= max_pixel_in_order(inst)-min_pixel_in_order(inst)+1
    @assert 1 <= size(λ,2) <= max_order(inst)-min_order(inst)+1
    Spectra2DEXPRES{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var),typeof(inst)}(λ,flux,var,blaze,continuum,tellurics,inst,metadata)
end
