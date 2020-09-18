"""
    Code for LineList types used by CCfs
Author: Eric Ford
Created: August 2020
"""

"""A struct implementing a line list should be a subtype of AbstractLineList. """
abstract type AbstractLineList end

""" A basic line list for passing to compute CCFs.
Contains (views into) arrays specifying the minimum and maximum wavelength range and weight for each line. """
struct BasicLineList{T<:Real, AA<:AbstractArray{T,1} } <: AbstractLineList
    λ::AA
    weight::AA
end

""" BasicLineList( λ, weight ) """
function BasicLineList{T}(λ::AA, w::AA) where { T<:Real, AA<:AbstractArray{T,1} }
    @assert length(λ) == length(w)
    @assert length(λ) >= 1
    @assert 0.0 .<= w .<= 1.0
    BasicLineList{eltype(w),typeof(w)}(λ,w)
end

function EmptyBasicLineList()
    return BasicLineList{Float64,Array{Float64,1}}(zeros(0),zeros(0))
end
#=  Not fully implemented/tested yet
""" A line list for passing to compute CCFs with variable line widths.
Contains (views into) arrays specifying the minimum and maximum wavelength range and weight for each line. """
struct VarWidthLineList{T<:Real, AA<:AbstractArray{T,1} } <: AbstractLineList
    λ_lo::AA
    λ_hi::AA
    weight::AA
end

""" VarWidthLineList( λ_lo, λ_hi, weight ) """
function VarWidthLineList{T}(lo::AA, hi::AA, w::AA) where { T<:Real, AA<:AbstractArray{T,1} }
    @assert length(lo) == length(hi) == length(w)
    @assert length(lo) >= 1
    @assert 0.0 .<= λ_hi.-λ_lo .<= 2.0
    @assert 0.0 .<= w .<= 1.0
    VarWidthLineList{eltype(w),typeof(w)}(lo,hi,w)
end
=#

import Base.length
""" Return length of line linst. """
length(ll::AbstractLineList) = length(ll.weight)
