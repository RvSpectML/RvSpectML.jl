"""
Authors: Various
Compiled by: Eric Ford (see each function for credits)
Created: August 2020
Contact: https://github.com/eford/
"""

# Constants
const speed_of_light_mps = 299792458.0 # TODO: Update value

"""
   calc_doppler_factor(rv; v_perp)

Return the Doppler boost factor (non-relativistic) for rv in km/s.
"""
function calc_doppler_factor end

calc_doppler_factor(rv::Real) = one(rv) + rv/speed_of_light_mps
calc_doppler_factor(rv::Real, v_perp::Real) = (one(rv) + rv/speed_of_light_mps)/(one(rv) - (rv^2+v_perp^2)/speed_of_light_mps^2)


"""
   absorption_line(x; mid, width, depth)

Return a Gaussian absorption line profile evaluated at x.
"""
function absorption_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:Real
    return one(T) - depth * exp(-0.5*((x-mid)/width)^2)
end

"""
Estimate line width based on stellar Teff (K) and optionally v_rot (km/s).  Output in km/s.
"""
function predict_intrinsic_stellar_line_width(Teff::Real; v_rot::Real=zero(Teff))
    @assert 3000 < Teff < 10000 # K
    @assert 0 <= v_rot <=100 # km/s
    line_width_thermal = 13*sqrt(Teff/1e4) # km/s
    line_width = sqrt(v_rot^2+line_width_thermal^2) # km/s
end


"""
   searchsortednearest(a<:AbstractVector, x::Real; assume_sorted = false )
   searchsortednearest(a<:AbstractVector, x<:AbstractVector; assume_sorted = false )

   Find the index of vector a where the value of a is closest to x.
   All vectors are assumed to already be sorted.  
   To turn off assertions, set assume_sorted to true.

Credit: traktofon @ https://discourse.julialang.org/t/findnearest-function/4143/4
Vector Vector version by Christian Gilbertson?
issorted assertion and optional assume_sorted added by Eric Ford
"""
function searchsortednearest end

function searchsortednearest(a::AbstractVector{T} where T<:Real, x::Real ; assume_sorted::Bool = false)
   if !assume_sorted
	   @assert issorted(a)
   end
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   #if (abs(a[idx]-x) < abs(a[idx-1]-x))
   if (abs2(a[idx]-x) < abs2(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end

function searchsortednearest(x::T, a::AbstractVector{T}; assume_sorted::Bool = false) where T
    @warn "Did you mean to reverse the order of x and a?"
    return searchsortednearest(a, x, assume_sorted=assume_sorted)
end

function searchsortednearest(a::AbstractVector{T1}, x::AbstractVector{T2}; assume_sorted::Bool = false) where { T1, T2 }
   if !assume_sorted
	   @assert issorted(a)
   	   @assert issorted(x)
   end
   len_x = length(x)
   len_a = length(a)
   idxs = zeros(Int64, len_x)
   idxs[1] = searchsortednearest(a, x[1])
   for i in 2:len_x
	   idxs[i] = idxs[i-1] + searchsortednearest(view(a, idxs[i-1]:len_a), x[i]) - 1
   end
   return idxs
end



"""A generalized version of the built in append!() function
By Christian Gilbertson?
# TODO:  Ask Christian what the purpose of this is relative to std append
"""
function multiple_append!(a::Vector{T}, b...) where {T<:Real}
    for i in 1:length(b)
        append!(a, b[i])
    end
    return a
end

""" Return true if all elements of array are equal to each other. """
@inline function allequal(x::AbstractArray{T,1}) where {T<:Real}
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end


""" findargminmax(a)
Return (argmin, min, argmax, max)
Adapapted from https://github.com/JuliaLang/julia/blob/697e782ab86bfcdd7fd15550241fe162c51d9f98/base/array.jl#L2191
"""
findargminmax(a) = _findminmax(a, :)
function _findminmax(a, ::Colon)
    p = pairs(a)
    y = iterate(p)
    if y === nothing
        throw(ArgumentError("collection must be non-empty"))
    end
    (mi, m), s = y
    i = mi
	argmin = mi
	valmin = m
	argmax = mi
	valmax = m
    while true
        y = iterate(p, s)
        y === nothing && break
        valmin != valmin && break
		valmax != valmax && break
        (i, ai), s = y
        if ai != ai || isless(ai, valmin)
            valmin = ai
            argmin = i
        end
		if ai != ai || isless(valmax, ai)
			valmax = ai
			argmax = i
		end

    end
    return (argmin=argmin, min=valmin, argmax=argmax, max=valmax)
end
