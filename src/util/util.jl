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

Return the Doppler boost factor (non-relativistic) for rv in m/s.
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
Estimate line width based on stellar Teff (K) and optionally v_rot (m/s).  Output in m/s.
"""
function predict_intrinsic_stellar_line_width(Teff::Real; v_rot::Real=zero(Teff))
    @assert 3000 < Teff < 10000 # K
    @assert 0 <= v_rot <=100e3 # m/s
    line_width_thermal = 13e3*sqrt(Teff/1e4) # m/s
    line_width = sqrt(v_rot^2+line_width_thermal^2) # m/s
end


default_Δλ_over_λ_threshold_check_if_line_match = 2.25e-5
""" check_if_line_match ( λ, list ; threshold )
Return true if list contains a wavelength differing from λ by no more than threshold (in units of Δλ/λ)
"""
function check_if_line_match( λ::Real, list::AbstractVector{T} ; threshold::Real = default_Δλ_over_λ_threshold_check_if_line_match ) where { T<:Real }
  idx = searchsortednearest(list, λ)
  Δ = abs(list[idx]-λ)/λ
  if Δ <= threshold
    return true
  else
    return false
  end
end

""" find_which_line_fits_in_line_list( fit_list, line_list; threshold )
Return list of Bools indicatin which line(s) from fit_list match a line in line_list to within threshold (in units of Δλ/λ)
Warning: Untested
"""
function find_which_line_fits_in_line_list(fit_list::AbstractVector{T1}, line_list::AbstractVector{T2} ; threshold::Real = default_Δλ_over_λ_threshold_check_if_line_match ) where { T1<:Real, T2<:Real }
  @assert issorted(line_list)
  perm = sortperm(fit_list)
  idx = RvSpectML.searchsortednearest(line_list, fit_list[perm])
  Δ = Array{Float64,1}(undef, size(df,1) )
  Δ[perm] .= abs.(line_list[idx] .- fit_list[perm])./fit_list[perm]
  line_matches = Δ .<= threshold
  return line_matches
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

""" calc_line_width(λ, flux; frac_depth )
Returns the line width (units of λ) for specified fractional line depth (default of 0.5).
Assumes continuum is the maximum flux provided.
"""
function calc_line_width(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}; frac_depth::Real = 0.5 ) where { T1<:Real, T2<:Real }
   @assert 0.05 <= frac_depth <= 0.99
   idx_min_flux = argmin(flux)
   min_flux = flux[idx_min_flux]
   continuum = maximum(flux)
   depth = 1.0 - min_flux/continuum
   target_flux = continuum*(1-frac_depth*depth)
   idxhi = idx_min_flux-1+searchsortedfirst(view(flux,idx_min_flux:length(flux)), target_flux)
   if !(min_flux<=idxhi<=length(flux))   idx = length(flux) end
   idxlo = idxhi-1
   λ_hi = RvSpectML.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
   idxlo = idx_min_flux+1-searchsortedfirst(view(flux,idx_min_flux:-1:1), target_flux )
   if !(1<=idxlo<=min_flux)   idx = 1 end
   idxhi = idxlo+1
   λ_lo = RvSpectML.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
   width = λ_hi-λ_lo
   return width
end

""" calc_line_bisector_at_frac_depth(λ, flux, frac_depth )
Returns the line average of wavelengths (units of λ) at specified fractional line depth.
Assumes continuum is the maximum flux provided.
"""
function calc_line_bisector_at_frac_depth(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, frac_depth::Real ) where { T1<:Real, T2<:Real }
   @assert 0.05 <= frac_depth <= 0.99
   idx_min_flux = argmin(flux)
   min_flux = flux[idx_min_flux]
   continuum = maximum(flux)
   depth = 1.0 - min_flux/continuum
   abs_depth = frac_depth*depth

   target_flux = continuum*(1-frac_depth*depth)
   idxhi = idx_min_flux-1+searchsortedfirst(view(flux,idx_min_flux:length(flux)), target_flux)
   idxlo = idxhi-1
   λ_hi = RvSpectML.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
   idxlo = idx_min_flux+1-searchsortedfirst(view(flux,idx_min_flux:-1:1), target_flux )
   idxhi = idxlo+1
   λ_lo = RvSpectML.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
   bisector = (λ_hi+λ_lo)/2
   return bisector
end

""" calc_line_bisector_at_frac_depth(λ, flux, abs_depth )
Returns the line average of wavelengths (units of λ) at specified absolute line depth.
Assumes continuum is the maximum flux provided.
"""
function calc_line_width(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, abs_depth::Real ) where { T1<:Real, T2<:Real }
   @assert 0.05 <= abs_depth <= 0.99
   idx_min_flux = argmin(flux)
   continuum = maximum(flux)
   if flux[idx_min_flux]/continuum > abs_depth  # Line isn't that deep!
	   return nothing
   end
   target_flux = continuum*(1-abs_depth)
   idxhi = idx_min_flux-1+searchsortedfirst(view(flux,idx_min_flux:length(flux)), target_flux)
   idxlo = idxhi-1
   λ_hi = RvSpectML.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
   idxlo = idx_min_flux+1-searchsortedfirst(view(flux,idx_min_flux:-1:1), target_flux )
   idxhi = idxlo+1
   λ_lo = RvSpectML.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
   width = λ_hi-λ_lo
   return width
end


""" interp_linear(;x1::T1,x2::T1,y1::T2,y2::T2,xpred::T1)
Return result of simple linear interpolant at xpred.
Does not test that xpred is between x1 and x2.
"""
function interp_linear(;x1::T1,x2::T1,y1::T2,y2::T2,xpred::T1) where { T1<:Real, T2<:Real }
  ypred = y1+(y2-y1)*((xpred-x1)/(x2-x1))
end
