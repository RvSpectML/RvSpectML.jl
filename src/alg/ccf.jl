"""
    Code to compute CCF
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactors and optimized by Eric Ford
"""

""" Module for computing CCFs """
module CCF

export AbstractCCFPlan, BasicCCFPlan
export AbstractCCFMaskShape, TopHatCCFMask
export AbstractLineList, BasicLineList

export calc_ccf_v_grid, ccf_1D, ccf_1D!
export calc_ccf_chunk, calc_ccf_chunklist, calc_ccf_chunklist_timeseries

import ..RvSpectML
import ..RvSpectML: AbstractChuckOfSpectrum, AbstractChunkList, AbstractChunkListTimeseries
import ..RvSpectML: AbstractInstrument #, default_ccf_v_width

using ThreadedIterables

#using ThreadTools

# physical constants
const c_ms = 2.99782458e8    # m/s
const c_kms = 2.99782458e5   # km/s

# default parameters
const default_v_center = 0.0    # m/s
const default_v_step = 250.0    # m/s
const default_v_max = 15.0e3    # m/s
const default_v_width = 410.0   # m/s

"""A struct implementing a specific plans describing where the CCF is to be evaluated should be a subtype of AbstractCCFPlan. """
abstract type AbstractCCFPlan end

""" Basic plan for computing the CCF roughly between v_center-v_max and v_center+v_max with step size v_step. """
struct BasicCCFPlan <: AbstractCCFPlan
    v_center::Float64
    v_step::Float64
    v_max::Float64

    function BasicCCFPlan(midpoint::Real,step::Real, max::Real)
        @assert 1.0e3 <= max <=100.0e3    # Reasonable range m/s, designed to prevent mistakes
        @assert 1 < step < 1000           # Reasonable range m/s
        @assert abs(midpoint) < max
        new(midpoint,step,max)
    end

end

"""   BasicCCFPlan
# Optional arguments:
- midpoint: (default_v_center)
- step: (default_v_step)
- max: (default_v_max)
"""
function BasicCCFPlan(;midpoint::Real=default_v_center,
            step::Real=default_v_step, max::Real=default_v_max )
    BasicCCFPlan(midpoint,step,max)
end

""" calc_ccf_v_grid( plan )
Return range with 2n+1 points between -v_max and v_max where CCF is to be evaluated.
Units based on those in plan.
"""
function calc_ccf_v_grid(p::AP where AP<:BasicCCFPlan )
    n = ceil(Int, p.v_max/p.v_step)
    range(p.v_center-n*p.v_step, p.v_center+n*p.v_step, length=2*n+1)
end


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

import Base.length
""" Return length of line linst. """
length(ll::AbstractLineList) = length(ll.weight)


"A struct implementing a specific mask shapes should be a subtype of AbstractCCFMaskShape."
abstract type AbstractCCFMaskShape  end

"""   TopHatCCFMask
The standard tophat mask with one parameter, it's width as a velocity in m/s.
Mask weights are stored separately in a line list.
"""
struct TopHatCCFMask <: AbstractCCFMaskShape
    half_width::Float64
end

""" TopHatCCFMask( ; half_width=default_v_width ) """
function TopHatCCFMask(w::Real=default_v_width)
    @assert 0 < w < default_v_max
    TopHatCCFMask(w/2)
end

""" TopHatCCFMask( inst  ) """
function TopHatCCFMask(inst::AbstractInstrument; scale_factor::Real = 1)
    w = scale_factor * RvSpectML.default_ccf_v_width(inst) # / c_ms
    TopHatCCFMask(w/2)
end

λ_min(m::TopHatCCFMask,λ::Real) = λ/calc_doppler_factor(m.half_width)
λ_max(m::TopHatCCFMask,λ::Real) = λ*calc_doppler_factor(m.half_width)

""" Functor for returning 1 for any Δv <= width.  """
function (m::TopHatCCFMask)(Δv::Real)
    return abs2(Δv)<=abs2(m.half_width) ? 1 : 0
end

"""
    project_mask!( output, λs, line_list; shift_factor, mask_shape )

Compute the projection of the mask onto the 1D array of wavelengths (λs) at a given shift factor (default: 1).
The mask is computed from a line_list and mask_shape (default: tophat).

TODO: Implement/test other mask_shapes.
"""
function project_mask_opt!(projection::A2, λs::A1, line_list::ALL; shift_factor::Real=one(T1),
            mask_shape::ACMS = TopHatCCFMask() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2},
                ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }

    # find bin edges
    λsle = find_bin_edges_opt(λs)  # Read as λ's left edges
    # TODO: OPT: Can get remove this allocation and compute these as needed, at least for tophat mask.  But might be useful to keep for more non-const masks.

    # allocate memory for mask projection
    nm = length(line_list) # size(mask, 1)
    nx = size(λsle, 1)-2
    ny = size(λsle, 2)
    @assert ny == 1              # TODO:  Ask Alex what purpose of this is.  Maybe running CCF with multiple masks from the same lines at once?
    @assert size(projection,1) == nx
    @assert size(projection,2) == ny
    #projection = zeros(nx, ny)
    projection .= zero(eltype(projection))  # Important to zero out projection before accumulating there!

    # set indexing variables
    p = 1
    λsle_cur = λsle[p]   # Current pixel's left edge
    λsre_cur = λsle[p+1] # Current pixel's right edge
    m = 1
    on_mask = false
    #mask_lo = line_list.λ_lo[m] * shift_factor   # current line's lower limit
    #mask_hi = line_list.λ_hi[m] * shift_factor   # current line's upper limit
    mask_lo = λ_min(mask_shape,line_list.λ[m]) * shift_factor   # current line's lower limit
    mask_hi = λ_max(mask_shape,line_list.λ[m]) * shift_factor   # current line's upper limit

    mask_weight = line_list.weight[m]            # current liine's weight
    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if λsre_cur > mask_lo
                if λsre_cur > mask_hi   # Pixel overshot this mask entry, so weight based on mask filling only a portion of the pixel,
                    projection[p] += mask_weight * (mask_hi - mask_lo) / (λsre_cur - λsle_cur)
                    m += 1
                    if m<=length(line_list)      # Move to next line
                        mask_lo = λ_min(mask_shape,line_list.λ[m]) * shift_factor
                        mask_hi = λ_max(mask_shape,line_list.λ[m]) * shift_factor
                        mask_weight = line_list.weight[m]
                    else         # We're out of lines, so exit

                        break
                    end
                else                       # Right edge of current pixel has entered the mask for this line, but left edge hasn't
                    projection[p] += mask_weight * (λsre_cur - mask_lo) / (λsre_cur - λsle_cur)
                    on_mask = true         # Indicate next pixel will hav something to contribute based on the current line
                    p+=1                   # Move to next pixel
                    λsle_cur = λsle[p]
                    λsre_cur = λsle[p+1]
                end
            else   # ALl of pixel is still to left of beginning of mask for this line.
                p+=1
                λsle_cur = λsle[p]
                λsre_cur = λsle[p+1]
            end
        else
            if λsre_cur > mask_hi               # Right edge of this pixel moved past the edge of the current line's mask
                projection[p] += mask_weight * (mask_hi - λsle_cur) / (λsre_cur - λsle_cur)
                on_mask = false                 # Indicate that we're done with this line
                m += 1                          # Move to next line
                if m<=length(line_list)
                    mask_lo = λ_min(mask_shape,line_list.λ[m]) * shift_factor
                    mask_hi = λ_max(mask_shape,line_list.λ[m]) * shift_factor
                    mask_weight =  line_list.weight[m]
                else                            # We're done with all lines, can return early
                    break
                end
            else                                # Mask window is entirely within this pixel
                projection[p] += mask_weight    # Add the full mask weight
                p += 1                          # move to next pixel
                λsle_cur = λsle[p]
                λsre_cur = λsle[p+1]
            end
        end
        if p>length(projection) break end       # We're done with all pixels in this chunk.
    end
    return projection
end



"""
    ccf_1D!(ccf_out, λs, fluxes, line_list; mask_shape, plan )

Compute the cross correlation function of a spectrum with a mask.
# Inputs:
- ccf_out: 1-d array of size length(calc_ccf_v_grid(plan)) to store output
- λs: 1-d array of wavelengths
- fluxes:  1-d array of fluxes
- line_linst:  Each mask entry consists of the left and right end of tophats and a weight.
# Optional Arguments:
- mask_shape: shape of mask to use (TopHatCCFMask())
- plan:  parameters for computing ccf (BasicCCFPlan())
"""
function ccf_1D!(ccf_out::A1, λ::A2, flux::A3, line_list::ALL, #mask_shape1::A3
                ; mask_shape::ACMS = TopHatCCFMask(), plan::AP = BasicCCFPlan() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape, AP<:AbstractCCFPlan }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    v_grid = calc_ccf_v_grid(plan)
    #ccf_out = zeros(size(v_grid))
    @assert size(ccf_out,1) == length(v_grid)

    mask_projections = zeros(length(λ),1)

    # loop through each velocity, project the mask and compute ccf at given v
    for i in 1:size(ccf_out,1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask_opt!(mask_projections, λ, line_list, shift_factor=doppler_factor)

        # compute the ccf value at the current velocity shift
        ccf_out[i] = sum(flux .* mask_projections)
    end
    return ccf_out
end


"""
    ccf_1D!(ccf_out, λs, fluxes, line_list; mask_shape, plan )

Compute the cross correlation function of a spectrum with a mask.
# Inputs:
- λs: 1-d array of wavelengths
- fluxes:  1-d array of fluxes
- line_linst:  Each mask entry consists of the left and right end of tophats and a weight.
# Optional Arguments:
- mask_shape: shape of mask to use (TopHatCCFMask())
- plan:  parameters for computing ccf (BasicCCFPlan())
# Returns:
- 1-d array of size length(calc_ccf_v_grid(plan))
"""
function ccf_1D(λ::A2, flux::A3, line_list::ALL, #mask_shape1::A3
                ; mask_shape::ACMS = TopHatCCFMask(), plan::AP = BasicCCFPlan() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape, AP<:AbstractCCFPlan }
    v_grid = calc_ccf_v_grid(plan)
    ccf_out = zeros(size(v_grid))
    ccf_1D!(ccf_out, λ, flux, line_list, mask_shape=mask_shape, plan=plan )
    return ccf_out
end


"""
    calc_doppler_factor(vel)

Compute the longitudinal relativistic doppler factor given a velocity
in meters per second.
"""
function calc_doppler_factor(vel::Real)
    one(vel) + vel/c_ms
end

#=
 Eric replaced this version with above.  Should double check that won't cause anyeone suprises.
function calc_doppler_factor_old(vel::Real)
    num = one(vel) + vel/c_ms
    den = one(vel) - vel/c_ms
    return sqrt(num/den)
end
=#

"""  calc_ccf_chunk( chunk, line_list::ALL )
Convenience function to compute CCF for one chunk of spectrum.
# Inputs:
- chunk
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
CCF for one chunk of spectrum, evaluated using mask_shape and plan
"""
function calc_ccf_chunk(chunk::AbstractChuckOfSpectrum,
                            line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan::AbstractCCFPlan = BasicCCFPlan() ) where {
                                    ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  v_grid = calc_ccf_v_grid(plan)
  ccf_out = zeros(size(v_grid))
  ccf_1D!(ccf_out, chunk.λ, chunk.flux, line_list, mask_shape=mask_shape, plan=plan)
  return ccf_out
end

"""  calc_ccf_chunk( chunklist, line_list )
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using mask_shape and plan.
"""
function calc_ccf_chunklist(chunk_list::AbstractChunkList,
                            line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan::AbstractCCFPlan = BasicCCFPlan() ) where {
                                    ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  mapreduce(chunk->calc_ccf_chunk(chunk, line_list,mask_shape=mask_shape, plan=plan), +, chunk_list.data)
end

"""  calc_ccf_chunk( chunklist_timeseries, line_list )
Convenience function to compute CCF for a timeseries of spectra, each with a chunklist.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using mask_shape and plan.
"""
function calc_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                            line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan::AbstractCCFPlan = BasicCCFPlan() ) where {
                                    ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  @threaded  mapreduce(cl->calc_ccf_chunklist(cl, line_list, mask_shape=mask_shape, plan=plan),hcat,clt.chunk_list)
end

"""  calc_ccf_chunk( chunklist_timeseries, line_list )
Convenience function to compute separate CCFs for each chunk in a spectrum.
CCF evaluated using mask_shape and plan.
# Inputs:
- chunklist_timeseries
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccfs_chunklist(chunk_list::AbstractChunkList, line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
    plan::AbstractCCFPlan = BasicCCFPlan() ) where {
        ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
    mapreduce(chunk->calc_ccf_chunk(chunk, line_list, mask_shape=mask_shape, plan=plan), hcat, chunk_list.data)
end

"""  calc_ccf_chunk( chunklist_timeseries, line_list )
Convenience function to compute separate CCFs for each chunk of each spectrum in a timeseries.
CCF is evaluated using mask_shape and plan.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
A 3-d array containing the CCF at each (velocity, chunk, spectrum)
"""
function calc_order_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                            line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan::AbstractCCFPlan = BasicCCFPlan() ) where {
                                    ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  nvs = length(calc_ccf_v_grid(plan))
  norders = length(clt.chunk_list[1].data)
  nobs =  length(clt.chunk_list)
  order_ccfs = zeros(nvs, norders, nobs)
  Threads.@threads for i in 1:nobs
      order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], line_list,  mask_shape=mask_shape, plan=plan)
  end
  return order_ccfs
end

"""
   find_bin_edges( pixel_centers )

Internal function used by project_mask!.
"""
function find_bin_edges_opt(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = Array{eltype(fws),1}(undef,length(fws)+2)
    fwse[2:end-2] = (fws[2:end] + fws[1:end-1]) ./ 2.0
    fwse[1] = 2.0 * fws[1] - fwse[2]
    fwse[end-1] = 2.0 * fws[end] - fwse[end-2]
    fwse[end] = zero(eltype(A))
    return fwse
end


function find_bin_edges_compare(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = (fws[2:end] + fws[1:end-1]) ./ 2.0
    le = 2.0 * fws[1] - fwse[1]
    re = 2.0 * fws[end] - fwse[end]
    fwsle = vcat(le, fwse)
    fwsre = vcat(fwse, re)
    return fwsle, fwsre
end



function project_mask_compare!(projection::A2, λs::A1, mask_in::ALL; shift_factor::Real=one(T1),
            mask_shape::ACMS = TopHatCCFMask() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2},
                ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }

#function project_mask_compare(fws::A1, mask::A2; shift_factor::Real=one(T1)) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    # find bin edges
    fwsle, fwsre = find_bin_edges_compare(λs)

    # allocate memory for line_list projection
    nm = length(mask_in)
    nx = size(fwsle, 1)
    ny = size(fwsle, 2)
    projection = zeros(nx, ny)
    #flush(stdout);  println("size(projection) = ", size(projection), " size(fws) = ", size(fws), " size(mask) = ", size(mask), " shift_factr = ", shift_factor)

    # shift the mask
    mask_shifted = zeros(length(mask_in),3)
    mask_shifted[:,1] = mask_in.λ_lo  .* shift_factor  # view(mask,:,1) .* shift_factor
    mask_shifted[:,2] = mask_in.λ_hi  .* shift_factor  # view(mask,:,2) .* shift_factor
    mask_shifted[:,3] = mask_in.weight                 # view(mask,:,3)
    local mask = mask_shifted

    # set indexing variables
    p = 1
    m = 1
    on_mask = false

    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if fwsre[p] > mask[m,1]
                if fwsre[p] > mask[m,2]
                    projection[p] += mask[m,3] * (mask[m,2] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    m += 1
                else
                    projection[p] += mask[m,3] * (fwsre[p] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    on_mask = true
                    p+=1
                end
            else
                p+=1
            end
        else
            if fwsre[p] > mask[m,2]
                projection[p] += mask[m,3] * (mask[m,2] - fwsle[p]) / (fwsre[p] - fwsle[p])
                on_mask = false
                m += 1
            else
                projection[p] += mask[m,3]
                p += 1
            end
        end
        if p>length(projection) break end
    end
    return projection
end


#=
# set some CCF parameters
# TODO: Put into a struct
#=
const CCF_width = 15e3       # m/s window
const dV = 250.0             # size of velocity step
const mask_width = 410.0     # m/s (width of mask entries)
=#



"""
    project_mask(fws, mask, shift_factor)

Compute the projection of the tophat mask onto the 1D array
of wavelengths (fws) at a given shift factor.
"""
function project_mask(fws::A1, mask::A2; shift_factor::Real=one(T1), projection::A2 = zeros(length(fws,1))) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    @assert size(mask,2) == 3
    #project_mask_compare(fws,mask,shift_factor=shift_factor)
    project_mask_opt(fws,mask,shift_factor=shift_factor,projection=projection)
end

function project_mask_opt(fws::A1, mask::A2; shift_factor::Real=one(T1),
            projection::A2 = zeros(length(fws),1) ) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    # find bin edges
    fwsle = find_bin_edges_opt(fws)

    # allocate memory for mask projection
    nm = size(mask, 1)
    nx = size(fwsle, 1)-2
    ny = size(fwsle, 2)
    @assert ny == 1              # TODO:  Ask Alex what purpose of this is.  Maybe running CCF with multiple masks from the same lines at once?
    @assert size(projection,1) == nx
    @assert size(projection,2) == ny
    #projection = zeros(nx, ny)
    projection .= zero(eltype(projection))

    # set indexing variables
    p = 1
    fwsle_cur = fwsle[p]
    fwsre_cur = fwsle[p+1]
    m = 1
    on_mask = false
    mask_lo = mask[m,1] * shift_factor
    mask_hi = mask[m,2] * shift_factor
    mask_weight = mask[m,3]
    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if fwsre_cur > mask_lo
                if fwsre_cur > mask_hi
                    projection[p] += mask_weight * (mask_hi - mask_lo) / (fwsre_cur - fwsle_cur)
                    m += 1
                    if m<=size(mask,1)
                        mask_lo = mask[m,1] * shift_factor
                        mask_hi = mask[m,2] * shift_factor
                        mask_weight =  mask[m,3]
                    else
                        break
                    end
                else
                    projection[p] += mask_weight * (fwsre_cur - mask_lo) / (fwsre_cur - fwsle_cur)
                    on_mask = true
                    p+=1
                    fwsle_cur = fwsle[p]
                    fwsre_cur = fwsle[p+1]
                end
            else
                p+=1
                fwsle_cur = fwsle[p]
                fwsre_cur = fwsle[p+1]
            end
        else
            if fwsre_cur > mask_hi
                projection[p] += mask_weight * (mask_hi - fwsle_cur) / (fwsre_cur - fwsle_cur)
                on_mask = false
                m += 1
                if m<=size(mask,1)
                    mask_lo = mask[m,1] * shift_factor
                    mask_hi = mask[m,2] * shift_factor
                    mask_weight = mask[m,3]
                else
                    break
                end
            else
                projection[p] += mask_weight
                p += 1
                fwsle_cur = fwsle[p]
                fwsre_cur = fwsle[p+1]
            end
        end
        if p>length(projection) break end
    end
    return projection
end

function project_mask_compare(fws::A1, mask::A2; shift_factor::Real=one(T1)) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    # find bin edges
    fwsle, fwsre = find_bin_edges(fws)

    # allocate memory for mask projection
    nm = size(mask, 1)
    nx = size(fwsle, 1)
    ny = size(fwsle, 2)
    projection = zeros(nx, ny)
    #flush(stdout);  println("size(projection) = ", size(projection), " size(fws) = ", size(fws), " size(mask) = ", size(mask), " shift_factr = ", shift_factor)

    # shift the mask
    mask_shifted = zeros(size(mask))
    mask_shifted[:,1] = view(mask,:,1) .* shift_factor
    mask_shifted[:,2] = view(mask,:,2) .* shift_factor
    mask_shifted[:,3] = view(mask,:,3)
    mask = mask_shifted

    # set indexing variables
    p = 1
    m = 1
    on_mask = false

    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if fwsre[p] > mask[m,1]
                if fwsre[p] > mask[m,2]
                    projection[p] += mask[m,3] * (mask[m,2] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    m += 1
                else
                    projection[p] += mask[m,3] * (fwsre[p] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    on_mask = true
                    p+=1
                end
            else
                p+=1
            end
        else
            if fwsre[p] > mask[m,2]
                projection[p] += mask[m,3] * (mask[m,2] - fwsle[p]) / (fwsre[p] - fwsle[p])
                on_mask = false
                m += 1
            else
                projection[p] += mask[m,3]
                p += 1
            end
        end
        if p>length(projection) break end
    end
    return projection
end

"""
    ccf_1D(ws, sp, mask1)

Compute the cross correlation function of a spectrum and a tophat mask given
1D arrays of wavelengths and flux (ws and sp) and a 2D mask with one entry
per line. One mask entry consists of the left and right end of tophats and a
weight.
"""
function ccf_1D(ws::A1, sp::A2, mask1::A3
                ; RVguess::Real=zero(T1),
                wccf::Real=CCF_width,
                res_factor::Real=one(T1)) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,2} }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(ws) == 1
    @assert ndims(sp) == 1

    ccfd = calc_ccf_Δv_grid(RVguess=RVguess, wccf=wccf, res_factor=res_factor)
    ccf = zeros(size(ccfd))

    # calculate shift factors ahead of time
    shift_factors = calc_doppler_factor.(ccfd)
    mask_projections = zeros(length(ws),1)

    # loop through each velocity, project the mask and compute ccf
    for i in 1:length(ccf)
        # project shifted mask on wavelength domain
        #mask_projections .= project_mask_compare(ws, mask1, shift_factor=shift_factors[i])
        project_mask_opt(ws, mask1, shift_factor=shift_factors[i], projection=mask_projections)

        # compute the ccf value at the current velocity shift
        ccf[i] = sum(sp .* mask_projections)
    end
    return ccf
end

function ccf_1D(ws::A1, sp::A2, mask1::A3
            ; RVguess::Real=zero(T1), wccf::Real=CCF_width, res_factor::Real=one(T1)) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2} }
    # get number of dims
    ntimes = size(sp, 2)
    nccfa = length(calc_ccf_Δv_grid(RVguess=RVguess, wccf=wccf, res_factor=res_factor))

    # allocate memory & loop over rest
    ccfa = zeros(nccfa, ntimes);
    Threads.@threads for i in 1:ntimes
        ccfa[:,i] .= ccf_1D(ws, view(sp,:,i), mask1, RVguess=RVguess, wccf=wccf, res_factor=res_factor)
    end
    return ccfa
end

"""
    calc_ccf_Δv_grid(; opt params)

Compute the cross correlation function of a spectrum and a tophat mask given
1D arrays of wavelengths and flux (ws and sp) and a 2D mask with one entry
per line. One mask entry consists of the left and right end of tophats and a
weight.
"""
function calc_ccf_Δv_grid(; RVguess::Real=0.0, wccf::Real=CCF_width, res_factor::Real=1.0)
    # set arbitrary params that AbstractFloatfect RV measurement result
    dv = dV / res_factor
    nccf = ceil(Int, wccf/dv)
    # calculate ccf domain and allocate memory for ccf
    ccfd = range(RVguess-nccf*dv, RVguess+nccf*dv, length=2*nccf+1)
end

=#

end # module
