"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Optimized by Eric Ford
"""

# packages
#using LsqFit
#using LinearAlgebra
#import Polynomials

module CCFTophat

export calc_ccf_Δv_grid, ccf_1D
export calc_ccf_chunk, calc_ccf_chunklist, calc_ccf_chunklist_timeseries

import ..RvSpectML: AbstractChuckOfSpectrum, AbstractChunkList, AbstractChunkListTimeseries

# physical constants
const c_ms = 2.99782458e8    # m/s
const c_kms = 2.99782458e5   # km/s

# set some CCF parameters
# TODO: Put into a struct
const CCF_width = 15e3       # m/s window
const dV = 250.0             # size of velocity step
const mask_width = 410.0     # m/s (width of mask entries)

function calc_ccf_chunk(chunk::AbstractChuckOfSpectrum, mask )
  res = ccf_1D(chunk.λ,chunk.flux, mask)
end

function calc_ccf_chunklist(chunk_list::AbstractChunkList, mask )
  mapreduce(chunk->calc_ccf_chunk(chunk, mask), +, chunk_list.data)
end

function calc_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries, mask )
  mapreduce(cl->calc_ccf_chunklist(cl, mask),hcat,clt.chunk_list)
end

function calc_order_ccfs_chunklist(chunk_list::AbstractChunkList, mask )
  mapreduce(chunk->calc_ccf_chunk(chunk, mask), hcat, chunk_list.data)
end

function calc_order_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries, mask )
  nΔvs = length(calc_ccf_Δv_grid())
  norders = length(clt.chunk_list[1].data)
  nobs =  length(clt.chunk_list)
  order_ccfs = zeros(nΔvs, norders, nobs)
  for i in 1:nobs
      order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], mask)
  end
  return order_ccfs
end

"""
   find_bin_edges(pixel_centers::AbstractArray)

TODO: Reduce allocations to one array
"""
function find_bin_edges(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = (fws[2:end] + fws[1:end-1]) ./ 2.0
    le = 2.0 * fws[1] - fwse[1]
    re = 2.0 * fws[end] - fwse[end]
    fwsle = vcat(le, fwse)
    fwsre = vcat(fwse, re)
    return fwsle, fwsre
end

function find_bin_edges_opt(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = Array{eltype(fws),1}(undef,length(fws)+2)
    fwse[2:end-2] = (fws[2:end] + fws[1:end-1]) ./ 2.0
    fwse[1] = 2.0 * fws[1] - fwse[2]
    fwse[end-1] = 2.0 * fws[end] - fwse[end-2]
    fwse[end] = zero(eltype(A))
    return fwse
end

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
    doppler_factor(vel)

Compute the longitudinal relativistic doppler factor given a velocity
in meters per second.
"""
function doppler_factor(vel::Real)
    num = one(vel) + vel/c_ms
    den = one(vel) - vel/c_ms
    return sqrt(num/den)
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
    shift_factors = doppler_factor.(ccfd)
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
    for i in 1:ntimes
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


end # module
