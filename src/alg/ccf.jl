"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
"""

# packages
#using LsqFit
#using LinearAlgebra
#import Polynomials

# physical constants
const c_ms = 2.99782458e8    # m/s
const c_kms = 2.99782458e5   # km/s

# set some CCF parameters
const CCF_width = 15e3       # m/s window
const dV = 250.0             # size of velocity step
const mask_width = 410.0     # m/s (width of mask entries)

"""

"""
function find_bin_edges(fws::AbstractArray{T,1}) where T<:AbstractFloat
    fwse = (fws[2:end] + fws[1:end-1]) ./ 2.0
    le = 2.0 * fws[1] - fwse[1]
    re = 2.0 * fws[end] - fwse[end]
    fwsle = vcat(le, fwse)
    fwsre = vcat(fwse, re)
    return fwsle, fwsre
end

"""
    project_mask(fws, mask, shift_factor)

Compute the projection of the tophat mask onto the 1D array
of wavelengths (fws) at a given shift factor.
"""
function project_mask(fws::AbstractArray{T,1}, mask::AbstractArray{T,2};
                      shift_factor::T=1.0) where T<:AbstractFloat
    # find bin edges
    fwsle, fwsre = find_bin_edges(fws)

    # allocate memory for mask projection
    nm = size(mask, 1)
    nx = size(fwsle, 1)
    ny = size(fwsle, 2)
    projection = zeros(nx, ny)

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
    end
    return projection
end


"""
    doppler_factor(vel)

Compute the longitudinal relativistic doppler factor given a velocity
in meters per second.
"""
function doppler_factor(vel::T) where T<:AbstractFloat
    num = one(T) + vel/c_ms
    den = one(T) - vel/c_ms
    return sqrt(num/den)
end


"""
    ccf_1D(ws, sp, mask1)

Compute the cross correlation function of a spectrum and a tophat mask given
1D arrays of wavelengths and flux (ws and sp) and a 2D mask with one entry
per line. One mask entry consists of the left and right end of tophats and a
weight.
"""
function ccf_1D(ws::AbstractArray{T,1}, sp::AbstractArray{T,1},
                mask1::AbstractArray{T,2}; RVguess::T=0.0,
                wccf::T=CCF_width, res_factor::T=1.0) where T<:AbstractFloat
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(ws) == 1
    @assert ndims(sp) == 1

    # set arbitrary params that AbstractFloatfect RV measurement result
    dv = dV / res_factor
    nccf = ceil(Int, wccf/dv)

    # calculate ccf domain and allocate memory for ccf
    ccfd = range(RVguess-nccf*dv, RVguess+nccf*dv, length=2*nccf+1)
    ccf = zeros(size(ccfd))

    # calculate shift factors ahead of time
    shift_factors = doppler_factor.(ccfd)

    # loop through each velocity, project the mask and compute ccf
    for i in 1:length(ccf)
        # project shifted mask on wavelength domain
        mask_projections = project_mask(ws, mask1, shift_factor=shift_factors[i])

        # compute the ccf value at the current velocity shift
        ccf[i] = sum(sp .* mask_projections)
    end
    return ccfd, ccf
end

function ccf_1D(ws::AbstractArray{T,1}, sp::AbstractArray{T,2},
                mask1::AbstractArray{T,2}; RVguess::T=0.0,
                wccf::T=CCF_width, res_factor::T=1.0) where T<:AbstractFloat
    # get number of dims
    ntimes = size(sp, 2)

    # do the first ccf
    ccfd, ccf = ccf_1D(ws, view(sp,:,1), mask1, RVguess=RVguess, wccf=wccf, res_factor=res_factor)

    # allocate memory & loop over rest
    ccfa = zeros(length(ccf), ntimes); ccfa[:,1] = ccf
    for i in 2:ntimes
        ccfd, ccftemp = ccf_1D(ws, view(sp,:,i), mask1, RVguess=RVguess, wccf=wccf, res_factor=res_factor)
        ccfa[:,i] = ccftemp
    end
    return ccfd, ccfa
end
