"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
"""

# set some constants as globals
c_ms = 2.99782458e8
c_kms = c_ms/1000.0

# set CCF parameters as globals
CCF_width = 15e3  # m/s window
dV = 820.0 / 10.0   # one tenth of one pixel (which instrument?)
mask_width = 1500.0 # a few pixels (which instrument?)

function find_bin_edges(wavs::AA{T,1}) where T<:Real
    fwse = (fws[2:end] + fws[1:end-1]) ./ 2.0
    le = 2.0 * fws[1] - fwse[1]
    re = 2.0 * fws[end] - fwse[end]
    fwsle = vcat(le, fwse)
    fwsre = vcat(fwse, re)
    return fwsle, fwsre
end

function project_mask(fws::AA{T,1}, mask::AA{T,2}; shift_factor::T=1.0) where T<:Real
    # find bin edges
    fwsle, fwsre = find_bin_edges

    # allocate memory for mask projection
    nm = size(mask, 1)
    nx = size(fwsle, 1)
    ny = size(fwsle, 2)
    projection = zeros(nx, ny)

    # shift the mask
    mask_shifted = zeros(size(mask))
    mask_shifted[:,1] = mask[:,1] .* shift_factor
    mask_shifted[:,2] = mask[:,2] .* shift_factor
    mask_shifted[:,3] = mask[:,3]
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

function ccf_1D(ws::AA{T,1}, sp::AA{T,1}, mask1::AA{T,2};
                RVguess::T=0.0, wccf::T=CCF_width, res_factor::T=1.0) where T<:Real
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(ws) == 1
    @assert ndims(sp) == 1

    # set arbitrary params that affect RV measurement result
    dv = dV / res_factor
    nccf = ceil(Int, wccf/dv)

    # calculate ccf domain and allocate memory for ccf
    ccfd = range(RVguess-nccf*dv, RVguess+nccf*dv, length=2*nccf+1)
    ccf = zeros(size(ccfd))

    # calculate shift factors ahead of time
    shift_factors = sqrt.((1.0 .+ ccfd./c_ms)./(1.0 .- ccfd./c_ms))

    # loop through each velocity, project the mask and compute ccf
    for i in 1:length(ccf)
        # project shifted mask on wavelength domain
        mask_projections = project_mask(ws, mask1, shift_factor=shift_factors[i])

        # compute the ccf value at the current velocity shift
        ccf[i] = sum(sp .* mask_projections)
    end
    return ccfd, ccf
end

function ccf_1D(ws::AA{T,1}, sp::AA{T,2}, mask1::AA{T,2};
                RVguess::T=0.0, wccf::T=CCF_width, res_factor::T=1.0) where T<:Real
    # get number of dims
    ntimes = size(sp, 2)

    # do the first ccf
    ccfd, ccf = ccf_1D(ws, sp[:,1], mask1, RVguess=RVguess, wccf=wccf, res_factor=res_factor)

    # allocate memory & loop over rest
    ccfa = zeros(length(ccf), ntimes); ccfa[:,1] = ccf
    for i in 2:ntimes
        ccfd, ccftemp = ccf_1D(ws, sp[:,i], mask1, RVguess=RVguess, wccf=wccf, res_factor=res_factor)
        ccfa[:,i] = ccftemp
    end
    return ccfd, ccfa
end

@. gaussian(x, p) = p[4] + p[3] * exp(-(x-p[1])^2.0/(2.0 * p[2]^2.0))

function rv_from_ccf(ccfd::AA{T,1}, ccf::AA{T,1}; fit_type="gauss_fit") where T<:Real
    # make initial guess parameters
    μ = ccfd[argmin(ccf)]
    σ = 500
    amp = minimum(ccf) - maximum(ccf)
    y0 = maximum(ccf)
    p0 = [μ, σ, amp, y0]

    if fit_type == "gauss_fit"
        fit = curve_fit(gaussian, ccfd, ccf, p0)
    end
    return -coef(fit)[1]
end

function rv_from_ccf(ccfd::AA{T,1}, ccf::AA{T,2}; fit_type="gauss_fit") where T<:Real
    RVs = zeros(size(ccf,2))
    for i in 1:length(RVs)
        RVs[i] = rv_from_ccf(ccfd, ccf[:,i], fit_type=fit_type)
    end
    return RVs
end

# function measure_rv(wavs, spec, mask)
