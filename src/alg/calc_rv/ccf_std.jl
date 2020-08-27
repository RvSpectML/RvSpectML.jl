"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
"""
# packages
using LsqFit
using LinearAlgebra
import Polynomials

"""
    measure_rv(wavs, spec, mask; fit_type="gaussian")

Compute the cross correlation function and fit to calculate a velocity
using the specified method. Valid arguments for fit_type include
"gaussian", "quadratic", and "centroid".
"""
function measure_rv(wavs::AbstractArray{T,1}, spec::AbstractArray{T,1},
                    mask::AbstractArray{T,2}; fit_type::String="gaussian") where T<:AbstractFloat
    ccfd, ccf = ccf_1D(wavs, spec, mask, res_factor=1.0)

    # do the fit for velocity
    @assert fit_type in ["quadratic", "gaussian", "centroid"]
    if fit_type == "quadratic"
        return rv_from_ccf_quadratic(ccfd, ccf)
    elseif fit_type == "gaussian"
        return rv_from_ccf_gaussian(ccfd, ccf)
    elseif fit_type == "centroid"
        return rv_from_ccf_centroid(ccfd, ccf)
    else
        return nothing
    end
end

function measure_rv(wavs::AbstractArray{T,1}, spec::AbstractArray{T,2},
                    mask::AbstractArray{T,2}; fit_type::String="gaussian") where T<:AbstractFloat
    RVs = zeros(size(spec,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv(wavs, spec[:,i], mask, fit_type=fit_type)
    end
    return RVs
end

"""

"""
function prelim_fit(ccfd::AbstractArray{T,1}, ccf::AbstractArray{T,1}) where T<:AbstractFloat
    dep = maximum(ccf) - minimum(ccf)
    val = (dep/2.0) + minimum(ccf)
    ind1 = findfirst(ccf .<= val)
    ind2 = findlast(ccf .<= val)
    return ccfd[ind2] - ccfd[ind1]
end

"""

"""
function rv_from_ccf_quadratic(ccfd::AbstractArray{T,1}, ccf::AbstractArray{T,1}) where T<:AbstractFloat
    # do a prelim fit to get the width
    wid = prelim_fit(ccfd, ccf)

    # find the min and fit only that
    amin = argmin(ccf)
    lend = ccfd[amin] - 0.5 * wid
    rend = ccfd[amin] + 0.5 * wid

    # get the indices
    lind = searchsortednearest(ccfd[1:amin], lend)
    rind = amin + searchsortednearest(ccfd[amin+1:end], rend)
    inds = lind:rind

    # do the polyfit
    pfit = Polynomials.fit(ccfd[inds], ccf[inds], 2)

    # get center from coeffs
    c, b, a = Polynomials.coeffs(pfit)
    return -b/(2*a)
end

"""

"""
function rv_from_ccf_gaussian(ccfd::AbstractArray{T,1}, ccf::AbstractArray{T,1}) where T<:AbstractFloat
    # do a prelim fit to get the width
    wid = prelim_fit(ccfd, ccf)

    # find the min and fit only that
    amin = argmin(ccf)
    lend = ccfd[amin] - 0.5 * wid
    rend = ccfd[amin] + 0.5 * wid

    # get the indices
    lind = searchsortednearest(ccfd[1:amin], lend)
    rind = amin + searchsortednearest(ccfd[amin+1:end], rend)
    inds = lind:rind

    # make initial guess parameters
    μ = ccfd[argmin(ccf)]
    σ = 500
    amp = minimum(ccf) - maximum(ccf)
    y0 = maximum(ccf)
    p0 = [μ, σ, amp, y0]

    # fit and return the mean of the distribution
    fit = curve_fit(gaussian_prof, ccfd[inds], ccf[inds], p0)
    return coef(fit)[1]
end

"""

"""
function rv_from_ccf_centroid(ccfd::AbstractArray{T,1}, ccf::AbstractArray{T,1}) where T<:AbstractFloat
    return -dot(ccf, ccfd) / sum(ccf)
end
