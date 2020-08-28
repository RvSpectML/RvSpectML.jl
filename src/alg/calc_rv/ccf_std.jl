"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
"""

module RVFromCCF

export measure_rv_from_ccf
import ..RvSpectML: searchsortednearest

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
function measure_rv(wavs::A1, spec::A2, mask::A3; fit_type::String="gaussian") where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,2}  }
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

function measure_rv(wavs::A1, spec::A2,
                    mask::A3; fit_type::String="gaussian") where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,2}  }
    RVs = zeros(size(spec,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv(wavs, spec[:,i], mask, fit_type=fit_type)
    end
    return RVs
end


"""
    measure_rv(wavs, spec, mask; fit_type="gaussian")

Compute the cross correlation function and fit to calculate a velocity
using the specified method. Valid arguments for fit_type include
"gaussian", "quadratic", and "centroid".
"""
function measure_rv_from_ccf(ccfd::A1, ccf::A2,
            ; fit_type::String="gaussian") where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }

    # do the fit for velocity
    @assert fit_type in ["quadratic", "gaussian", "centroid", "bestfit"]
    if fit_type == "quadratic"
        return rv_from_ccf_quadratic(ccfd, ccf)
    elseif fit_type == "gaussian"
        return rv_from_ccf_gaussian(ccfd, ccf)
    elseif fit_type == "centroid"
        return rv_from_ccf_centroid(ccfd, ccf)
    elseif fit_type == "bestfit"
        return rv_from_ccf_bestfit(ccfd, ccf)
    else
        return nothing
    end
end

function measure_rv_from_ccf(ccfd::A1, ccf::A2, mask::A3; fit_type::String="gaussian")  where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,2} }
    RVs = zeros(size(ccf,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv_from_ccf(ccfd, ccf[:,i], fit_type=fit_type)
    end
    return RVs
end


"""

"""
function prelim_fit(ccfd::A1, ccf::A2 ) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    dep = maximum(ccf) - minimum(ccf)
    val = (dep/2.0) + minimum(ccf)
    ind1 = findfirst(ccf .<= val)
    ind2 = findlast(ccf .<= val)
    return ccfd[ind2] - ccfd[ind1]
end

"""

"""
function rv_from_ccf_quadratic(ccfd::A1, ccf::A2) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
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

@. gaussian_prof(x, p) = p[4] + p[3] * exp(-(x-p[1])^2.0/(2.0 * p[2]^2.0))

"""

"""
function rv_from_ccf_gaussian(ccfd::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
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
function rv_from_ccf_centroid(ccfd::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    return -dot(ccf, ccfd) / sum(ccf)
end

"""

"""
function rv_from_ccf_bestfit(ccfd::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    ccfd[findmin(ccf)[2]]
end

end # module
