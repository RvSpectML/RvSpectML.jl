"""
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Optimizations by Eric Ford
"""

"""
Module for estimating the radial velocity based on the CCF
"""
module RVFromCCF

export measure_rv_from_ccf, measure_rvs_from_ccf
import ..RvSpectML: searchsortednearest

# packages
using LsqFit
using LinearAlgebra
import Polynomials

"""
    measure_rv(wavs, spec, mask; fit_type="gaussian")

    TODO: Need to update for CCF refactoring.

Compute the cross correlation function and fit to calculate a velocity
using the specified method. Valid arguments for fit_type include
"gaussian", "quadratic", and "centroid".

"""
function measure_rv(wavs::A1, spec::A2, mask::A3; fit_type::String="gaussian") where  {
                    T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,2}  }
    vels, ccf = ccf_1D(wavs, spec, mask, res_factor=1.0)

    # do the fit for velocity
    @assert fit_type in ["quadratic", "gaussian", "centroid"]
    if fit_type == "quadratic"
        return rv_from_ccf_quadratic(vels, ccf)
    elseif fit_type == "gaussian"
        return rv_from_ccf_gaussian(vels, ccf)
    elseif fit_type == "centroid"
        return rv_from_ccf_centroid(vels, ccf)
    else
        return nothing
    end
end

# TODO: Need to update for CCF refactoring.
function measure_rv(wavs::A1, spec::A2,
                    mask::A3; fit_type::String="gaussian") where {
                    T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    RVs = zeros(size(spec,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv(wavs, spec[:,i], mask, fit_type=fit_type)
    end
    return RVs
end


"""
    measure_rv_from_ccf(vels, ccf; fit_type="gaussian")

Fit a gaussian to the CCF to calculate a velocity
using the specified method. Valid arguments for fit_type include
"gaussian", "quadratic", and "centroid".
"""
function measure_rv_from_ccf(vels::A1, ccf::A2,
            ; fit_type::String="gaussian") where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }

    # do the fit for velocity
    @assert fit_type in ["quadratic", "gaussian", "centroid", "bestfit"]
    if fit_type == "quadratic"
        return rv_from_ccf_quadratic(vels, ccf)
    elseif fit_type == "gaussian"
        return rv_from_ccf_gaussian(vels, ccf)
    elseif fit_type == "centroid"
        return rv_from_ccf_centroid(vels, ccf)
    elseif fit_type == "bestfit"
        return rv_from_ccf_bestfit(vels, ccf)
    else
        return nothing
    end
end

function measure_rv_from_ccf(vels::A1, ccf::A2; fit_type::String="gaussian")  where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2} }
    RVs = zeros(size(ccf,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv_from_ccf(vels, ccf[:,i], fit_type=fit_type)
    end
    return RVs
end



"""
    measure_rvs_from_ccf(vels, ccf; fit_type="gaussian")

At each time, compute the cross correlation function and fit to calculate a velocity
using the specified method. Valid arguments for fit_type include
"gaussian", "quadratic", and "centroid".
"""
function measure_rvs_from_ccf(vels::A1, ccfs::A2,
            ; fit_type::String="gaussian") where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
 rvs = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(vels,ccfs[:,i],fit_type=fit_type) for i in 1:size(ccfs,2) ]
 return rvs
end


"""

"""
function prelim_fit(vels::A1, ccf::A2 ) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    dep = maximum(ccf) - minimum(ccf)
    val = (dep/2.0) + minimum(ccf)
    ind1 = findfirst(ccf .<= val)
    ind2 = findlast(ccf .<= val)
    return vels[ind2] - vels[ind1]
end

function find_idx_at_and_around_minimum(vels::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # do a prelim fit to get the width
    wid = prelim_fit(vels, ccf)

    # find the min and fit only that
    amin = argmin(ccf)
    lend = vels[amin] - 0.5 * wid
    rend = vels[amin] + 0.5 * wid

    # get the indices
    lind = searchsortednearest(vels[1:amin], lend)
    rind = amin + searchsortednearest(vels[amin+1:end], rend)
    inds = lind:rind

    return (amin, inds)
end

"""
Report RV based on fitting quadratic near minimum of CCF.
"""
function rv_from_ccf_quadratic(vels::A1, ccf::A2) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only that
    amin, inds = find_idx_at_and_around_minimum(vels, ccf)

    # do the polyfit
    pfit = Polynomials.fit(vels[inds], ccf[inds], 2)
    @assert length(Polynomials.coeffs(pfit)) >= 3   # just in case fails to fit a quadratic

    # get center from coeffs
    c, b, a = Polynomials.coeffs(pfit)
    return -b/(2*a)
end

@. gaussian_prof(x, p) = p[4] + p[3] * exp(-(x-p[1])^2.0/(2.0 * p[2]^2.0))


"""
Report RV base don fit Gaussian near minimum of CCF.
"""
function rv_from_ccf_gaussian(vels::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only that
    amin, inds = find_idx_at_and_around_minimum(vels, ccf)

    # make initial guess parameters
    μ = vels[amin]
    σ = 500                          # TODO: Make sure this is robust.  Probably need parameter
    amp = minimum(ccf) - maximum(ccf)
    y0 = maximum(ccf)
    p0 = [μ, σ, amp, y0]

    # fit and return the mean of the distribution
    fit = curve_fit(gaussian_prof, vels[inds], ccf[inds], p0)
    return coef(fit)[1]
end

"""
Report RV based on centroid near minimum of CCF.
WARNING:  Something's fishy with these.
"""
function rv_from_ccf_centroid(vels::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and centroid only that
    amin, inds = find_idx_at_and_around_minimum(vels, ccf)

    return -dot(ccf[inds], vels[inds]) / sum(ccf[inds])
end

"""
Report RV based on minimum of CCF
WARNING:  This is highly pixelated and should be used for special purposes (e.g., initial guesses).
"""
function rv_from_ccf_bestfit(vels::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    vels[findmin(ccf)[2]]
end

end # module
