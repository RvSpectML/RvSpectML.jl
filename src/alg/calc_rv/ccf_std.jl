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
export AbstractMeasureRvFromCCF, MeasureRvFromCCFCentroid, MeasureRvFromCCFQuadratic, MeasureRvFromCCFGaussian, MeasureRvFromCCFBestFit

import ..RvSpectML: searchsortednearest

# packages
using LsqFit
using LinearAlgebra
import Polynomials

#=
"""
    measure_rv(wavs, spec, mask; fit_type=:gaussian)

    TODO: Need to update for CCF refactoring.

Compute the cross correlation function and fit to calculate a velocity
using the specified method. Valid arguments for fit_type include
:gaussian, :quadratic", and :centroid.

"""
function measure_rv(wavs::A1, spec::A2, mask::A3; fit_type::Symbol=:gaussian) where  {
                    T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,2}  }
    vels, ccf = ccf_1D(wavs, spec, mask, res_factor=1.0)

    # do the fit for velocity
    @assert fit_type in [:quadratic, :gaussian, :centroid]
    if fit_type == :quadratic
        return rv_from_ccf_quadratic(vels, ccf)
    elseif fit_type == :gaussian
        return rv_from_ccf_gaussian(vels, ccf)
    elseif fit_type == :centroid
        return rv_from_ccf_centroid(vels, ccf)
    else
        return nothing
    end
end

function measure_rv(wavs::A1, spec::A2, mask::A3; fit_type::String="gaussian") where  {
                    T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    @warn "Please update call to measure_rv to use symbol instead of string for fit_type."
    if fit_type == "quadratic"
        return measure_rv(wavs,spec,mask,fit_type=:quadratic)
    elseif fit_type == "gaussian"
        return measure_rv(wavs,spec,mask,fit_type=:gaussian)
    elseif fit_type == "centroid"
        return measure_rv(wavs,spec,mask,fit_type=:centroid)
    else
        return nothing
    end
end

# TODO: Need to update for CCF refactoring.
function measure_rv(wavs::A1, spec::A2,
                    mask::A3; fit_type::Symbol=:gaussian) where {
                    T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, T3<:Real, A3<:AbstractArray{T3,2}  }
    RVs = zeros(size(spec,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv(wavs, spec[:,i], mask, fit_type=fit_type)
    end
    return RVs
end
=#

default_measure_width_at_frac_depth = 0.5
default_frac_of_width_to_fit = 0.5
default_init_guess_ccf_σ = 500.0

""" Abstract type for functors to estimate the raidal velocitiy from a CCF and its velocity grid.  """
abstract type AbstractMeasureRvFromCCF end

"""  Functor to estimate RV based on the centroid of the CCF.  """
struct MeasureRvFromCCFCentroid <: AbstractMeasureRvFromCCF
    frac_of_width_to_fit::Float64
    measure_width_at_frac_depth::Float64
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
* frac_of_width_to_fit: (0.5)
* measure_width_at_frac_depth: (0.5)
"""
function MeasureRvFromCCFCentroid(; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth )
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 2.0
    MeasureRvFromCCFCentroid(frac_of_width_to_fit,measure_width_at_frac_depth)
end

"""
Estimate RV based on centroid velocity of the CCF.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Only makes sense if the CCF has been variance normalized.  So don't use this yet.
"""
function (::AbstractMeasureRvFromCCF)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} } end

function (mrv::MeasureRvFromCCFCentroid)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and use only that part of the CCF for computing centroid
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

    weights = exp.(-0.5.*(ccf[inds].-ccf[amin]))
    v_centroid = sum(weights.*vels[inds]) / sum(weights)
    return v_centroid
end

"""  Functor to estimate RV based on fitting quadratic near minimum of CCF.
TODO: Revist the logic here and see if need to perform transformation first.
"""
struct MeasureRvFromCCFQuadratic <: AbstractMeasureRvFromCCF
    frac_of_width_to_fit::Float64
    measure_width_at_frac_depth::Float64
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
* frac_of_width_to_fit: (0.5)
* measure_width_at_frac_depth: (0.5)
"""
function MeasureRvFromCCFQuadratic(; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth )
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 2.0
    MeasureRvFromCCFQuadratic(frac_of_width_to_fit,measure_width_at_frac_depth)
end

function (mrv::MeasureRvFromCCFQuadratic)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only the part near the minimum of the CCF
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

    # do the polyfit
    pfit = Polynomials.fit(vels[inds], ccf[inds], 2)
    @assert length(Polynomials.coeffs(pfit)) >= 3   # just in case fails to fit a quadratic

    # get center from coeffs
    c, b, a = Polynomials.coeffs(pfit)
    v_at_min_of_quadratic = -b/(2*a)
    return v_at_min_of_quadratic
end

"""  Functor to estimate RV based on fitting a Gaussian quadratic near minimum of the CCF. """
struct MeasureRvFromCCFGaussian <: AbstractMeasureRvFromCCF
    frac_of_width_to_fit::Float64
    measure_width_at_frac_depth::Float64
    init_guess_ccf_σ::Float64
end

"""
Construct functor to estimate RV based on the CCF.
Optional Arguments:
* frac_of_width_to_fit: (0.5)
* measure_width_at_frac_depth: (0.5)
* init_guess_ccf_σ: (500m/s)
"""
function MeasureRvFromCCFGaussian(; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
                                    measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth,
                                    init_guess_ccf_σ::Real = default_init_guess_ccf_σ )
    @assert 0.25 <= measure_width_at_frac_depth <= 0.75
    @assert 0.1 <= frac_of_width_to_fit <= 2.0
    @assert 1 <= init_guess_ccf_σ <= 30e3
    MeasureRvFromCCFGaussian(frac_of_width_to_fit,measure_width_at_frac_depth,init_guess_ccf_σ)
end

@. gaussian_line_helper(x, p) = p[4] + p[3] * exp(-0.5*((x-p[1])/p[2])^2)

function (mrv::MeasureRvFromCCFGaussian)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
        # find the min and fit only the part near the minimum of the CCF
        amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=mrv.frac_of_width_to_fit, measure_width_at_frac_depth=mrv.measure_width_at_frac_depth)

        # make initial guess parameters
        μ = vels[amin]
        σ = mrv.init_guess_ccf_σ                          # TODO: Make sure this is robust.
        minccf, maxccf = extrema(ccf)
        amp = minccf - maxccf
        y0 = maxccf
        p0 = [μ, σ, amp, y0]

        # fit and return the mean of the distribution
        result = curve_fit(gaussian_line_helper, view(vels,inds), view(ccf,inds), p0)

        # center of line is first parameter to gaussian_line_helper
        rv = result.converged ?  coef(result)[1] : NaN
        return rv
end


"""  Functor to estimate RV based on velocity at minimum of CCF.
Warning:  Since a discrete grid of velocities are evaluated, this should only be used in limited situations (e.g., an initial guess). """
struct MeasureRvFromCCFBestFit <: AbstractMeasureRvFromCCF
end

function (mrv::MeasureRvFromCCFBestFit)(vels::A1, ccf::A2 ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    idx_min_ccf = findmin(ccf)[2]
    vels_min_ccf = vels[idx_min_ccf]
    return vels_min_ccf
end


"""
    measure_rv_from_ccf(vels, ccf; alg )
Return estimated RV based on the CCF using specified algorithm object.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* alg: Functor specifying how to measure the RV from the CCF.  Options include:
MeasureRvFromCCFGaussian (default), MeasureRvFromCCFQuadratic, MeasureRvFromCCFCentroid, and MeasureRvFromCCFBestFit.
"""
function measure_rv_from_ccf(vels::A1, ccf::A2, ; alg::AlgT=MeasureRvFromCCFGaussian() )  where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, AlgT<:AbstractMeasureRvFromCCF }
    return alg(vels,ccf)
end

"""   measure_rv_from_ccf(vels, ccf; alg )
At each time, return the estimated radial velocity based on the CCF using the specified algorithm.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* alg: Functor specifying how to measure the RV from the CCF.  Options include:
MeasureRvFromCCFGaussian (default), MeasureRvFromCCFQuadratic, MeasureRvFromCCFCentroid, and MeasureRvFromCCFBestFit.
"""
function measure_rvs_from_ccf(vels::A1, ccf::A2; alg::AlgT=MeasureRvFromCCFGaussian() )  where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, AlgT<:AbstractMeasureRvFromCCF }
    RVs = zeros(size(ccf,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv_from_ccf(vels, view(ccf,:,i), alg=alg)
    end
    return RVs
end


#=
"""
    measure_rvs_from_ccf(vels, ccf; alg)
At each time, compute the cross correlation function and fit to calculate a velocity using the specified method.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* fit_type:  What method to use for estimating the velocity.  Valid arguments for fit_type include :gaussian, :quadratic, and :centroid.  (:gaussian)
"""
function measure_rvs_from_ccf(vels::A1, ccfs::A2; alg::AlgT=MeasureRvFromCCFGaussian() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1},  AlgT<:AbstractMeasureRvFromCCF  }
   rvs = [ measure_rv_from_ccf(vels,view(ccfs,:,i),alg=alg) for i in 1:size(ccfs,2) ]
   return rvs
end
=#





#=
"""
    measure_rv_from_ccf(vels, ccf; fit_type=:gaussian)
Inputs:
    * vels: Array of velocites where CCF was evaluated.
    * ccf:  Array of values of CCF
Optional Arguements:
    * fit_type:  What method to use for estimating the velocity.  Valid arguments for fit_type include :gaussian, :quadratic, and :centroid.  (:gaussian)

Fit a gaussian to the CCF to calculate a velocity using the specified method.
"""
function measure_rv_from_ccf(vels::A1, ccf::A2; fit_type::Symbol=:gaussian) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }

    # do the fit for velocity
    @assert fit_type in [:quadratic, :gaussian, :centroid, :bestfit]
    if fit_type == :quadratic
        return rv_from_ccf_quadratic(vels, ccf)
    elseif fit_type == :gaussian
        return rv_from_ccf_gaussian(vels, ccf)
    elseif fit_type == :centroid
        return rv_from_ccf_centroid(vels, ccf)
    elseif fit_type == :bestfit
        return rv_from_ccf_bestfit(vels, ccf)
    else
        return nothing
    end
end

# Only kept around for compatability
function measure_rv_from_ccf(vels::A1, ccf::A2,
            ; fit_type::String) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    @warn "Please update call to measure_rv_from_ccf to use symbol instead of string for fit_type."
    # do the fit for velocity
    if fit_type == "quadratic"
        return rv_from_ccf_quadratic(vels, ccf, fit_type=:quadratic)
    elseif fit_type == "gaussian"
        return rv_from_ccf_gaussian(vels, ccf, fit_type=:gaussian)
    elseif fit_type == "centroid"
        return rv_from_ccf_centroid(vels, ccf, fit_type=:centroid)
    elseif fit_type == "bestfit"
        return rv_from_ccf_bestfit(vels, ccf, fit_type=:bestfit)
    else
        return nothing
    end
end

"""   measure_rv_from_ccf(vels, ccf; fit_type )
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* fit_type:  What method to use for estimating the velocity.  Valid arguments for fit_type include :gaussian, :quadratic, and :centroid.  (:gaussian)
"""
function measure_rv_from_ccf(vels::A1, ccf::A2; fit_type::Symbol=:gaussian)  where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} } #, T3<:Real, A3<:AbstractArray{T3,2} }
    RVs = zeros(size(ccf,2))
    for i in eachindex(RVs)
        RVs[i] = measure_rv_from_ccf(vels, ccf[:,i], fit_type=fit_type)
    end
    return RVs
end



"""
    measure_rvs_from_ccf(vels, ccf; fit_type=:gaussian)
At each time, compute the cross correlation function and fit to calculate a velocity using the specified method.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* fit_type:  What method to use for estimating the velocity.  Valid arguments for fit_type include :gaussian, :quadratic, and :centroid.  (:gaussian)
"""
function measure_rvs_from_ccf(vels::A1, ccfs::A2,
            ; fit_type::Symbol=:gaussian) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
   rvs = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(vels,ccfs[:,i],fit_type=fit_type) for i in 1:size(ccfs,2) ]
   return rvs
end


"""
   rv_from_ccf_quadratic(vels, ccf; frac_of_width_to_fit, measure_width_at_frac_depth )
Return RV based on fitting quadratic near minimum of CCF.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* frac_of_width_to_fit: How large of a range of velocities should be included in the fit (0.5)
* measure_width_at_frac_depth: What fractional CCF depth should be used for defining the width (0.5)

TODO: Revist the logic here and see if need to perform transformation first.
"""
function rv_from_ccf_quadratic(vels::A1, ccf::A2; frac_of_width_to_fit::Real = default_frac_of_width_to_fit, measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only that
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=frac_of_width_to_fit, measure_width_at_frac_depth=measure_width_at_frac_depth)

    # do the polyfit
    pfit = Polynomials.fit(vels[inds], ccf[inds], 2)
    @assert length(Polynomials.coeffs(pfit)) >= 3   # just in case fails to fit a quadratic

    # get center from coeffs
    c, b, a = Polynomials.coeffs(pfit)
    return -b/(2*a)
end





"""
Report RV based on fit Gaussian near minimum of CCF.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* frac_of_width_to_fit: How large of a range of velocities should be included in the fit (0.5)
* measure_width_at_frac_depth: What fractional CCF depth should be used for defining the width (0.5)
* init_guess_ccf_σ:  Initial guess for σ when performign itterative fit of Gaussian to CCF
"""
function rv_from_ccf_gaussian(vels::A1, ccf::A2; frac_of_width_to_fit::Real = default_frac_of_width_to_fit,
            measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth, init_guess_ccf_σ::Real = default_init_guess_ccf_σ ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and fit only that
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=frac_of_width_to_fit, measure_width_at_frac_depth=measure_width_at_frac_depth)

    # make initial guess parameters
    μ = vels[amin]
    σ = init_guess_ccf_σ                          # TODO: Make sure this is robust.
    amp = minimum(ccf) - maximum(ccf)
    y0 = maximum(ccf)
    p0 = [μ, σ, amp, y0]

    # fit and return the mean of the distribution
    fit = curve_fit(gaussian_prof, vels[inds], ccf[inds], p0)
    return coef(fit)[1]
end

"""   rv_from_ccf_centroid(vels, ccf )
Report RV based on centroid near minimum of CCF.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* frac_of_width_to_fit: How large of a range of velocities should be included in the fit (0.5)
* measure_width_at_frac_depth: What fractional CCF depth should be used for defining the width (0.5)
WARNING:  Something's fishy with these.
TODO: Revisit algorithm/logic
"""
function rv_from_ccf_centroid(vels::A1, ccf::A2; frac_of_width_to_fit::Real = default_frac_of_width_to_fit, measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth ) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # find the min and centroid only that
    amin, inds = find_idx_at_and_around_minimum(vels, ccf, frac_of_width_to_fit=frac_of_width_to_fit, measure_width_at_frac_depth=measure_width_at_frac_depth)

    weights = exp.(-0.5.*(ccf[inds].-ccf[amin]))
    v_centroid = sum(weights.*vels[inds]) / sum(weights)
    return v_centroid
end

"""   rv_from_ccf_bestfit(vels, ccf)
Report RV based on minimum of CCF
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
WARNING:  This is highly pixelated and should only be used for special purposes (e.g., initial guesses).
"""
function rv_from_ccf_bestfit(vels::A1, ccf::A2) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    vels[findmin(ccf)[2]]
end

=#


"""  prelim_fit(vels, ccf; measure_width_at_frac_depth = 0.5 )
Return rough estimate of a ccf (or line) full width at the specified fractional depth (i.e., fraction of total depth).
This is based on the velocities of the first/last pixel to have a value less than the target value, with no interpolation.
Assumes vels is sorted.  Depth is measured assuming to the maximum value of ccf represents the continuum.
Could be improved for noisy data and segments with a steep slope due to the blaze or neighboring lines.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* measure_width_at_frac_depth: What fractional CCF depth should be used for defining the width (0.5)
"""
function prelim_fit(vels::A1, ccf::A2; measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth ) where  {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    minccf, maxccf = extrema(ccf)
    depth = maxccf - minccf
    #depth = maximum(ccf) - minimum(ccf)
    #target_val = minimum(ccf) + measure_width_at_frac_depth * (1-depth)
    target_val = minccf + (1-measure_width_at_frac_depth) * depth
    ind1 = findfirst(ccf .<= target_val)
    ind2 = findlast(ccf .<= target_val)
    if isnothing(ind1) || isnothing(ind2)
        println("ccf = ",ccf)
        println("minccf= ", minccf, " maxccf= ", maxccf, " depth= ", depth, " measure_width_at_frac_depth= ", measure_width_at_frac_depth, " targetval= ",target_val, " ind1= ", ind1, " ind2= ", ind2)
    end
    return vels[ind2] - vels[ind1]
end

"""  find_idx_at_and_around_minimum(vels, ccf; frac_of_width, measure_width_at_frac_depth )
Return the a pair with the index of the lowest value of ccf and a range of pixels surrounding it.
The range is based on finding the pixels with velocities nearest to the
Assumes vels is sorted.
Inputs:
* vels: Array of velocites where CCF was evaluated.
* ccf:  Array of values of CCF
Optional Arguements:
* frac_of_width_to_fit: How large of a range of velocities should be included in the fit (0.5)
* measure_width_at_frac_depth: What fractional CCF depth should be used for defining the width (0.5)
"""
function find_idx_at_and_around_minimum(vels::A1, ccf::A2; frac_of_width_to_fit::Real = default_frac_of_width_to_fit, measure_width_at_frac_depth::Real = default_measure_width_at_frac_depth) where {T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1} }
    # do a prelim fit to get the width
    full_width = prelim_fit(vels, ccf, measure_width_at_frac_depth=measure_width_at_frac_depth)

    # find the min and fit only that
    amin = argmin(ccf)
    lend = vels[amin] - frac_of_width_to_fit * full_width
    rend = vels[amin] + frac_of_width_to_fit * full_width

    # get the indices
    lind = searchsortednearest(vels[1:amin], lend)
    rind = amin + searchsortednearest(vels[amin+1:end], rend)
    inds = lind:rind

    return (amin, inds)
end

end # module
