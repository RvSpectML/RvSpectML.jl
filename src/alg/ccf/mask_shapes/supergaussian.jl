"""
    Code for a super-Gaussian mask shape for use with CCF
Author: Eric Ford
Created: September 2019
"""

using QuadGK
#using StaticArrays

"""   SuperGaussianCCFMask
A truncated Gaussian mask with two parameters, it's stdandard deviation and where to truncate it, both as a velocity in m/s.
Mask weights are stored separately in a line list.

TODO: Repalce Gaussian with super-Gaussian
Warning:  Not implemented/tested yet.
"""
struct SuperGaussianCCFMask <: AbstractCCFMaskShape
    σ_sqrt2::Float64
    power::Float64
    half_width_truncation::Float64
    normalization::Float64

    """ SuperGaussianCCFMask( σ ; half_truncation_width_in_σ=2 ) """
    function SuperGaussianCCFMask(σ::Real, p::Real, w::Real=2 )
        @assert 0 < σ <= 300000   # 300km/s is arbitrary choice for an upper limit
        @assert 1 <= p <= 2
        @assert 0 < w <= 4
        function integrand(x)
            exp(-(0.5*(x/σ)^2)^p)
        end
        integral = quadgk(integrand, -w*σ/2, w*σ/2)[1]
        norm = 1.0/integral
        new(σ*sqrt(2.0),p,σ*w,norm)
    end

end

""" SuperGaussianCCFMask( inst ; scale_factor ) """
function SuperGaussianCCFMask(inst::InstT; power::Real = default_supergaussian_ccf_exponent, fwhm::Real = default_supergaussian_ccf_fwhm,
        σ_scale_factor::Real = 1, truncation_scale_factor::Real = default_supergaussian_ccf_truncation_scale_factor ) where { InstT<:AbstractInstrument }
    # TODO: Update default values
    σ = σ_scale_factor * fwhm/sqrt(8 * log(2)^(1/power))  # From Ryan Petersburg email 9/11/2020
    w = truncation_scale_factor
    SuperGaussianCCFMask(σ,power,w/2)
end

λ_min(m::SuperGaussianCCFMask,λ::Real) = λ/calc_doppler_factor(m.half_width_truncation)
λ_max(m::SuperGaussianCCFMask,λ::Real) = λ*calc_doppler_factor(m.half_width_truncation)

function integrate(m::SuperGaussianCCFMask, v_lo::Real,v_hi::Real)
    #quadgk(m, v_lo, v_hi)[1]
    quadgk(m, v_lo, v_hi, atol=1e-4)[1]
end

""" Functor for returning PSF for Δv <= half_width.  """
function (m::SuperGaussianCCFMask)(Δv::Real)
    if abs2(Δv) > abs2(m.half_width_truncation)
        return zero(Δv)
    else
        return m.normalization*exp(-((Δv/m.σ_sqrt2)^2)^m.power)
    end
end
