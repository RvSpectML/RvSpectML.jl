"""
    Code for a Gaussian mask shape for use with CCF
Author: Eric Ford
Created: August 2019
"""

default_σ = 5000

"""   GaussianCCFMask
A truncated Gaussian mask with two parameters, it's stdandard deviation and where to truncate it, both as a velocity in m/s.
Mask weights are stored separately in a line list.
"""
struct GaussianCCFMask <: AbstractCCFMaskShape
    σ_sqrt2::Float64
    half_width_truncation::Float64
    normalization::Float64

    """ GaussianCCFMask( σ ; half_truncation_width_in_σ=2 ) """
    function GaussianCCFMask(σ::Real, w::Real=2 )
        @assert 0 < σ <= 300000   # 300km/s is arbitrary choice for an upper limit
        @assert 0 < w <= 4
        norm = 1.0/(sqrt(2π)*σ*erf(w/(sqrt(2.0))))
        new(σ*sqrt(2.0),σ*w,norm)
    end

end

""" GaussianCCFMask( inst ; scale_factor ) """
function GaussianCCFMask(inst::AbstractInstrument; scale_factor::Real = 1)
    #=
    σ_sqrt2 = scale_factor * default_σ * sqrt(2.0)
    w = scale_factor * RvSpectML.default_ccf_mask_v_width(inst)
    norm = 1.0/(σ_sqrt2*sqrt(π)*erf(w/σ_sqrt2))
    GaussianCCFMask(σ_sqrt2,w/2) # ,norm)
    =#
    σ = scale_factor * default_gaussian_ccf_σ
    w = scale_factor * RvSpectML.default_ccf_mask_v_width(inst) / σ
    GaussianCCFMask(σ,w)
end

λ_min(m::GaussianCCFMask,λ::Real) = λ/calc_doppler_factor(m.half_width_truncation)
λ_max(m::GaussianCCFMask,λ::Real) = λ*calc_doppler_factor(m.half_width_truncation)

""" Functor for returning PSF for Δv <= half_width.  """
function (m::GaussianCCFMask)(Δv::Real)
    if abs2(Δv) > abs2(m.half_width_truncation)
        return zero(Δv)
    else
        return m.normalization*exp(-(Δv/m.σ_sqrt2)^2)
    end
end
