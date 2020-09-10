"""
    Code to specifing mask shape for CCF
Author: Eric Ford
Created: August 2019
"""

"A struct implementing a specific mask shapes should be a subtype of AbstractCCFMaskShape."
abstract type AbstractCCFMaskShape  end

"""   TopHatCCFMask
The standard tophat mask with one parameter, it's half width as a velocity in m/s.
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
    w = scale_factor * RvSpectML.default_ccf_mask_v_width(inst)
    TopHatCCFMask(w/2)
end

λ_min(m::TopHatCCFMask,λ::Real) = λ/calc_doppler_factor(m.half_width)
λ_max(m::TopHatCCFMask,λ::Real) = λ*calc_doppler_factor(m.half_width)

""" Functor for returning 1 for any Δv <= width.  """
function (m::TopHatCCFMask)(Δv::Real)
    # TODO: Change once implement generic CCF mask shapes (currently, I don't understand why the normalization differs for the two)
    return abs2(Δv)<=abs2(m.half_width)*1.001 ? 0.5/m.half_width : zero(Δv)  # Old for project_mask that assumes tophat
    #return abs2(Δv)<=abs2(m.half_width) ? 1.0 : zero(Δv)    # New version for project_mask that works with a probability density
end

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
    σ = scale_factor * default_σ
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

"""   CosCCFMask
Cosine mask with one parameter, it's quarter period, i.e., where to truncate it, as a velocity in m/s.
Mask weights are stored separately in a line list.
"""
struct CosCCFMask <: AbstractCCFMaskShape
    half_width::Float64
    # normalization::Float64 #  unneeded since norm=pi/(2*half_width_truncation)

    """ CosCCFMask( full_width ) """
    function CosCCFMask( w::Real )
        @assert 100 <= w <= 300000   # 0.1 and 300km/s are arbitrary choices for an upper limit
        new(w/2)
    end

end

""" CosCCFMask( inst  ) """
function CosCCFMask(inst::AbstractInstrument; scale_factor::Real = 1)
    w = scale_factor * RvSpectML.default_ccf_mask_v_width(inst)
    CosCCFMask(w)
end

λ_min(m::CosCCFMask,λ::Real) = λ/calc_doppler_factor(m.half_width)
λ_max(m::CosCCFMask,λ::Real) = λ*calc_doppler_factor(m.half_width)

""" Functor for returning PSF for Δv <= half_width.  """
function (m::CosCCFMask)(Δv::Real)
    if abs2(Δv) > abs2(m.half_width)
        return zero(Δv)
    else
        return cos(0.5*π*Δv/m.half_width)*π/(2*m.half_width)
    end
end
