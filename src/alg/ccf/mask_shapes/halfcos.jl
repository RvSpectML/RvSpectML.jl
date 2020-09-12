"""
    Code for a half-cosine mask shape for use with CCF
Author: Eric Ford
Created: August 2019
"""

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
