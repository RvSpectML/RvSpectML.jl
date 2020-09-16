"""
    Code for a half-cosine mask shape for use with CCF
Author: Eric Ford
Created: August 2020
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
        new(w)
    end

end

""" CosCCFMask( inst  ) """
function CosCCFMask(inst::AbstractInstrument; scale_factor::Real = 1)
    w = scale_factor * RvSpectML.default_ccf_mask_v_width(inst)
    CosCCFMask(w/2)
end

λ_min(m::CosCCFMask,λ::Real) = λ/calc_doppler_factor(m.half_width)
λ_max(m::CosCCFMask,λ::Real) = λ*calc_doppler_factor(m.half_width)

function integrate(m::CosCCFMask, v_lo::Real,v_hi::Real)
    #quadgk(m, v_lo, v_hi)[1]
    0.5*(sin(0.5*π*v_hi/m.half_width)-sin(0.5*π*v_lo/m.half_width))
end

""" Functor for returning PSF for Δv <= half_width.  """
function (m::CosCCFMask)(Δv::Real)
    if abs2(Δv) > abs2(m.half_width)
        return zero(Δv)
    else
        return cos(0.5*π*Δv/m.half_width)*π/(4*m.half_width)
    end
end
