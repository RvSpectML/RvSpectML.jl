"""
    Code for a top hat mask shape for use with CCF
Author: Eric Ford
Created: August 2019
"""

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
