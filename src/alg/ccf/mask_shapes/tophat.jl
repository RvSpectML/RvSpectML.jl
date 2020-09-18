"""
    Code for a top hat mask shape for use with CCF
Author: Eric Ford
Created: August 2020
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

#global num_error_msgs = 0

function integrate(m::TopHatCCFMask, v_lo::Real,v_hi::Real)
#=    if ! (-m.half_width*1.01 <= v_lo <= m.half_width*1.01) || !(-m.half_width*1.01 <= v_hi <= m.half_width*1.01)
        global num_error_msgs += 1
        if num_error_msgs < 100
        println("# v_lo = ", v_lo, "  v_hi = ", v_hi, "   half_width = ", m.half_width)
#        @assert -m.half_width <= v_hi <= m.half_width
#        @assert -m.half_width <= v_lo <= m.half_width
        end
    end
    =#
    #quadgk(m, v_lo, v_hi)[1]
    0.5*(v_hi-v_lo)/m.half_width
    #0.5*(v_hi-v_lo)/m.half_width * (1+m.half_width/RvSpectML.speed_of_light_mps)/(1+(m.half_width/RvSpectML.speed_of_light_mps)^2/2)
end

""" Functor for returning 1 for any Δv <= width.  """
function (m::TopHatCCFMask)(Δv::Real)
    # TODO: Change once implement generic CCF mask shapes (currently, I don't understand why the normalization differs for the two)
    #return abs2(Δv)<=abs2(m.half_width) ? 0.5/m.half_width : zero(Δv)  # Old for project_mask that assumes tophat
    return abs(Δv)<=abs(m.half_width) ? 0.5/m.half_width : zero(Δv)  # Old for project_mask that assumes tophat
    #return abs2(Δv)<=abs2(m.half_width) ? 1.0 : zero(Δv)    # New version for project_mask that works with a probability density
end
