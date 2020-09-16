"""
    Code for a GaussianMixture mask shape for use with CCF
Author: Eric Ford
Created: September 2020
"""

using StaticArrays

"""   GaussianMixtureCCFMask
A mixture of truncated GaussianMixtures mask with three parameters per mixture component: weight,  standard deviation and offset (in m/s)
and one parameter for the velocity width at which the whole thing is truncated.
Mask weights are stored separately in a line list.
TODO: Test before using
"""
struct GaussianMixtureMixtureCCFMask{NMix::Integer} <: AbstractCCFMaskShape
    weight::SVector{NMix,Float64}
    σ_sqrt2::SVector{NMix,Float64}
    v_offset::SVector{NMix,Float64}
    truncation_Δv::Float64
    normalization::SVector{NMix,Float64}
    cdf_normalization::SVector{NMix,Float64}


    """ GaussianMixtureCCFMask( σ ; half_truncation_width_in_σ=2 ) """
    function GaussianMixtureCCFMask(weight::AbstractVector{Float64}, σ::AbstractVector{Float64}, truncation_Δv::Real; v_offset::AbstractVector{Float64}=zeros(length(weight)) )
        n = legnth(weight)
        @assert 1 <= n <= 10   # Arbitrary upper limit
        @assert length(σ) == n
        @assert length(v_offset) == n
        @assert all(0 .<= weight .<= one(eltype(weight)))   # 300km/s is arbitrary choice for an upper limit
        @assert all(0 .< σ .<= 300000)   # 300km/s is arbitrary choice for an upper limit
        @assert all(-1000 .<= v_offset .<= 1000)   # 1km/s is arbitrary choice for an upper limit
        @assert 0 < truncation_Δv <= 300000

        pdf_norm = zeros(n)
        cdf_norm = zeros(n)
        for i in 1:n
            integral = ( erf((truncation_Δv-v_offset[i])/(sqrt(2.0)*σ[i])) - erf((-truncation_Δv-v_offset[i])/(sqrt(2.0*σ[i]))) )
            pdf_norm[i] = weight[i]/(sqrt(2π)*σ[i]*integral)
            cdf_norm[i] = weight[i]/integral
        end
        new(weight, σ.*sqrt(2.0),v_offset, truncation_Δv,pdf_norm,cdf_norm)
    end

end

""" GaussianMixtureCCFMask( inst ; scale_factor )
For now just makes a single Gaussian, so for testing purposes only.
"""
function GaussianMixtureCCFMask(inst::AbstractInstrument; σ_scale_factor::Real = 1, truncation_Δv::Real = default_gaussian_mixture_ccf_truncation_Δv )
    σ = σ_scale_factor * RvSpectML.default_ccf_mask_v_width(inst)
    GaussianMixtureCCFMask([1.0], [σ],truncation_Δv,v_offset=[0.])
end

λ_min(m::GaussianMixtureCCFMask,λ::Real) = λ/calc_doppler_factor(m.truncation_Δv)
λ_max(m::GaussianMixtureCCFMask,λ::Real) = λ*calc_doppler_factor(m.truncation_Δv)

function cdf(m::GaussianMixtureCCFMask, v_lo::Real,v_hi::Real, i::Integer)
    m.cdf_normalization[i]*( erf((v_hi-v_offset[i])/σ_sqrt2[i]) - erf((-v_lo-v_offset[i])/σ_sqrt2[i]) )
end

function pdf(m::GaussianMixtureCCFMask, Δv::Real, i::Integer)
    m.normalization[i]*exp(-(Δv-m.v_offset[i]/m.σ_sqrt2[i])^2)
end

function integrate(m::GaussianMixtureCCFMask, v_lo::Real,v_hi::Real)
    @assert -m.truncation_Δv <= v_lo < v_hi <= m.truncation_Δv
    #quadgk(m, v_lo, v_hi)[1]
    mapreduce(i->cdf(m,v_lo,v_hi,i),+,1:length(m.weights))
end

""" Functor for returning PSF for Δv <= half_width.  """
function (m::GaussianMixtureCCFMask)(Δv::Real)
    if abs2(Δv) > abs2(m.truncation_Δv)
        return zero(Δv)
    else
        return mapreduce(i->pdf(m,Δv,i),+,1:length(m.weights))
    end
end
