"""
Module for performing a Doppler-constrained PCA analysis

For algorithm information, see:
  - Jones, Stenning, Ford et al. https://arxiv.org/abs/1711.01318
  - See also: Gilbertson, Ford, Jones & Stenning 2020 https://arxiv.org/abs/2009.01085

Author: Eric Ford
Date:   September 2020
"""
module DCPCA

using Statistics, MultivariateStats

#import .RvSpectML : speed_of_light_mps
const speed_of_light_mps = 299792458.0 # TODO: Update value

export doppler_constrained_pca, compute_spectra_perp_doppler_shift

function compute_spectra_perp_doppler_shift(spectra::AA, deriv::V1, rvs::V2) where {
            T1<:Real, AA<:AbstractArray{T1,2}, T2<:Real, V1<:AbstractVector{T2}, T3<:Real, V2<:AbstractVector{T3} }
   @assert size(spectra,1) == length(deriv)
   @assert size(spectra,2) == length(rvs)
   fm_perp = spectra .- rvs' .* deriv./speed_of_light_mps
end

function doppler_constrained_pca(flux::AbstractArray{T1,2}, deriv::AbstractVector{T2}, rvs::AbstractVector{T3} )  where { T1<:Real, T2<:Real, T3<:Real  }#, times::AbstractVector, rvs::AbstractVector)
  rv_shift = rvs .- mean(rvs)
  fm_perp = compute_spectra_perp_doppler_shift(flux,deriv, rv_shift )
  idx_good_obs = 1:length(rvs)
  M = fit(PCA, fm_perp[:,idx_good_obs]; maxoutdim=10)
  pca_out = MultivariateStats.transform(M,fm_perp[:,idx_good_obs])
  return pca_out, M
end

end # DCPCA
