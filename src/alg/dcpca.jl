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
#using LinearAlgebra #?

import RvSpectMLBase: speed_of_light_mps
import RvSpectMLBase: AbstractChunkListTimeseries, num_times, set_rv_est!
import RvSpectML: make_grids_for_chunklist_timeseries, pack_shifted_chunk_list_timeseries_to_matrix

#const speed_of_light_mps = speed_of_light_mps # 299792458.0 # TODO: Update value

export doppler_constrained_pca, compute_spectra_perp_doppler_shift
export clean_rvs_dcpca, calc_sigma_pca_scores

function compute_spectra_perp_doppler_shift(spectra::AA, deriv::V1, rvs::V2) where {
            T1<:Real, AA<:AbstractArray{T1,2}, T2<:Real, V1<:AbstractVector{T2}, T3<:Real, V2<:AbstractVector{T3} }
   @assert size(spectra,1) == length(deriv)
   @assert size(spectra,2) == length(rvs)
   fm_perp = spectra .- rvs' .* deriv./speed_of_light_mps
end

function doppler_constrained_pca_orig(flux::AbstractArray{T1,2}, deriv::AbstractVector{T2}, rvs::AbstractVector{T3} ; num_basis_vectors_pca::Integer = 1)  where { T1<:Real, T2<:Real, T3<:Real  }#, times::AbstractVector, rvs::AbstractVector)
  rv_shift = rvs .- mean(rvs)
  fm_perp = compute_spectra_perp_doppler_shift(flux,deriv, rv_shift )
  idx_good_obs = 1:length(rvs)
  M = fit(PCA, fm_perp[:,idx_good_obs]; maxoutdim=num_basis_vectors_pca)
  pca_out = MultivariateStats.transform(M,fm_perp[:,idx_good_obs])
  return pca_out, M
end

function calc_sigma_pca_scores(flux::AbstractArray{T1,2}, var::AbstractArray{T2,1}, model::PCA{T3}, n::Integer) where { T1<:Real, T2<:Real, T3<:Real }
  println("# indim = ", indim(model), "  outdim = ", outdim(model), " size(flux) = ", size(flux))
  mean_scores = MultivariateStats.transform(model, flux)
  println("# size(mean_scores) = ", size(mean_scores))
  samples = zeros(size(mean_scores,1),size(mean_scores,2),n)
  for i in 1:n
    samples[:,:,i] .= MultivariateStats.transform(model, flux .+ var .* randn(size(flux))) - mean_scores
  end
  std_scores = std(samples,dims=3)
  std_scores = reshape(std_scores,size(std_scores,1),size(std_scores,2))'
  return std_scores
end

#=
function calc_sigma_pca_scores(flux::AbstractArray{T1,2}, var::AbstractArray{T2,2}, model::PCA{T3}, n::Integer) where { T1<:Real, T2<:Real, T3<:Real }
  mean_scores = MultivariateStats.transform(model, flux)
  samples = zeros(size(mean_scores,1),size(mean_scores,2),n)
  for i in 1:n
    samples[:,:,i] .= MultivariateStats.transform(model, flux .+ var .* randn(size(flux))) - mean_scores
  end
  std_scores = sqrt.(std(samples,dims=3))
  std_scores = reshape(std_scores,size(std_scores,1),size(std_scores,2))'
  return std_scores
end
=#
function clean_rvs_dcpca(clt::ACLT, rvs::AbstractVector{T1}; num_basis_vectors_pca::Integer = 1,
              num_samples_σ_scores::Integer = 100,
              idx_use_pca::Union{Integer,AbstractVector{T2}} = 1:num_times(clt),
              chunk_grids::AbstractVector{AR} = make_grids_for_chunklist_timeseries(clt) ) where {
              ACLT<:AbstractChunkListTimeseries, T1<:Real, T2<:Integer, AR<:AbstractRange }
  @assert num_times(clt) == length(rvs)
  @assert 1 <= num_basis_vectors_pca <= num_times(clt)
  set_rv_est!(clt, rvs )
  result = pack_shifted_chunk_list_timeseries_to_matrix(clt, chunk_grids,oversample_factor=2.0,remove_rv_est = true)
  M = fit(PCA, result.matrix.flux[:,idx_use_pca]; maxoutdim=num_basis_vectors_pca)

  pca_scores = MultivariateStats.transform(M, result.matrix.flux)'
  σ_pca_scores = calc_sigma_pca_scores(result.matrix.flux, result.mean_var, M, num_samples_σ_scores)

  sol = llsq(pca_scores, rvs)
  a, b = sol[1:end-1], sol[end]
  pred =  pca_scores * a  .+ b
  if num_samples_σ_scores > 1
    σ_pred = σ_pca_scores * a
  else
    σ_pred = fill(NaN,length(rvs))
  end
  rvs_dirty = pred.-mean(pred)
  rvs_clean = rvs .- rvs_dirty

  return (rvs_clean=rvs_clean, rvs_dirty=rvs_dirty, σ_rvs_dirty=σ_pred, scores=pca_scores, σ_scores=σ_pca_scores, model=M)
end

end # DCPCA
