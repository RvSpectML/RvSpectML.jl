
"""
`calc_depth_and_expected_rv_precission(spectrum, pixels_index, order_index; smooth_factor)`
Calculate expected RV uncertainty for one portion of spectrum given by pixels and order indicies.
Assumes only photon noise with given variances.
Returns (depth, `exp_σ_rv`)
"""
function calc_depth_and_expected_rv_precission(spectrum::ST, pixels::AR, order::Integer; smooth_factor::Real = 2.0) where { ST<:AbstractSpectra2D, AR<:AbstractRange{Int64}, AA1<:AbstractArray{Int64,1} }
   logλ = log.(view(spectrum.λ,pixels,order))
   flux = view(spectrum.flux,pixels,order)
   var = view(spectrum.var,pixels,order)
   (f_mean, f_deriv) = RvSpectML.TemporalGPInterpolation.predict_mean_and_deriv(logλ, flux, logλ;sigmasq_obs=var, use_logx=false, use_logy=false, smooth_factor=smooth_factor)
   depth = 1.0 .- minimum(f_mean)/maximum(f_mean)
   exp_sigma_rv = RvSpectMLBase.speed_of_light_mps / sqrt( sum(f_deriv.^2 ./ var) )
   return (depth=depth, exp_σ_rv=exp_sigma_rv)
 end

""" `calc_formal_rv_precission(spectrum, chunklist; smooth_factor)`
Calculate expected RV uncertainty for one spectrum and corresponding chunklist,
Assumes only photon noise with given variances.
"""
function calc_formal_rv_precission(spectrum::ST, chunklist::ACL; smooth_factor::Real = 2.0)  where {  ST<:AbstractSpectra2D, ACL<:AbstractChunkList }
  σ_rv_combo = 0
  for chid in 1:length(chunklist)
   (depth, σ_rv) = calc_depth_and_expected_rv_precission(spectrum,chunklist[chid].flux.indices[1],chunklist[chid].flux.indices[2], smooth_factor=smooth_factor)
   σ_rv_combo += 1/σ_rv^2
  end
  σ_rv_combo = sqrt(1/σ_rv_combo)
end

""" `calc_formal_rv_precission(spectra, chunklist_timeseries; smooth_factor)`
Calculate expected RV uncertainty for each spectrum in an array of spectra and corresponding chunklist_timeseries,
Assumies only photon noise with given variances.
"""
function calc_formal_rv_precission(spectra::AS, chunk_list_timeseries::ACLT; smooth_factor::Real = 2.0, verbose::Bool = false) where { ST<:AbstractSpectra, AS<:AbstractArray{ST,1}, ACLT<:AbstractChunkListTimeseries }
  @assert length(spectra) == length(chunk_list_timeseries)
  map(obs->calc_formal_rv_precission(spectra[obs], chunk_list_timeseries[obs], smooth_factor=smooth_factor), 1:length(spectra))
end
