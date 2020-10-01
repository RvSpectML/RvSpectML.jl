"""
Code for estimating the radial velocity based on projecting the flux onto the derivative of a template spectrum.

Author: Eric Ford
"""

""" Return mean flux (averaging over observations at different times, variance weighted) based on a common set of wavelengths.
   Inputs: flux & var (2d: pixel, time)
"""
function calc_mean_spectrum(flux::AbstractArray{T1,2}, var::AbstractArray{T2,2} ) where { T1<:Real, T2<:Real }
    flux_mean = vec(sum(flux./var,dims=2)./sum(1.0./var,dims=2))
end

""" Return mean numerical derivative (dflux/dlnλ) based on a common set of wavelengths.
Inputs:
- flux (2d)
- var (2d)
- λ (1d)
- chunk_map: Array of ranges specifying how each chunk maps into indexes of output derivative
Output:
- dflux/dlnλ: (1d)
 """
function calc_mean_dfluxdlnlambda(flux::AbstractArray{T1,2}, var::AbstractArray{T1,2}, λ::AbstractArray{T3,1},
        chunk_map::AbstractArray{URT,1}) where
    { T1<:Real, T2<:Real, T3<:Real, URT<:AbstractUnitRange} #, RT<:AbstractRange }
    flux_mean = calc_mean_spectrum(flux,var)
    deriv = Array{T1,1}(undef,length(flux_mean))
    #map(c->deriv[c] .= calc_dfluxdlnlambda(flux_mean[c],λ[c]),chunk_map )
    for c in chunk_map
         calc_dfluxdlnlambda!(deriv[c],flux_mean[c],λ[c])
    end
    return deriv
end

""" Return mean numerical derivative (d²flux/dlnλ²) based on a common set of wavelengths.
Inputs:
- flux (2d)
- var (2d)
- λ (1d)
- chunk_map: Array of ranges specifying how each chunk maps into indexes of output derivative
Output:
- d²flux/dlnλ²: (1d)
 """
 function calc_mean_d2fluxdlnlambda2(flux::AbstractArray{T1,2}, var::AbstractArray{T1,2}, λ::AbstractArray{T3,1},
        chunk_map::AbstractArray{URT,1}) where
    { T1<:Real, T2<:Real, T3<:Real, URT<:AbstractUnitRange} #, RT<:AbstractRange }
    flux_mean = calc_mean_spectrum(flux,var)
    deriv2 = Array{T1,1}(undef,length(flux_mean))
    #map(c->deriv2[c] .= calc_d2fluxdlnlambda2(flux_mean[c],λ[c]),chunk_map )
    for c in chunk_map
         calc_d2fluxdlnlambda2!(deriv2[c],flux_mean[c],λ[c])
    end
    return deriv2
end

""" `calc_rvs_from_taylor_expansion( spectra )`
Inputs:
- spectra: SpectralTimeSeriesCommonWavelengths
Optional Arguments:
- mean: mean spectrum
- deriv: dmeanflux/dlnλ
- idx: range of pixel indcies to use for calculation
- equal_weight:  For now, spectra are equal weighted. (true)
Output named pair:
- rv: Vector of estimated radial velocities
- σ_rv: Vector of uncertainties in rv estimates
"""
function calc_rvs_from_taylor_expansion(spectra::STS; mean::MT = calc_mean_spectrum(spectra),
                deriv::DT = calc_mean_dfluxdlnlambda(spectra), idx::RT = 1:length(mean),
                equal_weight::Bool = true ) where {
                    STS<:AbstractSpectralTimeSeriesCommonWavelengths, T1<:Real, MT<:AbstractVector{T1},
                    T2<:Real, DT<:AbstractVector{T2}, RT<:AbstractUnitRange }
   @assert length(mean) == length(deriv)
   @assert size(spectra.flux,1) == length(mean)
   @assert minimum(idx) >=1
   @assert maximum(idx) <= length(mean)

   if equal_weight # Pixels equally-weighted
      norm = sum(abs2.(deriv[idx]))
      rv = sum((spectra.flux[idx,:].-mean[idx]).*deriv[idx],dims=1).*(speed_of_light_mps/norm)
      # TODO: WARN: Uncertinaties ignores correlations between pixels, particularly problematic when oversample pixels
      σ_rv = sqrt.(sum(spectra.var[idx,:].*abs2.(deriv[idx]),dims=1)) .*(speed_of_light_mps/norm)
  else # Pixels inverse variance weighted
      # TODO:  CHECK/FIX? math for inverse variance weighting.
      # Or maybe the math is right, but it's just a bad idea to have different weighting within one line/chunk
      @info "Warning: I think this is either wrong or a bad idea.  Need to check."
      norm = sum(abs2.(deriv[idx])./spectra.var[idx,:])
      rv = sum((spectra.flux[idx,:].-mean[idx]).*deriv[idx]./spectra.var[idx,:],dims=1).*(speed_of_light_mps/norm)
      σ_rv = sqrt.(sum(abs2.(deriv)./spectra.var[idx,:],dims=1)).*(speed_of_light_mps/norm)
   end
   return (rv=-vec(rv), σ_rv=vec(σ_rv))
end

""" `calc_rvs_from_taylor_expansion_alt( spectra )`
Experimental version of [calc_rvs_from_taylor_expansion](@ref).
"""
function calc_rvs_from_taylor_expansion_alt(spectra::STS; mean::MT = calc_mean_spectrum(spectra),
                deriv::DT = calc_mean_dfluxdlnlambda(spectra),
                deriv2::DT = calc_d2fluxdlnlambda2(spectra), idx::RT = 1:length(mean),
                equal_weight::Bool = true ) where {
                    STS<:AbstractSpectralTimeSeriesCommonWavelengths, T1<:Real, MT<:AbstractVector{T1},
                    T2<:Real, DT<:AbstractVector{T2}, RT<:AbstractUnitRange }
   @assert length(mean) == length(deriv)
   @assert size(spectra.flux,1) == length(mean)
   @assert minimum(idx) >=1
   @assert maximum(idx) <= length(mean)

   if equal_weight # Pixels equally-weighted
      #norm = sum(abs2.(deriv[idx]))
      norm = sum(deriv[idx].^2 .- (spectra.flux[idx,:].-mean[idx]).*deriv2[idx],dims=1)
      rv = sum(((spectra.flux[idx,:].-mean[idx]).*deriv[idx]),dims=1)./norm
      rv *= speed_of_light_mps
      # TODO: WARN: Uncertinaties ignores correlations between pixels, particularly problematic when oversample pixels
      σ_rv = sqrt.(sum(spectra.var[idx,:].*abs2.(deriv[idx]),dims=1)./norm)
      σ_rv *= speed_of_light_mps
  else # Pixels inverse variance weighted
      # TODO:  CHECK/FIX? math for inverse variance weighting.
      # Or maybe the math is right, but it's just a bad idea to have different weighting within one line/chunk
      @info "Warning: I think this is either wrong or a bad idea.  Need to check."
      norm = sum(abs2.(deriv[idx]./spectra.var[idx,:]))  # TODO comparable to norm_old, if new norm works, then check updating this also helps
      rv = sum((spectra.flux[idx,:].-mean[idx]).*deriv[idx]./spectra.var[idx,:],dims=1).*(speed_of_light_mps/norm)
      σ_rv = sqrt.(sum(abs2.(deriv)./spectra.var[idx,:],dims=1)).*(speed_of_light_mps/norm)
   end
   return (rv=-vec(rv), σ_rv=vec(σ_rv))
end

""" `calc_chunk_rvs_from_taylor_expansion( spectra )`
Inputs:
- spectra: SpectralTimeSeriesCommonWavelengths
Optional Arguments:
- mean: mean spectrum
- deriv: dmeanflux/dlnλ
- equal_weight:  For now, spectra are equal weighted. (true)
Output named pair:
- rv: Vector of vector of estimated radial velocities for each chunk
- σ_rv: Vector of vector of uncertainties in rv estimates for each chunk
"""
function calc_chunk_rvs_from_taylor_expansion(spectra::STS; mean::MT = calc_mean_spectrum(spectra),
                deriv::DT = calc_mean_dfluxdlnlambda(spectra),
                equal_weight::Bool = false ) where {
                    STS<:AbstractSpectralTimeSeriesCommonWavelengths, T1<:Real, MT<:AbstractVector{T1},
                    T2<:Real, DT<:AbstractVector{T2}, RT<:AbstractUnitRange }
   @assert length(mean) == length(deriv)
   @assert size(spectra.flux,1) == length(mean)

   map(idx->calc_rvs_from_taylor_expansion(spectra,mean=mean,deriv=deriv,idx=idx),spectra.chunk_map)
end
