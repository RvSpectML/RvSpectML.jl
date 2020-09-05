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

""" Estimate numerical derivative of fluxes given wavelengths. """
function calc_dfluxdlnlambda(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
    @assert size(flux) == size(λ)
    @assert length(flux) >= 3
    dfdlogλ = Array{T1,1}(undef,length(flux))
    calc_dfluxdlnlambda!(dfdlogλ,flux,λ)
    #dfdlogλ[1] = 0.5 * (flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    #dfdlogλ[2:end-1] .= 0.5 .* (flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[1:end-2])
    #dfdlogλ[end] = 0.5 * (flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

""" Estimate numerical derivative of fluxes given wavelengths. """
function calc_dfluxdlnlambda!(dfdlogλ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, λ::AbstractArray{T3,1}) where { T1<:Real, T2<:Real, T3<:Real }
    @assert size(flux) == size(λ)
    @assert size(dfdlogλ) == size(flux)
    @assert length(flux) >= 3
    #dfdlogλ = Array{T1,1}(undef,length(flux))
    dfdlogλ[1] = 0.5 * (flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    dfdlogλ[2:end-1] .= 0.5 .* (flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[1:end-2])
    dfdlogλ[end] = 0.5 * (flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

""" Estimate numerical second derivative of fluxes given wavelengths. """
function calc_d2fluxdlnlambda2(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
    @assert size(flux) == size(λ)
    @assert length(flux) >= 3
    logλ = log.(λ)
    d2fdlogλ2 = Array{T1,1}(undef,length(flux))
    #d2fdlogλ2[2:end-1] .= 0.25*(flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1])./(λ[3:end].+λ[1:end-2].-2.0*λ[2:end-1]).*(λ[3:end].+λ[end-2]).^2
    #d2fdlogλ2[2:end-1] .= 0.5 * (flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1]).* ((λ[3:end].+λ[1:end-2])./(logλ[3:end].-logλ[1:end-2])).^2
    #d2fdlogλ2[end] = d2fdlogλ2[end-1]
    #d2fdlogλ2[1] = d2fdlogλ2[2]
    calc_d2fluxdlnlambda2!(d2fdlogλ2,flux,λ)
    return d2fdlogλ2
end

""" Estimate numerical second derivative of fluxes given wavelengths. """
function calc_d2fluxdlnlambda2!(d2fdlogλ2::AbstractArray{T1,1}, flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real, T3<:Real }
    @assert size(flux) == size(λ)
    @assert size(d2fdlogλ2) == size(flux)
    @assert length(flux) >= 3
    #logλ = log.(λ)
    #d2fdlogλ2 = Array{T1,1}(undef,length(flux))
    #d2fdlogλ2[2:end-1] .= 0.25*(flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1])./(λ[3:end].+λ[1:end-2].-2.0*λ[2:end-1]).*(λ[3:end].+λ[end-2]).^2
    d2fdlogλ2[2:end-1] .= 0.5 * (flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1]).* ((λ[3:end].+λ[1:end-2]).^2 ./ (2.0.*log.(λ[3:end]./λ[1:end-2])))
    d2fdlogλ2[end] = d2fdlogλ2[end-1]
    d2fdlogλ2[1] = d2fdlogλ2[2]
    return d2fdlogλ2
end

""" Return mean numerical derivative of fluxes based on a common set of wavelengths.
    Inputs: flux & var (2d) and λ (1d)
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


function calc_chunk_rvs_from_taylor_expansion(spectra::STS; mean::MT = calc_mean_spectrum(spectra),
                deriv::DT = calc_mean_dfluxdlnlambda(spectra),
                equal_weight::Bool = false ) where {
                    STS<:AbstractSpectralTimeSeriesCommonWavelengths, T1<:Real, MT<:AbstractVector{T1},
                    T2<:Real, DT<:AbstractVector{T2}, RT<:AbstractUnitRange }
   @assert length(mean) == length(deriv)
   @assert size(spectra.flux,1) == length(mean)

   map(idx->calc_rvs_from_taylor_expansion(spectra,mean=mean,deriv=deriv,idx=idx),spectra.chunk_map)
end
