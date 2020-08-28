
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
    dfdlogλ[1] = 0.5*(flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    dfdlogλ[2:end-1] .= 0.5*(flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[end-2])
    dfdlogλ[end] = 0.5*(flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

""" Return mean numerical derivative of fluxes based on a common set of wavelengths.
    Inputs: flux & var (2d) and λ (1d)
 """
function calc_mean_dfluxdlnlambda(flux::AbstractArray{T1,2}, var::AbstractArray{T1,2}, λ::AbstractArray{T3,1},
        chunk_map::AbstractArray{URT,1}) where
    { T1<:Real, T2<:Real, T3<:Real, URT<:AbstractUnitRange} #, RT<:AbstractRange }
    flux_mean = calc_mean_spectrum(flux,var)
    deriv = Array{T1,1}(undef,length(flux_mean))
    map(c->deriv[c] .= calc_dfluxdlnlambda(flux_mean[c],λ[c]),chunk_map )
    return deriv
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
      @info "Warning: I think this is either wrong or a bad idea."
      norm = sum(abs2.(deriv[idx])./spectra.var[idx,:])
      rv = sum((spectra.flux[idx,:].-mean[idx]).*deriv[idx]./spectra.var[idx,:],dims=1).*(speed_of_light_mps/norm)
      σ_rv = sqrt.(sum(abs2.(deriv)./spectra.var[idx,:],dims=1)).*(speed_of_light_mps/norm)
   end
   return (rv=vec(rv), σ_rv=vec(σ_rv))
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


function compute_spectra_perp_doppler_shift(spectra::AA, deriv::V1, rvs::V2) where {
            T1<:Real, AA<:AbstractArray{T1,2}, T2<:Real, V1<:AbstractVector{T2}, T3<:Real, V2<:AbstractVector{T3} }
   @assert size(spectra,1) == length(deriv)
   @assert size(spectra,2) == length(rvs)
   fm_perp = spectra .- rvs' .* deriv./RvSpectML.speed_of_light_mps
end