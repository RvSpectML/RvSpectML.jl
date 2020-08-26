using Interpolations

""" Return interpolator for fluxes in spectra. """
function make_interpolator_linear_flux(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectrum}
    LinearInterpolation(spectra.λ, spectra.flux)
end

""" Return interpolator for variances in spectra. """
function make_interpolator_linear_var(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectrum}
    LinearInterpolation(spectra.λ, spectra.var)
end

"""
 interp_chunk_to_grid_linear!( flux_out, var_out, chunk_of_spectrum, wavelengths )
Return spectra interpolated onto a grid of points using linear interpolation.
"""
function interp_chunk_to_grid_linear!( flux_out::AbstractArray{T1,1}, var_out::AbstractArray{T2,1}, chunk::AC, grid::AR ) where {
    T1<:Real, T2<:Real, AC<:AbstractChuckOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    lin_interp_flux = make_interpolator_linear_flux(chunk)
    lin_interp_var = make_interpolator_linear_var(chunk)
    flux_out .= lin_interp_flux(grid)
    var_out .= lin_interp_var(grid)
    return flux_out
end

#using Stheno
#using TemporalGPs
# NOTE:  Using own GP code for now, since include predicting derivatives and can minimize unnecessary dependancies
# TODO:  Revisit if can use TemporalGPs for improved speed
"""
 interp_chunk_to_grid_gp!( flux_out, var_out, chunk_of_spectrum, wavelengths )
Return spectra interpolated onto a grid of points using Gaussian Process interpolation.
"""
function interp_chunk_to_grid_gp!( flux_out::AbstractArray{T1,1}, var_out::AbstractArray{T2,1}, chunk::AC, grid::AR ) where { T1<:Real, T2<:Real, AC<:AbstractChuckOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    @warn "Placeholder... doesn't use GP yet."
    lin_interp_flux = make_interpolator_linear_flux(chunk)
    lin_interp_var = make_interpolator_linear_var(chunk)
    flux_out .= lin_interp_flux(grid)
    var_out .= lin_interp_var(grid)
    return flux_out
end


function pack_chunk_list_timeseries_to_matrix_linear(timeseries::ACLT, chunk_grids::Union{AR,AAV}; alg::Symbol = :Linear) where {
        ACLT<:AbstractChunkListTimeseries, RT<:AbstractRange, AR<:AbstractArray{RT,1}, T<:Real, AV<:AbstractVector{T}, AAV<:AbstractArray{AV,1} }
    @assert alg == :Linear || alg == :GP  # TODO: Eventually move to traits-based system?
    num_obs = length(timeseries)
    num_λ = sum(length.(chunk_grids))
    flux_matrix = Array{Float64,2}(undef,num_λ,num_obs)
    var_matrix = Array{Float64,2}(undef,num_λ,num_obs)
    λ_vec = Array{Float64,1}(undef,num_λ)
    chunk_map = Array{UnitRange{Int64}}(undef,length(chunk_grids))


   for t in 1:num_obs
       idx_start = 0
       for c in 1:length(chunk_grids)
           idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
           if alg == :Linear
               interp_chunk_to_grid_linear!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c])
           elseif alg == :GP
               interp_chunk_to_grid_gp!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c])
           end
           if t == 1
               λ_vec[idx] .= chunk_grids[c]
               chunk_map[c] = idx
           end
           idx_start += length(chunk_grids[c])
       end
   end

   return SpectralTimeSeriesCommonWavelengths(λ_vec,flux_matrix,var_matrix,chunk_map, Generic1D() )
end


#=
using Stheno

function make_interpolator_gp(spectra::Union{AS,AC}; length_scale::Real = 0.1, σ_scale::Real = 1.0) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectrum}
    xobs = spectra.λ
    yobs = copy(spectra.flux)
    obs_var = spectra.var
    med_y = mean(yobs)

    # Choose the length-scale and variance of the process.
    σ² = (σ_scale * med_y)^2
    # Construct a kernel with this variance and length scale.
    k = σ² * stretch(Matern52(), 1 / length_scale)

    # Specify a zero-mean GP with this kernel. Don't worry about the GPC object.
    f = GP(med_y, k, GPC())
    fx = f(xobs, obs_var)
    f_posterior = f | Obs(fx, yobs)
end

function interp_to_grid(spectra::Union{AS,AC}, grid::AR) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectrum, AR<:AbstractRange}
   #grid = chunk_grids[c]
   #make_interpolator_linear(spectra).(grid)
   f_posterior = make_interpolator_gp(spectra,length_scale=6e-5*mean(spectra.λ))
   (mean=mean(f_posterior(grid)),  std=std.(marginals(gp_interp(grid))))
end


function pack_chunks_into_matrix(timeseries::ACLT, chunk_grids::AR) where { ACLT<:AbstractChunkListTimeseries, RT<:AbstractRange, AR<:AbstractArray{RT,1} }
   num_obs = length(timeseries)
   num_λ = sum(length.(chunk_grids))
   flux_matrix = Array{Float64,2}(undef,num_λ,num_obs)
   var_matrix = Array{Float64,2}(undef,num_λ,num_obs)
   λ_vec = Array{Float64,1}(undef,num_λ)
   chunk_map = Array{UnitRange{Int64}}(undef,length(chunk_grids))

   for t in 1:num_obs
        idx_start = 0
        for c in 1:length(chunk_grids)
            if length(chunk_grids[c]) <= 512
                gp_interp = make_interpolator_gp(timeseries.chunk_list[t].data[c],length_scale=1e-4*mean(chunk_grids[c]))
                idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
                flux_matrix[idx,t] .= mean(gp_interp(chunk_grids[c]))
                var_matrix[idx,t] .= var.(marginals(gp_interp(chunk_grids[c])))
                if t == 1
                    λ_vec[idx] .= chunk_grids[c]
                    chunk_map[c] = idx
                    #= if 8500<idx_start <8700
                        println(c,": ",idx)
                    end
                    =#
                end
                idx_start += length(chunk_grids[c])
            else
                lin_interp_flux = make_interpolator_linear_flux(timeseries.chunk_list[t].data[c])
                lin_interp_var = make_interpolator_linear_var(timeseries.chunk_list[t].data[c])
                idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
                flux_matrix[idx,t] .= lin_interp_flux.(chunk_grids[c])
                var_matrix[idx,t] .= lin_interp_var.(chunk_grids[c])
                if t == 1
                    λ_vec[idx] .= chunk_grids[c]
                    chunk_map[c] = idx
                    #= if 8500<idx_start <8700
                        println(c,": ",idx)
                    end
                    =#
                end
                idx_start += length(chunk_grids[c])
            end
        end
    end

    #return flux_matrix
    return (flux=flux_matrix, var=var_matrix, λ=λ_vec, chunk_map=chunk_map)
end

=#

""" Return mean flux (averaging over observations at different times, variance weighted) based on a common set of wavelengths.
   Inputs: flux & var (2d: pixel, time)
"""
function calc_mean_spectrum(flux::AbstractArray{T1,2}, var::AbstractArray{T2,2} ) where { T1<:Real, T2<:Real }
    flux_mean = vec(sum(flux./var,dims=2)./sum(1.0./var,dims=2))
end

""" Estimate numerical derivative of fluxes given wavelengths. """
function calc_deriv(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
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
function calc_mean_deriv(flux::AbstractArray{T1,2}, var::AbstractArray{T1,2}, λ::AbstractArray{T3,1},
        chunk_map::AbstractArray{URT,1}) where
    { T1<:Real, T2<:Real, T3<:Real, URT<:AbstractUnitRange} #, RT<:AbstractRange }
    flux_mean = calc_mean_spectrum(flux,var)
    deriv = Array{T1,1}(undef,length(flux_mean))
    map(c->deriv[c] .= calc_deriv(flux_mean[c],λ[c]),chunk_map )
    return deriv
end

function calc_rvs_from_taylor_expansion(spectra::STS; mean::MT = calc_mean_spectrum(spectra),
                deriv::DT = calc_mean_deriv(spectra), idx::RT = 1:length(mean),
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
                deriv::DT = calc_mean_deriv(spectra),
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
