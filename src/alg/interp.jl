""" Return interpolator for fluxes in spectra. """
function make_interpolator_linear_flux(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChunkOfSpectra}
    LinearInterpolation(spectra.λ, spectra.flux)
end

""" Return interpolator for variances in spectra. """
function make_interpolator_linear_var(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChunkOfSpectra}
    LinearInterpolation(spectra.λ, spectra.var)
end

#=
using Stheno

function make_interpolator_gp(spectra::Union{AS,AC}; length_scale::Real = 0.1, σ_scale::Real = 1.0) where { AS<:AbstractSpectra, AC<:AbstractChunkOfSpectra}
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

function interp_to_grid(spectra::Union{AS,AC}, grid::AR) where { AS<:AbstractSpectra, AC<:AbstractChunkOfSpectra, AR<:AbstractRange}
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
    flux_mean = sum(flux./vm,dims=2)./sum(1.0./vm,dims=2)
    deriv = Array{T1,1}(undef,length(flux_mean))
    map(c->deriv[c] .= calc_deriv(flux_mean[c],λv[c]),chunk_map )
    return deriv
end
