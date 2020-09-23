"""
Delegates loading of code for different interpolation aglorithms (in their own modules).
Also wraps these modules, so they're called easily by `pack_chunk_list_timeseries_to_matrix`.

Author: Eric Ford
Created: August 2020
"""

# First include modules for different interpolation algoritihms
# File called module will include whatever packages it depends on, and also
# provide at least one function with a common interface (currently named) interp_chunk_to_grid_ALGNAME!.

include("linear.jl")
import .LinearInterpolation
export interp_chunk_to_grid_linear!

include("sinc.jl")
using .SincInterpolation
export interp_chunk_to_grid_sinc!

#=
include("gp/gp_brute_force.jl")             # Basic, brute-force GPs with no dependancies
import .GPInterpolation
export interp_chunk_to_grid_gp!
=#

include("gp/temporalgps.jl")    # Dependenacies on Stheno, Temporal GPs and Static Arrays
import .TemporalGPInterpolation
export interp_chunk_to_grid_temporalgp!

"""
   interp_chunk_to_grid_linear!( flux_out, var_out, chunk_of_spectrum, wavelengths )
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- flux_out::AbstractArray
- var_out::AbstractArray
- chunk::AbstractChunkOfSpectrum
- grid::AbstractRange or AbstractArray
# Alters
- flux_out
- var_out::AbstractArray
# Returns
- flux_out
"""
function interp_chunk_to_shifted_grid_linear!( flux_out::AbstractArray{T1,1}, var_out::AbstractArray{T2,1}, chunk::AC, grid::AR, boost_factor::Real ) where {
    T1<:Real, T2<:Real, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    lin_interp_flux = #LinearInterpolation.extrapolate(
        LinearInterpolation.make_interpolator_linear_flux(chunk) #, Flat() )
    lin_interp_var = #LinearInterpolation.extrapolate(
        LinearInterpolation.make_interpolator_linear_var(chunk)  #, Flat() )
    flux_out .= lin_interp_flux(grid./boost_factor)
    var_out .= lin_interp_var(grid./boost_factor)
    return flux_out
end

function interp_chunk_to_shifted_grid_linear( chunk::AC, grid::AR, boost_factor::Real ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_shifted_grid_linear!(flux_out, var_out, chunk, grid, boost_factor)
    return (flux=flux_out, var=var_out)
end

function interp_chunk_to_grid_linear!( flux_out::AbstractArray{T1,1}, var_out::AbstractArray{T2,1}, chunk::AC, grid::AR ) where {
    T1<:Real, T2<:Real, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    lin_interp_flux = LinearInterpolation.make_interpolator_linear_flux(chunk)
    lin_interp_var = LinearInterpolation.make_interpolator_linear_var(chunk)
    flux_out .= lin_interp_flux(grid)
    var_out .= lin_interp_var(grid)
    return flux_out
end

function interp_chunk_to_grid_linear( chunk::AC, grid::AR ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_linear!(flux_out, var_out, chunk, grid)
    return (flux=flux_out, var=var_out)
end

"""
   interp_chunk_to_grid_sinc!( flux_out, var_out, chunk_of_spectrum, wavelengths )
Return spectra interpolated onto a grid of points using sinc interpolation.
# Arguments:
- flux_out::AbstractArray
- var_out::AbstractArray
- chunk::AbstractChunkOfSpectrum
- grid::AbstractRange or AbstractArray
# Alters
- flux_out
- var_out::AbstractArray
# Returns
- flux_out
"""
function interp_chunk_to_shifted_grid_sinc!( flux_out::AbstractArray{T1,1}, var_out::AbstractArray{T2,1}, chunk::AC, grid::AR, boost_factor::Real; Filter::AbstractVector = zeros(0) ) where {
    T1<:Real, T2<:Real, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    if length(Filter) >= 1
        flux_out .= SincInterpolation.spectra_interpolate(grid./ boost_factor,chunk.λ,chunk.flux, Filter=Filter)
        var_out .= SincInterpolation.spectra_interpolate(grid./ boost_factor,chunk.λ ,chunk.var, Filter=Filter)
    else
        flux_out .= SincInterpolation.spectra_interpolate(grid./ boost_factor,chunk.λ ,chunk.flux)
        var_out .= SincInterpolation.spectra_interpolate(grid ./ boost_factor,chunk.λ,chunk.var)
    end
    return flux_out
end

function interp_chunk_to_shifted_grid_sinc( chunk::AC, grid::AR, boost_factor::Real; Filter::AbstractVector = zeros(0) ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    if length(Filter) >= 1
        interp_chunk_to_shifted_grid_sinc!(flux_out, var_out, chunk, grid, boost_factor, Filter=Filter)
    else
        interp_chunk_to_shifted_grid_sinc!(flux_out, var_out, chunk, grid, boost_factor)
    end
    return (flux=flux_out, var=var_out)
end

function interp_chunk_to_grid_sinc!( flux_out::AbstractArray{T1,1}, var_out::AbstractArray{T2,1}, chunk::AC, grid::AR; Filter::AbstractVector = zeros(0) ) where {
    T1<:Real, T2<:Real, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    if length(Filter) >= 1
        flux_out .= SincInterpolation.spectra_interpolate(grid,chunk.λ,chunk.flux,Filter=Filter)
        var_out .= SincInterpolation.spectra_interpolate(grid,chunk.λ,chunk.var,Filter=Filter)
    else
        flux_out .= SincInterpolation.spectra_interpolate(grid,chunk.λ,chunk.flux)
        var_out .= SincInterpolation.spectra_interpolate(grid,chunk.λ,chunk.var)
    end
    return flux_out
end

function interp_chunk_to_grid_sinc( chunk::AC, grid::AR ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_sinc!(flux_out, var_out, chunk, grid)
    return (flux=flux_out, var=var_out)
end

"""
   interp_chunk_to_grid_gp!( flux_out, var_out, chunk_of_spectrum, wavelengths )
Return spectra interpolated onto a grid of points using Gaussian Process interpolation.
# Arguments:
- flux_out::AbstractArray
- var_out::AbstractArray
- chunk::AbstractChunkOfSpectrum
- grid::AbstractRange or AbstractArray
# Alters
- flux_out
- var_out::AbstractArray
# Returns
- flux_out

NOTE:  Using own GP code for now, since include predicting derivatives and can minimize unnecessary dependancies.  We may need to
revisit this if we want improved speed
"""

function interp_chunk_to_shifted_grid_gp!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR, boost_factor::Real ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    rho = 2*5000/speed_of_light_mps * mean(chunk.λ)
    sigmasq_kernel = 0.25 # Float64(mean(chunk.var))
    println(" rho = ", rho, "  σ²_kernel = ", sigmasq_kernel)
    flux_out .= GPInterpolation.predict_mean(chunk.λ./boost_factor, chunk.flux, grid, sigmasq_obs = chunk.var, kernel = GPs.matern52_sparse_kernel, rho=rho, sigmasq_cor=sigmasq_kernel ) # 	sigmasq_cor=1.0, rho=1
    # TODO: Update var_out to actually use the right GP or at least do something more sensible
    #var_out .= GPInterpolation.predict_mean(chunk.λ, chunk.var, grid, sigmasq_obs = chunk.var, kernel = GPs.matern52_sparse_kernel, rho=5000/speed_of_light_mps) # 	sigmasq_cor=1.0, rho=1
    return flux_out
end

function interp_chunk_to_shifted_grid_gp( chunk::AC, grid::AR, boost_factor::Real ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp!(flux_out, var_out, chunk, grid, boost_factor=boost_factor)
    return (flux=flux_out, var=var_out)
end

function interp_chunk_to_grid_gp!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    rho = 2*5000/speed_of_light_mps * mean(chunk.λ)
    sigmasq_kernel = 0.25 # Float64(mean(chunk.var))
    println(" rho = ", rho, "  σ²_kernel = ", sigmasq_kernel)
    flux_out .= GPInterpolation.predict_mean(chunk.λ, chunk.flux, grid, sigmasq_obs = chunk.var, kernel = GPs.matern52_sparse_kernel, rho=rho, sigmasq_cor=sigmasq_kernel ) # 	sigmasq_cor=1.0, rho=1
    # TODO: Update var_out to actually use the right GP or at least do something more sensible
    #var_out .= GPInterpolation.predict_mean(chunk.λ, chunk.var, grid, sigmasq_obs = chunk.var, kernel = GPs.matern52_sparse_kernel, rho=5000/speed_of_light_mps) # 	sigmasq_cor=1.0, rho=1
    return flux_out
end

function interp_chunk_to_grid_gp( chunk::AC, grid::AR ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp!(flux_out, var_out, chunk, grid, )
    return (flux=flux_out, var=var_out)
end




function interp_chunk_to_shifted_grid_gp_temporal!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR, boost_factor::Real; use_logx::Bool = true,  use_logy::Bool = false, smooth_factor::Real=1 ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    y = eltype(chunk.flux) == Float64  ? chunk.flux : convert.(Float64,chunk.flux)
    var = eltype(chunk.var) == Float64  ? chunk.var : convert.(Float64,chunk.var)

    if !use_logy
        mean_y = mean(y)
        y_to_use = y./mean_y .- 1.0
    else
        y_to_use = y
    end
    #y_to_use = y
    flux_out .= TemporalGPInterpolation.predict_mean(chunk.λ./boost_factor, y_to_use, grid, sigmasq_obs = var, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    flux_out .= mean_y.*(flux_out.+1.0)
    # TODO: Update var_out to actually use the right GP or at least do something more sensible
    var_out .= TemporalGPInterpolation.predict_mean(chunk.λ./boost_factor, var, grid, sigmasq_obs = var, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    return flux_out
end

function interp_chunk_to_shifted_grid_gp_temporal( chunk::AC, grid::AR, boost_factor::Real; use_logx::Bool = true,  use_logy::Bool = false,smooth_factor::Real=1 ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp!(flux_out, var_out, chunk, grid, boost_factor=boost_factor, smooth_factor=smooth_factor)
    return (flux=flux_out, var=var_out)
end


function interp_chunk_to_grid_gp_temporal!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR ; use_logx::Bool = true,  use_logy::Bool = false, smooth_factor::Real=1 ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    y = eltype(chunk.flux) == Float64  ? chunk.flux : convert.(Float64,chunk.flux)
    var = eltype(chunk.var) == Float64  ? chunk.var : convert.(Float64,chunk.var)
    if !use_logy
        mean_y = mean(y)
        y_to_use = y./mean_y .- 1.0
    else
        y_to_use = y
    end
    y_to_use = y
    flux_out .= TemporalGPInterpolation.predict_mean(chunk.λ, y_to_use, grid, sigmasq_obs = var, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    flux_out .= mean_y.*(flux_out.+1.0)
    # TODO: Update var_out to actually use the right GP or at least do something more sensible
    var_out .= flux_out
    var_out .= TemporalGPInterpolation.predict_mean(chunk.λ, var, grid, sigmasq_obs = var, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    return flux_out
end

function interp_chunk_to_grid_gp_temporal( chunk::AC, grid::AR ; use_logx::Bool = true,  use_logy::Bool = false, smooth_factor::Real=1 ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp!(flux_out, var_out, chunk, grid, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor)
    return (flux=flux_out, var=var_out)
end



# Wrapper code to deal with weird data structures
function pack_chunk_list_timeseries_to_matrix(timeseries::ACLT, chunk_grids::Union{AR,AAV}; alg::Symbol = :TemporalGP,
    oversample_factor::Real = 1, smooth_factor::Real=1, verbose::Bool = false ) where {
        ACLT<:AbstractChunkListTimeseries, RT<:AbstractRange, AR<:AbstractArray{RT,1}, T<:Real, AV<:AbstractVector{T}, AAV<:AbstractArray{AV,1} }
    @assert alg == :Linear || alg == :GP  || alg == :Sinc || alg == :TemporalGP # TODO: Eventually move to traits-based system?
    num_obs = length(timeseries)
    num_λ = sum(length.(chunk_grids))
    flux_matrix = Array{Float64,2}(undef,num_λ,num_obs)
    var_matrix = Array{Float64,2}(undef,num_λ,num_obs)
    λ_vec = Array{Float64,1}(undef,num_λ)
    chunk_map = Array{UnitRange{Int64}}(undef,length(chunk_grids))
    mean_flux = zeros(num_λ)
    mean_var = zeros(num_λ)
    dfluxdlnλ = zeros(num_λ)
    d2fluxdlnλ2 = zeros(num_λ)
    if (alg == :GP) && (maximum(length.(chunk_grids))>1024)
        @error "To use a GP with more than 1024 points in a chunk, use a more efficient GP package, e.g., :TemporalGP."
    end
    if alg == :Sinc    # Setup workspace for Sync.  TODO: Put into functor
        filter_size=23
        kaiserB=13
        sinc_filter = SincInterpolation.create_filter_curve(filter_size*21; filter_size=filter_size, kaiserB=kaiserB)
    end

    for t in 1:num_obs
       idx_start = 0
       for c in 1:length(chunk_grids)
           idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
           if verbose && !issorted(timeseries.chunk_list[t].data[c].λ)
               flush(stdout);  println("t= ",t, " c= ",c," idx= ", idx, " size(flux)= ",size(flux_matrix), " issorted(λ)=", issorted(timeseries.chunk_list[t].data[c].λ))
           end
           if t == 1
               λ_vec[idx] .= chunk_grids[c]
               chunk_map[c] = idx
           end
           mean_flux_in_chunk = mean(timeseries.chunk_list[t].data[c].flux)
           if mean_flux_in_chunk > 0
               if alg == :Linear
                   interp_chunk_to_grid_linear!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c])
               elseif alg == :Sinc
                   interp_chunk_to_grid_sinc!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c], Filter=sinc_filter)
               elseif alg == :GP
                   interp_chunk_to_grid_gp!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c])
               elseif alg == :TemporalGP
                   interp_chunk_to_grid_gp_temporal!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c] , use_logx=true, use_logy=false, smooth_factor=smooth_factor)
               end
               mean_flux[idx] .+= flux_matrix[idx,t]
               mean_var[idx]  .+= var_matrix[idx,t]
               dfluxdlnλ[idx] .+= calc_dfluxdlnlambda(view(flux_matrix,idx,t),view(λ_vec,idx))
               d2fluxdlnλ2[idx] .+= calc_d2fluxdlnlambda2(view(flux_matrix,idx,t),view(λ_vec,idx))

               flux_matrix[idx,t] ./= mean_flux_in_chunk
               var_matrix[idx,t]  ./= mean_flux_in_chunk^2
           else  # no flux in chunk?!?
               @warn "No flux in chunk " * string(c) * " at time " * string(t) * "."
           end
           idx_start += length(chunk_grids[c])
       end # for c
   end # for t

    # TODO Think if can propagate uncertainties better
    var_matrix .*= sqrt(oversample_factor)

    #flush(stdout)
    #println("# mean flux pre-normalize = ", mean(mean_flux))
    idx_start = 0
    for c in 1:length(chunk_grids)
        idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
        mean_flux_in_chunk = mean(mean_flux[idx])
        #println("Normalizing chunk ", c, " (",idx,") by ", mean_flux_in_chunk)
        #dmeanfluxdlnλ[idx] .= calc_dfluxdlnlambda(vec(sum(view(flux_matrix,idx,:),dims=2)),vec(sum(view(var_matrix,idx,:), dims=2)) ) ./mean_flux_in_chunk
        mean_flux[idx] ./= mean_flux_in_chunk
        mean_var[idx] ./= mean_flux_in_chunk^2
        #dfluxdlnλ[idx] ./= mean_flux_in_chunk
        #d2fluxdlnλ2[idx] ./= mean_flux_in_chunk
        dfluxdlnλ[idx] .= calc_dfluxdlnlambda(mean_flux[idx],view(λ_vec,idx))
        d2fluxdlnλ2[idx] .= calc_d2fluxdlnlambda2(mean_flux[idx],view(λ_vec,idx))

        #dmeanfluxdlnλ[idx] .= calc_dfluxdlnlambda(vec(sum(view(flux_matrix,idx,:),dims=2)),view(λ_vec,idx))  ./mean_flux_in_chunk
        idx_start += length(chunk_grids[c])
    end
    #println("# mean flux post-normalize = ", mean(mean_flux))
    #flush(stdout)
    #mean_flux /= num_obs
    mean_var .*= sqrt(oversample_factor)
    #dmeanfluxdlnλ ./= num_obs
    return ( matrix=SpectralTimeSeriesCommonWavelengths(λ_vec,flux_matrix,var_matrix,chunk_map, Generic1D() ), mean_flux=mean_flux, mean_var=mean_var, deriv=dfluxdlnλ, deriv2=d2fluxdlnλ2 )
end


# Wrapper code to deal with weird data structures
function pack_shifted_chunk_list_timeseries_to_matrix(timeseries::ACLT, chunk_grids::Union{AR,AAV}; alg::Symbol = :Linear,
    oversample_factor::Real = 1, smooth_factor::Real=1, remove_rv_est::Bool = true, verbose::Bool = false ) where {
        ACLT<:AbstractChunkListTimeseries, RT<:AbstractRange, AR<:AbstractArray{RT,1}, T<:Real, AV<:AbstractVector{T}, AAV<:AbstractArray{AV,1} }
    @assert alg == :Linear || alg == :GP  || alg == :Sinc || alg == :TemporalGP # TODO: Eventually move to traits-based system?
    num_obs = length(timeseries)
    num_λ = sum(length.(chunk_grids))
    flux_matrix = Array{Float64,2}(undef,num_λ,num_obs)
    var_matrix = Array{Float64,2}(undef,num_λ,num_obs)
    λ_vec = Array{Float64,1}(undef,num_λ)
    chunk_map = Array{UnitRange{Int64}}(undef,length(chunk_grids))
    mean_flux = zeros(num_λ)
    mean_var = zeros(num_λ)
    dfluxdlnλ = zeros(num_λ)
    d2fluxdlnλ2 = zeros(num_λ)
    if remove_rv_est
        @assert haskey(first(timeseries.metadata),:rv_est)
    end
    if (alg == :GP) && (maximum(length.(chunk_grids))>1024)
        @error "To use a GP with more than 1024 points in a chunk, use a more efficient GP package, e.g., :TemporalGP."
    end
    if alg == :Sinc    # Setup workspace for Sync.  TODO: Put into functor
        filter_size=23
        kaiserB=13
        sinc_filter = SincInterpolation.create_filter_curve(filter_size*21; filter_size=filter_size, kaiserB=kaiserB)
    end
    for t in 1:num_obs
       idx_start = 0
       boost_factor = remove_rv_est ? calc_doppler_factor(timeseries.metadata[t][:rv_est]) : 1.0
       for c in 1:length(chunk_grids)
           idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
           if verbose
               flush(stdout);  println("t= ",t, " c= ",c," idx= ", idx, " size(flux)= ",size(flux_matrix))
           end
           if t == 1
               λ_vec[idx] .= chunk_grids[c]
               chunk_map[c] = idx
           end
           mean_flux_in_chunk = mean(timeseries.chunk_list[t].data[c].flux)
           if mean_flux_in_chunk > 0
               #println(" mean flux in chunk ", c, " at time ", t, " = ", mean_flux_in_chunk)
               if alg == :Linear
                   interp_chunk_to_shifted_grid_linear!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c], boost_factor )
               elseif alg == :Sinc
                   interp_chunk_to_shifted_grid_sinc!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c], boost_factor, Filter=sinc_filter )
               elseif alg == :GP
                   interp_chunk_to_shifted_grid_gp!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c], boost_factor)
               elseif alg == :TemporalGP
                   interp_chunk_to_shifted_grid_gp_temporal!(view(flux_matrix,idx,t), view(var_matrix,idx,t), timeseries.chunk_list[t].data[c], chunk_grids[c], boost_factor, use_logx=true, use_logy=false, smooth_factor=smooth_factor)
               end
               if t == 1
                  λ_vec[idx] .= chunk_grids[c]
                  chunk_map[c] = idx
               end
            else  # no flux in chunk?!?
                @warn "No flux in chunk " * string(c) * " at time " * string(t) * "."
            end

           mean_flux[idx] .+= view(flux_matrix,idx,t)
           mean_var[idx]  .+= view(var_matrix,idx,t)
           dfluxdlnλ[idx] .+= calc_dfluxdlnlambda(view(flux_matrix,idx,t),view(λ_vec,idx))
           d2fluxdlnλ2[idx] .+= calc_d2fluxdlnlambda2(view(flux_matrix,idx,t),view(λ_vec,idx))

           flux_matrix[idx,t] ./= mean_flux_in_chunk
           var_matrix[idx,t]  ./= mean_flux_in_chunk^2

           idx_start += length(chunk_grids[c])
       end
    end

    # TODO Think if can propagate uncertainties better
    var_matrix .*= sqrt(oversample_factor)

    #flush(stdout)
    #println("# mean flux pre-normalize = ", mean(mean_flux))
    idx_start = 0
    for c in 1:length(chunk_grids)
        idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
        mean_flux_in_chunk = mean(mean_flux[idx])
        #println("Normalizing chunk ", c, " (",idx,") by ", mean_flux_in_chunk)
        #dmeanfluxdlnλ[idx] .= calc_dfluxdlnlambda(vec(sum(view(flux_matrix,idx,:),dims=2)),vec(sum(view(var_matrix,idx,:), dims=2)) ) ./mean_flux_in_chunk
        #sum_in_chunk = vec(sum(view(flux_matrix,idx,:),dims=2))
        mean_flux[idx] ./= mean_flux_in_chunk
        mean_var[idx] ./= mean_flux_in_chunk^2
        #dfluxdlnλ[idx] ./= mean_flux_in_chunk
        #d2fluxdlnλ2[idx] ./= mean_flux_in_chunk
        dfluxdlnλ[idx] .= calc_dfluxdlnlambda(mean_flux[idx],view(λ_vec,idx))
        d2fluxdlnλ2[idx] .= calc_d2fluxdlnlambda2(mean_flux[idx],view(λ_vec,idx))
        idx_start += length(chunk_grids[c])
    end

    #println("# mean flux post-normalize = ", mean(mean_flux))
    #flush(stdout)
    #mean_flux /= num_obs
    mean_var .*= sqrt(oversample_factor)
    #dmeanfluxdlnλ ./= num_obs
    return ( matrix=SpectralTimeSeriesCommonWavelengths(λ_vec,flux_matrix,var_matrix,chunk_map, Generic1D() ), mean_flux=mean_flux, mean_var=mean_var, deriv=dfluxdlnλ, deriv2=d2fluxdlnλ2 )
end

"""   repack_flux_vector_to_chunk_matrix(λ, flux, var, chunk_map, λc; alg )
Warning:  This doesn't work yet
"""
function repack_flux_vector_to_chunk_matrix(λ::AA1, flux::AA2, var::AA3, chunk_map::CMT, λc::AA4; alg::Symbol = :TemporalGP,
    oversample_factor::Real = 1, smooth_factor::Real=1, verbose::Bool = false ) where {
        T1<:Real, AA1<:AbstractVector{T1}, T2<:Real, AA2<:AbstractVector{T2}, T3<:Real, AA3<:AbstractVector{T3}, T4<:Real, AA4<:AbstractVector{T4},
        CMT<:AbstractArray{UnitRange{Int64}} } #RT<:AbstractRange, AR<:AbstractArray{RT,1}, T<:Real, AV<:AbstractVector{T}, AAV<:AbstractArray{AV,1} }
    @assert alg == :Linear || alg == :GP  || alg == :Sinc || alg == :TemporalGP # TODO: Eventually move to traits-based system?
    num_chunks = length(chunk_map)
    num_λ_per_chunk = minimum(length.(chunk_map))
    flux_matrix = zeros(num_λ_per_chunk,num_chunks)
    var_matrix = zeros(num_λ_per_chunk,num_chunks)
    if (alg == :GP) && (num_λ_per_chunk>1024)
        @error "To use a GP with more than 1024 points in a chunk, use a more efficient GP package, e.g., :TemporalGP."
    end
    if alg == :Sinc    # Setup workspace for Sync.  TODO: Put into functor
        filter_size=23
        kaiserB=13
        sinc_filter = SincInterpolation.create_filter_curve(filter_size*21; filter_size=filter_size, kaiserB=kaiserB)
    end

   idx_start = 0
   for c in 1:length(chunk_map)
       idx_in = (idx_start+1):(idx_start+length(chunk_map[c]))
       idx_out_begin = floor(Int,(idx_in[1]+idx_in[end])//2)-floor(Int,num_λ_per_chunk//2)
       idx_out_end = idx_out_begin+num_λ_per_chunk-1
       idx_out = idx_out_begin:idx_out_end
       Δλ_in_chunk = λ[idx_out].-λc[c]
       mean_flux_in_chunk = mean(flux[idx_out])
       chunk_of_spectrum = ChunkOfSpectrum{eltype(λ),eltype(flux),eltype(var)}(view(λ,idx_out), view(flux,idx_out), view(var,idx_out))
       Δv = 400  # TODO: Figure out good choice for Δv
       Δλ = Δv /RvSpectML.speed_of_light_mps*λc[c]
       if mean_flux_in_chunk > 0
           if alg == :Linear
               interp_chunk_to_grid_linear!(view(flux_matrix,:,c), view(var_matrix,:,c), chunk_of_spectrum, chunk_grid)
           elseif alg == :Sinc
               chunk_grid = range(λc[c]-(num_λ_per_chunk-1)//2*Δλ, stop=λc[c]+(num_λ_per_chunk-1)//2*Δλ, length=num_λ_per_chunk )
               interp_chunk_to_grid_sinc!(view(flux_matrix,:,c), view(var_matrix,:,c), chunk_of_spectrum, chunk_grid, Filter=sinc_filter)
           elseif alg == :GP
               interp_chunk_to_grid_gp!(view(flux_matrix,:,c), view(var_matrix,:,c), chunk_of_spectrum, chunk_grid)
           elseif alg == :TemporalGP
               interp_chunk_to_grid_gp_temporal!(view(flux_matrix,:,c), view(var_matrix,:,c), chunk_of_spectrum, chunk_grid , use_logx=true, use_logy=false, smooth_factor=smooth_factor)
           end
       else  # no flux in chunk?!?
           @warn "No flux in chunk " * string(c) * " at time " * string(t) * "."
       end
       idx_start += length(chunk_map[c])
   end # for c
   return flux_matrix
    #flush(stdout)
    #println("# mean flux pre-normalize = ", mean(mean_flux))
    #=
    idx_start = 0
    for c in 1:length(chunk_grids)
        idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
        mean_flux_in_chunk = mean(mean_flux[idx])
        #println("Normalizing chunk ", c, " (",idx,") by ", mean_flux_in_chunk)
        #dmeanfluxdlnλ[idx] .= calc_dfluxdlnlambda(vec(sum(view(flux_matrix,idx,:),dims=2)),vec(sum(view(var_matrix,idx,:), dims=2)) ) ./mean_flux_in_chunk
        mean_flux[idx] ./= mean_flux_in_chunk
        mean_var[idx] ./= mean_flux_in_chunk^2
        #dfluxdlnλ[idx] ./= mean_flux_in_chunk
        #d2fluxdlnλ2[idx] ./= mean_flux_in_chunk
        dfluxdlnλ[idx] .= calc_dfluxdlnlambda(mean_flux[idx],view(λ_vec,idx))
        d2fluxdlnλ2[idx] .= calc_d2fluxdlnlambda2(mean_flux[idx],view(λ_vec,idx))

        #dmeanfluxdlnλ[idx] .= calc_dfluxdlnlambda(vec(sum(view(flux_matrix,idx,:),dims=2)),view(λ_vec,idx))  ./mean_flux_in_chunk
        idx_start += length(chunk_grids[c])
    end
    =#
    #println("# mean flux post-normalize = ", mean(mean_flux))
    #flush(stdout)
    #mean_flux /= num_obs
    #mean_var .*= sqrt(oversample_factor)
    #dmeanfluxdlnλ ./= num_obs
    #return flux ( matrix=SpectralTimeSeriesCommonWavelengths(λ_vec,flux_matrix,var_matrix,chunk_map, Generic1D() ), mean_flux=mean_flux, mean_var=mean_var, deriv=dfluxdlnλ, deriv2=d2fluxdlnλ2 )
end
