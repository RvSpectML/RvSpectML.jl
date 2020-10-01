"""
   interp_chunk_to_grid_gp_brute_force!( flux_out, var_out, chunk_of_spectrum, wavelengths )
Return spectra interpolated onto a grid of points using Gaussian Process interpolation.
# Arguments:
- flux_out: (results stored into this array)
- var_out: (results stored into this array)
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
# Returns
- flux_out

NOTE:  Using own GP code for now, since include predicting derivatives and can minimize unnecessary dependancies.  We may need to
revisit this if we want improved speed
"""

function interp_chunk_to_shifted_grid_gp_brute_force!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR, boost_factor::Real ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
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

function interp_chunk_to_shifted_grid_gp_brute_force( chunk::AC, grid::AR, boost_factor::Real ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp_brute_force!(flux_out, var_out, chunk, grid, boost_factor=boost_factor)
    return (flux=flux_out, var=var_out)
end

function interp_chunk_to_grid_gp_brute_force!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
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

function interp_chunk_to_grid_gp_brute_force( chunk::AC, grid::AR ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp_brute_force!(flux_out, var_out, chunk, grid, )
    return (flux=flux_out, var=var_out)
end
