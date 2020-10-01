"""
   `interp_chunk_to_grid_gp_temporal!( flux_out, var_out, chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- flux_out: (results stored into this array)
- var_out: (results stored into this array)
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
- boost_factor: divide wavelengths by boost_factor
Optional Arguments:
- Filter: Vector with pre-allocated workspace (if length>=1)
# Returns
- flux_out
"""
function interp_chunk_to_shifted_grid_gp_temporal!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR, boost_factor::Real; use_logx::Bool = true,  use_logy::Bool = false, smooth_factor::Real=1 ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    y = eltype(chunk.flux) == Float64  ? chunk.flux : convert.(Float64,chunk.flux)
    var = eltype(chunk.var) == Float64  ? chunk.var : convert.(Float64,chunk.var)
    # TODO: Could reduce allocations with Lazy Arrays?
    mean_y = mean(y)
    if !use_logy
        y_to_use .= y_to_use./mean_y .- 1.0
        var_to_use .= var_to_use ./ mean_y^2
    else
        var_to_use .= var_to_use ./ y_to_use.^2
        y_to_use .= log.(y_to_use).-log(mean_y)
    end

    flux_out .= TemporalGPInterpolation.predict_mean(chunk.λ./boost_factor, y_to_use, grid, sigmasq_obs = var, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    # TODO: Update var_out to actually use the right GP
    var_out .= TemporalGPInterpolation.predict_mean(chunk.λ./boost_factor, var, grid, sigmasq_obs = var, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    if !use_logy
        flux_out .= mean_y.*(flux_out.+1.0)
        var_out .*= mean_y^2
    else
        flux_out .= exp.(flux_out.+log(mean_y))
        var_out .*= flux_out.^2
    end
    return flux_out
end

"""
   `interp_chunk_to_shifted_grid_gp_temporal( chunk_of_spectrum, wavelengths, boost_factor )
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
- boost_factor: divide wavelengths by boost_factor
Optional Arguments:
- Filter: Vector with pre-allocated workspace (if length>=1)
# Returns
- flux_out
"""
function interp_chunk_to_shifted_grid_gp_temporal( chunk::AC, grid::AR, boost_factor::Real; use_logx::Bool = true,  use_logy::Bool = false,smooth_factor::Real=1 ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp!(flux_out, var_out, chunk, grid, boost_factor=boost_factor, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor)
    return (flux=flux_out, var=var_out)
end

"""
   `interp_chunk_to_grid_gp_temporal( chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
# Returns
- flux_out
"""
function interp_chunk_to_grid_gp_temporal!( flux_out::AA1, var_out::AA2, chunk::AC, grid::AR ; use_logx::Bool = true,  use_logy::Bool = false, smooth_factor::Real=1 ) where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, AC<:AbstractChunkOfSpectrum, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    @assert size(flux_out) == size(var_out)
    @assert size(flux_out) == size(grid)
    y_to_use = eltype(chunk.flux) == Float64  ? copy(chunk.flux) : convert.(Float64,chunk.flux)
    var_to_use = eltype(chunk.var) == Float64  ? copy(chunk.var) : convert.(Float64,chunk.var)

    # TODO: Could reduce allocations with Lazy Arrays?
    mean_y = mean(y_to_use)
    if !use_logy
        y_to_use .= y_to_use./mean_y .- 1.0
        var_to_use .= var_to_use ./ mean_y^2
    else
        var_to_use .= var_to_use ./ y_to_use.^2
        y_to_use .= log.(y_to_use).-log(mean_y)
    end

    flux_out .= TemporalGPInterpolation.predict_mean(chunk.λ, y_to_use, grid, sigmasq_obs = var_to_use, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    # TODO: Update var_out to actually use the right GP
    var_out .= TemporalGPInterpolation.predict_mean(chunk.λ, var_to_use, grid, sigmasq_obs = var_to_use, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
    if !use_logy
        flux_out .= mean_y.*(flux_out.+1.0)
        var_out .*= mean_y^2
    else
        flux_out .= exp.(flux_out.+log(mean_y))
        var_out .*= flux_out.^2
    end
    return flux_out
end

"""
   `interp_chunk_to_grid_gp_temporal( chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
# Returns
- flux_out
"""
function interp_chunk_to_grid_gp_temporal( chunk::AC, grid::AR ; use_logx::Bool = true,  use_logy::Bool = false, smooth_factor::Real=1 ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_gp!(flux_out, var_out, chunk, grid, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor)
    return (flux=flux_out, var=var_out)
end
