"""
   `interp_chunk_to_grid_linear!( flux_out, var_out, chunk_of_spectrum, wavelengths )`
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

"""
   `interp_chunk_to_shifted_grid_linear( chunk_of_spectrum, wavelengths, boost_factor )
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
function interp_chunk_to_shifted_grid_linear( chunk::AC, grid::AR, boost_factor::Real ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_shifted_grid_linear!(flux_out, var_out, chunk, grid, boost_factor)
    return (flux=flux_out, var=var_out)
end

"""
   `interp_chunk_to_grid_linear!( flux_out, var_out, chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- flux_out: (results stored into this array)
- var_out: (results stored into this array)
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
Optional Arguments:
- Filter: Vector with pre-allocated workspace (if length>=1)
# Returns
- flux_out
"""
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

"""
   `interp_chunk_to_grid_linear( chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using linear interpolation.
# Arguments:
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
# Returns
- flux_out
"""
function interp_chunk_to_grid_linear( chunk::AC, grid::AR ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_linear!(flux_out, var_out, chunk, grid)
    return (flux=flux_out, var=var_out)
end
