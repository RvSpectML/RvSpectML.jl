"""
   `interp_chunk_to_shifted_grid_sinc!( flux_out, var_out, chunk_of_spectrum, wavelengths, boost_factor )`
Return spectra interpolated onto a grid of points using sinc interpolation.
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

"""
   `interp_chunk_to_shifted_grid_sinc( chunk_of_spectrum, wavelengths, boost_factor )`
Return spectra interpolated onto a grid of points using sinc interpolation.
# Arguments:
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
- boost_factor: divide wavelengths by boost_factor
Optional Arguments:
- Filter: Vector with pre-allocated workspace (if length>=1)
# Returns
- flux_out
"""
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

"""
   `interp_chunk_to_grid_sinc!( flux_out, var_out, chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using sinc interpolation.
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

"""
   `interp_chunk_to_grid_sinc( chunk_of_spectrum, wavelengths )`
Return spectra interpolated onto a grid of points using sinc interpolation.
# Arguments:
- chunk_of_spectrum
- wavelengths: AbstractRange or AbstractArray of locations where chunk is to be interpolated to
# Returns
- flux_out
"""
function interp_chunk_to_grid_sinc( chunk::AC, grid::AR ) where {  AC<:AbstractChunkOfSpectrum, T2<:Real, AR<:Union{AbstractRange,AbstractArray{T2,1}} }
    flux_out = Array{Float64,1}(undef,length(grid))
    var_out = Array{Float64,1}(undef,length(grid))
    interp_chunk_to_grid_sinc!(flux_out, var_out, chunk, grid)
    return (flux=flux_out, var=var_out)
end
