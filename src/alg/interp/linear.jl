"""
Module with wrappers for performing linear interpolation..

Author: Eric Ford
"""
module LinearInterpolation

export make_interpolator_linear_flux, make_interpolator_linear_var

using Interpolations
using ..RvSpectML

""" Return interpolator for fluxes in spectra. """

function make_interpolator_linear_flux(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChunkOfSpectrum}
    Interpolations.LinearInterpolation(spectra.λ, spectra.flux)
end
#; boost_λ::Real = 1
""" Return interpolator for variances in spectra. """
function make_interpolator_linear_var(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChunkOfSpectrum}
    Interpolations.LinearInterpolation(spectra.λ, spectra.var)
end

end  # module
