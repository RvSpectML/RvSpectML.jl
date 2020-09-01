"""
Various utilities for manipulating spectra

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

"""
    apply_doppler_boost!(spectrum, doppler_factor) -> typeof(spectrum)
    apply_doppler_boost!(spectra, df) -> typeof(spectra)

Apply Doppler boost to spectra's λ's and update its metadata[:doppler_factor], so it will know how to undo the transform.
# Arguments:
* `spectrum::AbstractSpectra`: spectrum to be boosted
* `doppler_factor::Real`: boost factor (1 = noop)
or
* `spectra::AbstractArray{<:AbstractSpectra}`: spectra to be boosted
* `df::DataFrame`: provides `:drift` and `:ssb_rv` (in m/s) for calculating the Doppler boost for each spectrum

TODO: Improve documentation formatting.  This can serve as a template.
"""
function apply_doppler_boost! end

function apply_doppler_boost!(spectra::AS,doppler_factor::Real) where {AS<:AbstractSpectra}
    if doppler_factor == one(doppler_factor) return spectra end
    #println("# t= ",time, " doppler_factor= ",doppler_factor)
    spectra.λ .*= doppler_factor
    if hasproperty(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor
    end
    return spectra
end

function apply_doppler_boost!(spectra::AbstractArray{AS}, df::DataFrame ) where { AS<:AbstractSpectra }
    @assert size(spectra,1) == size(df,1)
    local doppler_factor = ones(size(spectra))
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :drift to apply."   end
    if  hasproperty(df,:drift)        doppler_factor .*= calc_doppler_factor.(df[!,:drift])          end
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :ssb_rv to apply."  end
    if  hasproperty(df,:ssb_rv)       doppler_factor   .*= calc_doppler_factor.(df[!,:ssb_rv])       end
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :diff_ext_rv to apply."  end
    if  hasproperty(df,:diff_ext_rv)  doppler_factor   .*= calc_doppler_factor.(df[!,:diff_ext_rv])  end
    map(x->apply_doppler_boost!(x[1],x[2]), zip(spectra,doppler_factor) );
end

"""  Calculate total SNR in (region of) spectra. """
function calc_snr(flux::AbstractArray{T1},var::AbstractArray{T2}) where {T1<:Real, T2<:Real}
    @assert size(flux) == size(var)
    sqrt(sum(flux.^2 ./ var))   # TODO: Test when weights are far from equal
end

function calc_snr(flux::Real,var::Real)
    flux / sqrt(var)   # TODO: Test when weights are far from equal
end

""" Calc normalization of spectra based on average flux in a ChunkList. """
function calc_normalization(chunk_list::ACL) where { ACL<:AbstractChunkList}
    total_flux = sum(sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end

""" Normalize spectrum, multiplying fluxes by scale_fac. """
function normalize_spectrum!(spectrum::ST, scale_fac::Real) where { ST<:AbstractSpectra }
    @assert 0 < scale_fac < Inf
    @assert !isnan(scale_fac^2)
    spectrum.flux .*= scale_fac
    spectrum.var .*= scale_fac^2
    return spectrum
end


""" Normalize each spectrum based on sum of fluxes in chunk_timeseries region of each spectrum. """
function normalize_spectra!(chunk_timeseries::ACLT, spectra::AS) where { ACLT<:AbstractChunkListTimeseries, ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    @assert length(chunk_timeseries) == length(spectra)
    for t in 1:length(chunk_timeseries)
        scale_fac = calc_normalization(chunk_timeseries.chunk_list[t])
        # println("# t= ",t, " scale_fac= ", scale_fac)
        normalize_spectrum!(spectra[t], scale_fac)
    end
    return chunk_timeseries
end
