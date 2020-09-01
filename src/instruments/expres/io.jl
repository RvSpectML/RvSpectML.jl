"""
   IO functions for the EXPRES spectrograph
      http://exoplanets.astro.yale.edu/expresBlog/
      https://ui.adsabs.harvard.edu/abs/2016SPIE.9908E..6TJ/abstract
Author: Eric Ford
Created: August 2020
"""

using DataFrames, FITSIO
using CSV, Interpolations

have_issued_expres_bjd_warning = false

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest(data_path::String)
    dir_filelist = readdir(data_path,join=true)
    idx_spectra = map(fn->occursin(r"^\d+_\d+\.\d+\.fits$", last(split(fn,'/')) ),dir_filelist)
    spectra_filelist = dir_filelist[idx_spectra]
    @assert length(spectra_filelist) >= 1
    df_files = DataFrame(read_metadata(spectra_filelist[1]))
    if length(spectra_filelist) >= 2
        map(fn->add_metadata_from_fits!(df_files,fn),spectra_filelist[2:end])
    end
    df_files
end

function jd2mjd(jd::Real)
    @assert jd > 2400000.5  # There's no EPRV data from that long ago!
    mjd = jd - 2400000.5
    return mjd
end


"""Create Dict containing filename and default metadata from file."""
function read_metadata(fn::String)
    dict1 = read_metadata_from_fits(fn,hdu=1,fields=metadata_symbols_default(EXPRES2D()),fields_str=metadata_strings_default(EXPRES2D()))
    global have_issued_expres_bjd_warning
    if !have_issued_expres_bjd_warning
        @warn "Currently, bjd field contains modified julian date of geometric midpoint of exposure for EXPRES observations, and is NOT corrected to be at solar system barycenter.\n  We need to fix this at some point."
        have_issued_expres_bjd_warning = true
    end
    dict1[:bjd] = jd2mjd(datetime2julian(DateTime(dict1[:midpoint])))
    dict1[:Filename] = fn
    dict2 = read_metadata_from_fits(fn,hdu=2,fields=metadata_hdu2_symbols_default(EXPRES2D()),fields_str=metadata_hdu2_strings_default(EXPRES2D()))
    dict = merge(dict1,dict2)
    return dict
end

function add_metadata_from_fits!(df::DataFrame, fn::String)
    metadata_keep = read_metadata(fn)
    push!(df, metadata_keep)
    return df
end

""" Read EXPRES data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(f::FITS, metadata::Dict{Symbol,Any} )
    λ, spectrum, uncertainty = FITSIO.read(f["optimal"],"bary_wavelength"), FITSIO.read(f["optimal"],"spectrum"), FITSIO.read(f["optimal"],"uncertainty")
    # For pixels where a presumably more accurate wavelength is avaliable, overwrite it.
    excalibur_mask, λ_excalibur = FITSIO.read(f["optimal"],"excalibur_mask"), FITSIO.read(f["optimal"],"bary_excalibur")
    λ[excalibur_mask] .= λ_excalibur[excalibur_mask]
    # If/when we want to read in additional information avaliable from EXPRES pipeline
    blaze = FITSIO.read(f["optimal"],"blaze")
    # blaze, continuum, tellurics = FITSIO.read(f["optimal"],"blaze"), FITSIO.read(f["optimal"],"continuum"), FITSIO.read(f["optimal"],"tellurics")
    # Restore fluxes to include the blaze function and scale uncertainties appropriately
    flux = spectrum.*blaze
    # Since EXPRES pipeline returns standard deviation rather than variance
    var = (uncertainty.*blaze).^2
    Spectra2DBasic(λ, flux, var, EXPRES2D(), metadata=metadata)
end

function read_data(fn::String, metadata = Dict{Symbol,Any}() )
    f = FITS(fn)
    hdr1 = FITSIO.read_header(f[1])
    hdr2 = FITSIO.read_header(f[2])
    metadata1 = Dict(zip(map(k->Symbol(k),hdr1.keys),hdr1.values))
    metadata2 = Dict(zip(map(k->Symbol(k),hdr2.keys),hdr2.values))
    metadata_from_data = Dict{Symbol,Any}()
    blaze, continuum, tellurics = FITSIO.read(f["optimal"],"blaze"), FITSIO.read(f["optimal"],"continuum"), FITSIO.read(f["optimal"],"tellurics")
    pixel_mask, excalibur_mask = FITSIO.read(f["optimal"],"pixel_mask"), FITSIO.read(f["optimal"],"excalibur_mask")
    metadata_from_data[:pixel_mask] = pixel_mask
    metadata_from_data[:excalibur_mask] = excalibur_mask
    metadata_from_data[:blaze] = blaze
    metadata_from_data[:continuum] = continuum
    metadata_from_data[:tellurics] = tellurics
    metadata_combo = merge(metadata,metadata1,metadata2,metadata_from_data)
    read_data(f,metadata_combo)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata)
end

function read_data_only(fn::String)
    f = FITS(fn)
    metadata = Dict{Symbol,Any}()
    read_data(f, metadata)
end
