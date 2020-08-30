"""
   IO functions for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

using DataFrames, FITSIO
using CSV, Interpolations

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest(data_path::String)
    dir_filelist = readdir(data_path,join=true)
    idx_spectra = map(fn->occursin(r"^neid\w+\.fits$", last(split(fn,'/')) ),dir_filelist)
    spectra_filelist = dir_filelist[idx_spectra]

    df_files = DataFrame(Filename = String[], target = String[], bjd = Float64[], ssbz=Float64[] )
    map(fn->add_metadata_from_fits!(df_files,fn),spectra_filelist)
    df_files
end

"""Create Dict containing filename and default metadata from file."""
function read_metadata(fn::String)
    dict = read_metadata_from_fits(fn,fields=metadata_symbols_default(NEID2D()),fields_str=metadata_strings_default(NEID2D()))
    dict[:Filename] = fn
    return dict
end

function add_metadata_from_fits!(df::DataFrame, fn::String)
    metadata_keep = read_metadata(fn)
    #if metadata_keep[:target] != "Solar" return df end
    push!(df, metadata_keep)
    return df
end

""" Read NEID (non-solar) data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] != "Solar"
    λ, flux, var  = FITSIO.read(f["SCIWAVE"]), FITSIO.read(f["Sci Flux"]), FITSIO.read(f["Sci Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_data(fn::String)
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] != "Solar"
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    λ, flux, var  = FITSIO.read(f["SCIWAVE"]), read(f["Sci Flux"]), read(f["Sci Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata)
end

""" Read NEID Solar data from FITS file, and return in a Spectra2DBasic object."""
function read_solar_data
end

function read_solar_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] == "Solar"
    λ, flux, var  = read(f["SKYWAVE"]), read(f["Sky Flux"]), read(f["Sky Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_solar_data(fn::String)
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] == "Solar"
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    λ, flux, var  = FITSIO.read(f["SKYWAVE"]), FITSIO.read(f["Sky Flux"]), FITSIO.ead(f["Sky Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_solar_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_solar_data(fn,metadata)
end


""" Read CSV of NEID drift corrections, interpolate to bjd's in df and insert into df[:,drift]. """
function read_drift_corrections!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    drift_corrections = CSV.read(fn, DataFrame, header=["bjd", "sci_drift", "cal_drift"]);
    @assert any(isequal(:bjd),propertynames(drift_corrections))
    @assert any(isequal(:cal_drift),propertynames(drift_corrections))
    drift_interp = LinearInterpolation(drift_corrections[!,:bjd],drift_corrections[!,:cal_drift])
    df[!,:drift] = -drift_interp.(df[!,df_time_col])
    return df
end

""" Read CSV of NEID barycentric corrections, interpolate to bjd's in df and insert into df[:,ssb_rv]. """
function read_barycentric_corrections!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    ssb_corrections = CSV.read(fn, DataFrame, header=["bjd","rv_ssb"], select=[1,2], types=types=[Float64,Float64], datarow=2, silencewarnings=true);
    @assert any(isequal(:bjd),propertynames(ssb_corrections))
    @assert any(isequal(:rv_ssb),propertynames(ssb_corrections))
    ssb_interp = LinearInterpolation(ssb_corrections.bjd, ssb_corrections.rv_ssb)
    df[!,:ssb_rv] = ssb_interp.(df[!,df_time_col])
    return df
end

""" Read space delimited file with differential extinction corrections, interpolate to bjd's in df and insert into df[:,diff_ext_rv]. """
function read_differential_extinctions!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    df_diff_ext = CSV.read(fn, DataFrame, delim='\t')
    @assert any(isequal(:JD),propertynames(df_diff_ext))
    @assert any(isequal(:delta_vr),propertynames(df_diff_ext))
    diff_ext = LinearInterpolation(df_diff_ext.JD, df_diff_ext.delta_vr)
    df[!,:diff_ext_rv] = -diff_ext(df[!,df_time_col])/100 # cm/s -> m/s
    return df
end
