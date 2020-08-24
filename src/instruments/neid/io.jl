#using DataFrames, FITSIO

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest(data_path::String)
    dir_filelist = readdir(data_path,join=true)
    idx_spectra = map(fn->occursin(r"^neid\w+\.fits$", last(split(fn,'/')) ),dir_filelist)
    spectra_filelist = dir_filelist[idx_spectra]

    df_files = DataFrame(Filename = String[], target = String[], bjd = Float64[], ssbz=Float64[] )
    map(fn->add_metadata_from_fits!(df_files,fn),spectra_filelist)
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
function read_data
end

function read_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] != "Solar"
    λ, flux, var  = read(f["SCIWAVE"]), read(f["Sci Flux"]), read(f["Sci Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_data(fn::String)
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] != "Solar"
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    λ, flux, var  = read(f["SCIWAVE"]), read(f["Sci Flux"]), read(f["Sci Variance"])
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
    λ, flux, var  = read(f["SKYWAVE"]), read(f["Sky Flux"]), read(f["Sky Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_solar_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_neid_solar_data(fn,metadata)
end


#=
""" """
function (fn::String)
end
=#
