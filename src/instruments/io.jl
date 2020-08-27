using DataFrames, CSV, FITSIO

"""Read manifest containing filename, bjd, target, and optionally additional metadata from CSV file. """
function read_manifest(fn::String)
    df = CSV.read(fn,DataFrame,threaded=false)
    @assert hasproperty(df,:filename)
    @assert hasproperty(df,:bjd)
    @assert hasproperty(df,:target)
    @assert size(df,1) >= 1
    return df
end

""" Read header from FITS file and return Dict with contents. """
function read_header(fn::String; header_idx::Integer = 1)
    println("# Reading: ",fn, " hdu= ",header_idx)
    f = FITS(fn)
    @assert 1<=header_idx<=length(f)
    #@assert read_key(f[header_idx],"SKY-OBJ")[1] == "Solar"
    hdr = FITSIO.read_header(f[header_idx])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
end


""" Read metadata in FITS header and return data for keys in fields_str/fields as a Dict. """
function read_metradata_from_fits
end

function read_metadata_from_fits(fn::String, fields::Array{Symbol,1})
    fields_str=string.(fields)
    read_metadata_from_fits(fn,fields=fields,fields_str=fields_str)
end

function read_metadata_from_fits(fn::String, fields_str::AbstractArray{AS,1} )  where { AS<:AbstractString }
    fields = map(f->Symbol(f),fields_str)
    read_metadata_from_fits(fn,fields=fields,fields_str=fields_str)
end

function read_metadata_from_fits(fn::String; fields::Array{Symbol,1}, fields_str::AbstractArray{AS,1} )  where { AS<:AbstractString }
    @assert length(fields) == length(fields_str)
    f = FITS(fn)
    hdr = FITSIO.read_header(f[1])
    # Check that header has all expected fields
    for field in fields_str
        @assert findfirst(isequal(field),keys(hdr)) != nothing
    end
    values = map(s->hdr[s],fields_str)
    df = Dict{Symbol,Any}(zip(fields,values))
    return df
end
