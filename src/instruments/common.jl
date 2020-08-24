# Trait-based functions defaults (can be overwritten by instrument-specific versions)
orders_all(inst::AbstractInstrument2D) = min_order(inst):max_order(inst)
pixels_all(inst::AbstractInstrument2D) = min_pixels_in_order(inst):max_pixels_in_order(inst)
pixels_all(inst::AbstractInstrument1D) = min_pixel(inst):max_pixel(inst)
max_pixels_in_spectra(inst::AbstractInstrument1D) = length(pixels_all(inst))
max_pixels_in_spectra(inst::AbstractInstrument2D) = (max_order(inst)-min_order(inst)+1) * (max_pixel_in_order(inst)-min_pixel_in_order(inst)+1)
min_pixels_in_chunk(inst::AbstractInstrument1D) = 6

using CSV

"""Read manifest containing filename, bjd, target, and optionally additional metadata from CSV file. """
function read_manifest(fn::String)
    CSV.read(fn,threaded=false)
end

""" Read header from FITS file and return Dict with contents. """
function read_header(fn::String; header_idx::Integer = 1)
    f = FITS(fn)
    @assert 1<=header_idx<=length(f)
    #@assert read_key(f[1],"SKY-OBJ")[1] == "Solar"
    hdr = read_header(f[header_idx])
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
    hdr = read_header(f[1])
    # Check that header has all expected fields
    for field in fields_str
        @assert findfirst(isequal(field),keys(hdr)) != nothing
    end
    values = map(s->hdr[s],fields_str)
    df = Dict{Symbol,Any}(zip(fields,values))
    return df
end

""" Read mask in ESPRESSO csv format.
   ESPRESSO format: two columns, lambda and weight.
"""
function read_mask_espresso(fn::String)
    CSV.read(fn,threaded=false,header=["lambda","weight"],silencewarnings=true)
end

""" Read mask in VALD csv format.
   VALD format: lambda_lo, lambdaa_hi and weight.
 """
function read_mask_vald(fn::String)
    CSV.read(fn,threaded=false,header=["lambda_lo","lambda_hi","depth"])
end
