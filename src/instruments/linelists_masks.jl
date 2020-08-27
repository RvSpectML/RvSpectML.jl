using DataFrames, CSV

abstract type AbstractCalcChunkWidth end

""" Functor to return a constant scale factor for chunk widths """
struct ChunkWidthFixedΔlnλ <: AbstractCalcChunkWidth
   Δlnλ::Float64
end
(f::ChunkWidthFixedΔlnλ)(λ::Real) = return f.Δlnλ

#function (f::ChunkWidthFixedΔlnλ)(λ::Real)
#    return Δlnλ
#end

"""
   Functor to return a scale factor for chunk widths that depends on the chunk's central wavelength.

# TODO: Replace with something useful, e.g., uses actual PSF behavior as function of wavelength
"""
struct ChunkWidthFuncOfλ <: AbstractCalcChunkWidth
   Δlnλ::Float64
end

function (f::ChunkWidthFuncOfλ)(λ::Real; chunk_size_factor::Real = default_chunk_size_factor,
    line_width::Real = default_line_width, min_chunk_Δv::Real = default_min_chunk_Δv )
    chunk_half_width = max(chunk_size_factor*line_width, min_chunk_Δv)/speed_of_light_mps
    psf_width = 1e-6*(1+(λ-6000)/60000)   # TODO: Repalce, just some random number for demo purposes
    Δlnλ = sqrt(chunk_half_width^2+psf_width^2)
    return Δlnλ
end

# Default values shared across instruments
default_chunk_size_factor = 3       # TODO: Figure out what value to use.  Ask Alex
default_min_chunk_Δv = 20          # km/s
default_line_width = RvSpectML.predict_line_width(5780,v_rot=1.8)  # km/s
default_calc_chunk_width = ChunkWidthFixedΔlnλ(default_chunk_size_factor*RvSpectML.predict_line_width(5780,v_rot=1.8)/speed_of_light_mps)

"""
    Read mask in ESPRESSO csv format.
ESPRESSO format: lambda and weight.
Warning: ESPRESSO masks don't provide line depth and sometimes include one entry for a blend of lines.
"""
function read_mask_espresso(fn::String; calcΔ::CCWT = default_calc_chunk_width) where { CCWT<:AbstractCalcChunkWidth }
    df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","weight"],delim=' ',ignorerepeated=true)
    df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda])
    Δ = calcΔ.(df[!,:lambda])
    df[!,:lambda_lo] = df[!,:lambda]./(1 .+ Δ)
    df[!,:lambda_hi] = df[!,:lambda].*(1 .+ Δ)
    df[!,:depth] = df[!,:weight]  # TODO: Decide out what we want to do about tracking depths and weights sepoarately
    return df
end

""" Read mask in VALD csv format.
   VALD format: lambda_lo, lambdaa_hi and depth.
"""
function read_mask_vald(fn::String; calcΔ::CCWT = default_calc_chunk_width) where { CCWT<:AbstractCalcChunkWidth }
    df = CSV.read(fn,DataFrame,threaded=false,header=["lambda_lo","lambda_hi","depth"])
    df[!,:lambda_lo] .= λ_air_to_vac.(df[!,:lambda_lo])
    df[!,:lambda_hi] .= λ_air_to_vac.(df[!,:lambda_hi])
    df[!,:lambda] = sqrt.(df[!,:lambda_lo].*df[!,:lambda_hi])
    Δ = calcΔ.(df[!,:lambda])
    df[!,:lambda_lo] = df[!,:lambda]./(1 .+ Δ)
    df[!,:lambda_hi] = df[!,:lambda].*(1 .+ Δ)
    df[!,:weight] = df[!,:depth] # TODO: Decide out what we want to do about tracking depths and weights sepoarately
    return df
end

""" Return indices of any chunks in df that have overlapping lambda_hi[i] and lambda_lo[i+1].  """
function find_overlapping_chunks(df::DataFrame; verbose::Bool = true)
    @assert hasproperty(df,:lambda_lo)
    @assert hasproperty(df,:lambda_hi)
   if ! any(df.lambda_hi[1:end-1] .>= df.lambda_lo[2:end])   # Is there any point to keeping the if return early bit?
       return
   end
   idx_overlap = findall(df.lambda_hi[1:end-1] .>= df.lambda_lo[2:end])
   if verbose
       println("# Number of overlapping chunks: ",length(idx_overlap))
   end
   return idx_overlap
end
