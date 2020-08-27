using DataFrames, CSV

abstract type AbstractCalcChunkWidth end

""" Functor to return a constant scale factor for chunk widths (regardless of wavelength passed)"""
struct ChunkWidthFixedΔlnλ <: AbstractCalcChunkWidth
   Δlnλ::Float64
end
(f::ChunkWidthFixedΔlnλ)(λ::Real) = f.Δlnλ

"""
   Functor to return a scale factor for chunk widths that depends on the chunk's central wavelength.

# TODO: Replace with something useful, e.g., uses actual PSF behavior as function of wavelength
"""
struct ChunkWidthFuncOfλ <: AbstractCalcChunkWidth
   Δlnλ::Float64
end

function (f::ChunkWidthFuncOfλ)(λ::Real; chunk_size_factor::Real = default_chunk_size_factor,
    line_width::Real = default_line_width_kmps, min_chunk_Δv::Real = default_min_chunk_Δv )
    chunk_half_width = max(chunk_size_factor*line_width, min_chunk_Δv)*1000/speed_of_light_mps
    psf_width = λ*1e-6*(1+(λ-6000)/60000)   # TODO: Repalce, just some random number for demo purposes
    Δlnλ = sqrt(chunk_half_width^2+psf_width^2)
    return Δlnλ
end


"""
    Read mask in ESPRESSO csv format.
ESPRESSO format: lambda and weight.
Warning: ESPRESSO masks don't provide line depth and sometimes include one entry for a blend of lines.
"""
function read_mask_espresso(fn::String; calcΔ::CCWT = default_calc_chunk_width) where { CCWT<:AbstractCalcChunkWidth }
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","weight"],delim=' ',ignorerepeated=true)
    df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda])
    local Δ = calcΔ.(df[!,:lambda])
    @assert all(Δ.>0)
    df[!,:lambda_lo] = df[!,:lambda]./(1 .+ Δ)
    df[!,:lambda_hi] = df[!,:lambda].*(1 .+ Δ)
    df[!,:depth] = df[!,:weight]  # TODO: Decide out what we want to do about tracking depths and weights sepoarately
    return df
end

""" Read mask in VALD csv format.
   VALD format: lambda_lo, lambdaa_hi and depth.
"""
function read_mask_vald(fn::String; calcΔ::CCWT = default_calc_chunk_width) where { CCWT<:AbstractCalcChunkWidth }
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda_lo","lambda_hi","depth"])
    df[!,:lambda_lo] .= λ_air_to_vac.(df[!,:lambda_lo])
    df[!,:lambda_hi] .= λ_air_to_vac.(df[!,:lambda_hi])
    df[!,:lambda] = sqrt.(df[!,:lambda_lo].*df[!,:lambda_hi])
    local Δ = calcΔ.(df[!,:lambda])
    @assert all(Δ.>0)
    df[!,:lambda_lo] = df[!,:lambda]./(1 .+ Δ)
    df[!,:lambda_hi] = df[!,:lambda].*(1 .+ Δ)
    df[!,:weight] = df[!,:depth] # TODO: Decide out what we want to do about tracking depths and weights sepoarately
    return df
end

"""  Convert vacuum wavelength (in Å) to air wavelength
Ref: Donald Morton (2000, ApJ. Suppl., 130, 403) via
     https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
"""
function λ_vac_to_air(λ_vac::Real)
    @assert 3500 < λ_vac < 13000  # Making sure in Å for optical/NIR spectra.
    local s = 10000/λ_vac
    local n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2)
    return λ_vac/n
end

""" Convert air wavelength (in Å) to vacuum wavelength
Ref: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
     VALD3 tools use the following solution derived by N. Piskunov
"""
function λ_air_to_vac(λ_air::Real)
    @assert 3500 < λ_air < 13000  # Making sure in Å for optical/NIR spectra.
    local s = 10000/λ_air
    local n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^22) + 0.0001599740894897 / (38.92568793293 - s^22)
    λ_air*n
end

"""Estimate line width based on stellar Teff (K) and optionally v_rot (km/s).  Output in km/s."""
function predict_line_width(Teff::Real; v_rot::Real=zero(Teff))
    @assert 3000 < Teff < 10000 # K
    @assert 0 <= v_rot <=100 # km/s
    line_width_thermal = 13*sqrt(Teff/1e4) # km/s
    line_width = sqrt(v_rot^2+line_width_thermal^2) # km/s
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

""" Return DataFrame with specifications for each chunk which will contain one or more lines.
    Input:  line_list a DataFrame with:
    -  lambda_lo, lambda_hi, lambda, depth
    Output:  DataFrame with
    - lambda_lo & lambda_hi: boundaries for chunk
    - lambda & line_depths: arrays with info about each line
"""
function merge_lines(line_list::DataFrame)
    @assert hasproperty(line_list,:lambda_lo)
    @assert hasproperty(line_list,:lambda_hi)
    @assert hasproperty(line_list,:lambda)
    @assert hasproperty(line_list,:depth)
    chunk_list_df = DataFrame(:lambda_lo=>Float64[],:lambda_hi=>Float64[],
                                :lambda=>Array{Float64,1}[],:depth=>Array{Float64,1}[])
    num_lines = size(line_list,1)
    @assert num_lines >= 2
    lambda_lo_last = line_list[1,:lambda_lo]
    lambda_hi_last = line_list[1,:lambda_hi]
    lambda = [line_list[1,:lambda]]
    line_depths = [line_list[1,:depth]]
    for i in 2:num_lines
        lambda_lo = line_list[i,:lambda_lo]
        lambda_hi = line_list[i,:lambda_hi]
        if lambda_lo>lambda_hi_last
            push!(chunk_list_df, (lambda_lo_last, lambda_hi_last, lambda, line_depths))
            (lambda_lo_last, lambda_hi_last) = (lambda_lo, lambda_hi)
            lambda = [line_list[i,:lambda]]
            line_depths = [line_list[i,:depth]]
        else
            lambda_hi_last = lambda_hi
            push!(lambda,line_list[i,:lambda])
            push!(line_depths,line_list[i,:depth])
        end
    end
    if chunk_list_df[end,:lambda_hi] != lambda_hi_last
        #push!(chunk_list_df, (lambda_lo_last, lambda_hi_last))
        push!(chunk_list_df, (lambda_lo_last, lambda_hi_last, lambda, line_depths))
    end
    return chunk_list_df
end
