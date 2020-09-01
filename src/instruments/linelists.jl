"""
Code for creating, loading and maniuplating line lists and masks.

Author: Eric Ford
Created: August 2020
"""

using DataFrames, CSV

"""
    Read line list in ESPRESSO csv format.
ESPRESSO format: lambda and weight.
Warning: ESPRESSO masks don't provide line depth and sometimes include one entry for a blend of lines.
"""
function read_linelist_espresso(fn::String) #; calcΔ::CCWT = default_calc_chunk_width) where { CCWT<:AbstractCalcChunkWidth }
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","weight"],delim=' ',ignorerepeated=true)
    @assert hasproperty(df, :lambda)
    @assert hasproperty(df, :weight)
    df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda])
    return df
end

""" Read line list in VALD csv format.
   VALD format: lambda_lo, lambdaa_hi and depth.
"""
function read_linelist_vald(fn::String) # ; calcΔ::CCWT = default_calc_chunk_width) where { CCWT<:AbstractCalcChunkWidth }
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda_lo","lambda_hi","depth"])
    @assert hasproperty(df, :lambda_lo)
    @assert hasproperty(df, :lambda_hi)
    @assert hasproperty(df, :depth)
    df[!,:lambda_lo] .= λ_air_to_vac.(df[!,:lambda_lo])
    df[!,:lambda_hi] .= λ_air_to_vac.(df[!,:lambda_hi])
    df[!,:lambda] = sqrt.(df[!,:lambda_lo].*df[!,:lambda_hi])
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
    local n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) + 0.0001599740894897 / (38.92568793293 - s^2)
    λ_air*n
end
