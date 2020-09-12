"""
    Code to loading mask shapes for CCF
Author: Eric Ford
Created: September 2019
"""

"A struct implementing a specific mask shapes should be a subtype of AbstractCCFMaskShape."
abstract type AbstractCCFMaskShape  end

default_gaussian_ccf_σ = 5000.0
default_supergaussian_ccf_fwhm = 448.0*4.5

include("tophat.jl")
include("gaussian.jl")
include("halfcos.jl")
include("supergaussian.jl")
