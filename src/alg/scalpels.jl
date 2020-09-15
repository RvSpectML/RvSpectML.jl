"""
Module for performing Scalpels (Self-Correlation Analysis of Line Profiles for Extracting Low-amplitude Shifts)
based on a CCF timeseries.

For algorithm information, see Collier-Cameron, Ford, Shahaf et al. 2020

Author: Eric Ford
Date:   September 2020
"""
module Scalpels

using Statistics, StatsBase
using LinearAlgebra
#using MultivariateStats
#using Plots, ColorSchemes

export clean_rvs_scalpels, scalpels_rms_vs_num_basis

function clean_rvs_scalpels(rvs::AbstractVector{T1}, ccfs::AbstractArray{T2,2}
                ; σ_rvs::AbstractVector{T3} = ones(length(rvs)),
                num_basis::Integer = 3 ) where { T1<:Real, T2<:Real, T3<:Real }
    @assert length(rvs) == length(σ_rvs)
    #@assert length(v_grid) == size(ccfs,1)
    @assert length(rvs) == size(ccfs,2)
    if num_basis == 0   return rvs   end
    @assert 0 <= num_basis < length(rvs)
    mean_rv = mean(rvs, weights(1.0 ./ σ_rvs.^2 ))
    rvs_centered = rvs .- mean_rv

    acfs = autocor(ccfs,0:size(ccfs,1)-1)
    acfs_minus_mean = acfs .- mean(acfs,dims=2)
    #Δv_grid = convert(Float64,v_grid.step).*(0:size(acfs,1)-1)
    svd_acfs = svd(acfs_minus_mean')
    alpha = svd_acfs.U'*rvs_centered

    idx = sortperm(abs.(alpha),rev=true)
    U_keep = view(svd_acfs.U,:,idx[1:num_basis])
    Δrv_shape = U_keep*U_keep'*rvs_centered
    rvs_clean = rvs .- Δrv_shape
    return rvs_clean
end

function scalpels_rms_vs_num_basis(rvs::AbstractVector{T1}, ccfs::AbstractArray{T2,2}
                ; σ_rvs::AbstractVector{T3} = ones(length(rvs)),
                max_num_basis::Integer = length(rvs)-1 ) where { T1<:Real, T2<:Real, T3<:Real }
    rms_scalels = map(b->std(clean_rvs_scalpels(rvs,ccfs,num_basis=b)), 0:max_num_basis)
    #scatter(0:max_num_basis,rms_scalels,xlabel="# Basis vectors",ylabel="RMS RV (m/s)", label=:none)
end

end # module Scalpels
