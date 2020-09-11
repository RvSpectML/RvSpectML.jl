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

export clean_rvs_scalpels, make_plots_scalpels

function make_plots_scalpels(rvs::AbstractVector{T1}, ccfs::AbstractArray{T2,2}
                ; σ_rvs::AbstractVector{T3} = ones(length(rvs)),
                num_basis::Integer = 3,
                v_grid::AbstractRange = 1:size(ccfs,1),
                times::AbstractVector{T4} = collect(1:length(rvs)) ) where
                { T1<:Real, T2<:Real, T3<:Real, T4<:Real  }

    @assert length(rvs) == length(σ_rvs)
    @assert length(v_grid) == size(ccfs,1)
    @assert length(rvs) == size(ccfs,2)
    @assert 1 <= num_basis < length(rvs)

    mean_rv = mean(rvs, weights(1.0 ./ σ_rvs.^2 ))
    rvs_centered = rvs .- mean_rv

    ccfs_minus_mean = ccfs .- mean(ccfs,dims=2)
    heatmap(v_grid, 1:size(ccfs,2),ccfs_minus_mean')
    title!("CCF(v)-<CCF(v)>")
    xlabel!("v (m/s)")
    ylabel!("Obs #")
    savefig("ccf_heatmap.png")

    acfs = autocor(ccfs,0:size(ccfs,1)-1)
    acfs_minus_mean = acfs .- mean(acfs,dims=2)
    Δv_grid = convert(Float64,v_grid.step).*(0:size(acfs,1)-1)

    heatmap(Δv_grid, 1:size(acfs,2),acfs_minus_mean')
    title!("ACF(δv)-<ACF(δv)>")
    xlabel!("Δv (m/s)")
    ylabel!("Obs #")
    savefig("acf_heatmap.png")

    svd_ccfs = svd(ccfs_minus_mean')
    plot(v_grid,svd_ccfs.V[:,1:num_basis])
    title!("CCF(v) basis functions")
    xlabel!("v (m/s)")
    savefig("ccf_basis.png")

    svd_acfs = svd(acfs_minus_mean')
    plot(Δv_grid,svd_acfs.V[:,1:num_basis])
    title!("ACF(v) basis functions")
    xlabel!("Δv (m/s)")
    savefig("acf_basis.png")

    alpha = svd_acfs.U'*rvs_centered
    idx = sortperm(abs.(alpha),rev=true)
    U_keep = view(svd_acfs.U,:,idx[1:num_basis])
    Δrv_shape = U_keep*U_keep'*rvs_ccf_gauss
    rvs_clean = rvs .- Δrv_shape
    plt = scatter(times,rvs, label="RVs orig")
    scatter!(plt, times,rvs_clean, label="RVs cleaned", legend=:bottomright)
    xlabel!("Time (d)")
    ylabel!("RV (m/s)")
    title!("Cleaning RVs via Scalpels (" * string(num_basis) * " basis vectors)")
    savefig("scalpels_rvs_" * string(num_basis) * ".png")
    return rvs_clean
end

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
    alpha = svd_acfs.U'*rvs_centered

    idx = sortperm(abs.(alpha),rev=true)
    U_keep = view(svd_acfs.U,:,idx[1:num_basis])
    #Δrv_shape = svd_acfs.U[:,idx[1:num_basis]]*svd_acfs.U[:,idx[1:num_basis]]'*rvs_ccf_gauss
    Δrv_shape = U_keep*U_keep'*rvs_ccf_gauss
    rvs_clean = rvs .- Δrv_shape
    return rvs_clean
end

function scalpels_rms_vs_num_basis(rvs::AbstractVector{T1}, ccfs::AbstractArray{T2,2}
                ; σ_rvs::AbstractVector{T3} = ones(length(rvs)),
                max_num_basis::Integer = length(rvs)-1 ) where { T1<:Real, T2<:Real, T3<:Real }
    rms_scalels = map(b->std(clean_rvs_scalpels(rvs,ccfs,num_basis=b)), 0:max_num_basis)
    scatter(0:max_num_basis,rms_scalels,xlabel="# Basis vectors",ylabel="RMS RV (m/s)", label=:none)
end

end # module Scalpels
