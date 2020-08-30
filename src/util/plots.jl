"""
Convenience code for plotting spectra, associated data types and results of analysis

Author: Eric Ford
Created: August 2020
"""

using Plots
export plot_spectrum_chunks

""" plot_spectrum_chunks( chunklist, chunks )
Returns a plot based of flux versus wavelength for selected chunks from the spectra in chunklist.
"""
function plot_spectrum_chunks(clt::CLT, chunks::Union{Integer,AbstractUnitRange};
    time_idx::Union{Integer,AbstractUnitRange}=1, color = nothing,
    plt::PT = Plots.plot(legend=:none)) where {CLT<:AbstractChunkListTimeseries, PT<:AbstractPlot}
    xmin = Inf # minimum(order_list_df.lambda_lo[chunk_idx])
    xmax = 0   # maximum(order_list_df.lambda_hi[chunk_idx])
    for c in chunks
        t = 1
        if color == nothing
            plot!(plt,clt.chunk_list[t].data[c].λ ,clt.chunk_list[t].data[c].flux)
        else
            colorval = typeof(color)<:Union{Integer,Symbol} ? color : color[c]
            plot!(plt,clt.chunk_list[t].data[c].λ ,clt.chunk_list[t].data[c].flux, linecolor=colorval)
        end
        xmin = min(xmin,minimum(clt.chunk_list[t].data[c].λ))
        xmax = max(xmax,maximum(clt.chunk_list[t].data[c].λ))
    end
    xlabel!(plt,"λ (Å)")
    ylabel!(plt,"Flux)")
    #xlims!(xmin,xmax)
    plt
end

"""" get_λs( grids::Vector{AbstractRange}, r::AbstractUnitRange )
Return range of wavelengths corresponding to indices in r from a vector of ranges.
Useful for extracting wavelengths from a chunklist.
"""
function get_λs(grids::V1, r::AR) where { T1<:AbstractRange, V1<:AbstractVector{T1}, AR<:AbstractUnitRange }
  local idx_offset = 1
  #println("first(r) = ",first(r), "   last(r) = ",last(r) )
  for chunk_range in grids
    #println("c: ",c," chunk_range=", chunk_range, " idx_offset=",idx_offset," - ", idx_offset+length(chunk_range)-1)
    if (first(r)>=idx_offset) && (last(r)<=idx_offset+length(chunk_range)-1)
       idx = (first(r)-idx_offset+1):(last(r)-idx_offset+1)
       return chunk_range[idx]
    else
      idx_offset += length(chunk_range)
    end
  end
  @error "get_λs: Failed to find range (",r,")."
  return
end

# Move these to PCA file/module?

function plot_basis_vectors(λ::VR1, f_mean::V1, deriv::V2, proj::A3;
                              num_basis::Integer = 4, idx_plt::AR = 1:length(f_mean) ) where {
                              R1<:AbstractRange, VR1<:AbstractVector{R1}, T1<:Real, V1<:AbstractVector{T1},
                              T2<:Real, V2<:AbstractVector{T2}, T3<:Real, A3<:AbstractArray{T3,2}, AR<:AbstractUnitRange }
  @assert num_basis >=1
  @assert sum(length.(λ)) == length(f_mean)
  @assert length(f_mean) == length(deriv)
  @assert size(proj,1) == length(f_mean)
  @assert minimum(idx_plt) >= 1
  @assert maximum(idx_plt) <= length(f_mean)
  if num_basis != 4   @warn "plot_basis_vectors hasn't been generalized for num_basis!=4."   end
  λ_plt = get_λs(λ,idx_plt)
  local plt0 = plot(λ_plt, (f_mean[idx_plt].-1.0)./std(f_mean[idx_plt]),linecolor=:black, label="Std Mean")
  plt0 = plot!(λ_plt, deriv[idx_plt]./std(deriv[idx_plt]),linecolor=:red,label="Std Deriv")
  plt0 = plot!(λ_plt, proj[idx_plt,1]./std(proj[idx_plt,1]),linecolor=:blue,label=:none)
  plt0 = plot!(λ_plt, proj[idx_plt,2]./std(proj[idx_plt,2]),linecolor=:green,label=:none)
  local plt1 = plot(λ_plt, proj[idx_plt,1],linecolor=:blue,label=:none) #"PC 1")
  ylabel!(plt1,"PC1")
  local plt2 = plot(λ_plt, proj[idx_plt,2],linecolor=:green,label=:none) #"PC 2")
  ylabel!(plt2,"PC2")
  local plt3 = plot(λ_plt, proj[idx_plt,3],linecolor=:cyan,label=:none) #"PC 3")
  ylabel!(plt3,"PC3")
  local plt4 = plot(λ_plt, proj[idx_plt,4],linecolor=:magenta,label=:none) #"PC 4")
  ylabel!(plt4,"PC4")
  local pltall = plot(plt0,plt1,plt2,plt3,plt4, layout = (5,1) )
  xlabel!(pltall,"λ (Å)")
  pltall
end

function plot_basis_scores(times::V1, rvs::V2, scores::A3;
                              num_basis::Integer = 4, idx_plt::AR = 1:length(times) ) where {
                              T1<:Real, V1<:AbstractVector{T1}, T2<:Real, V2<:AbstractVector{T2}, T3<:Real, A3<:AbstractArray{T3,2}, AR<:AbstractUnitRange }
  @assert length(times) >=1
  @assert length(times) == length(rvs)
  @assert length(times) == size(scores,2)
  @assert size(scores,1) >= num_basis
  @assert minimum(idx_plt) >= 1
  @assert maximum(idx_plt) <= length(times)
  if num_basis != 4   @warn "plot_basis_scores hasn't been generalized for num_basis!=4."   end
  local plt0 = scatter(times[idx_plt],rvs[idx_plt],label=:none,color=:black)
  ylabel!(plt0,"RV (m/s)")
  local plt1 = scatter(times[idx_plt],scores[1,idx_plt],label=:none,color=:blue)
  local plt2 = scatter(times[idx_plt],scores[2,idx_plt],label=:none,color=:green)
  local plt3 = scatter(times[idx_plt],scores[3,idx_plt],label=:none,color=:cyan)
  local plt4 = scatter(times[idx_plt],scores[4,idx_plt],label=:none,color=:magenta)
  local pltall = plot(plt0,plt1,plt2,plt3,plt4, layout = (5,1) )
  xlabel!(pltall,"Time")
  pltall
end


function plot_basis_scores_cor(rvs::V2, scores::A3;
                              num_basis::Integer = 4, idx_plt::AR = 1:length(rvs) ) where {
                              T2<:Real, V2<:AbstractVector{T2}, T3<:Real, A3<:AbstractArray{T3,2}, AR<:AbstractUnitRange }
  @assert length(rvs) >=1
  @assert length(rvs) == size(scores,2)
  @assert size(scores,1) >= num_basis
  @assert minimum(idx_plt) >= 1
  @assert maximum(idx_plt) <= length(rvs)
  if num_basis != 4   @warn "plot_basis_scores_cor hasn't been generalized for num_basis!=4."   end
  local plt1 = scatter(rvs[idx_plt],vec(scores[1,idx_plt]),xlabel="RV",ylabel="PC1",legend=:none)
  local plt2 = scatter(rvs[idx_plt],vec(scores[2,idx_plt]),xlabel="RV",ylabel="PC2",legend=:none)
  local plt3 = scatter(rvs[idx_plt],vec(scores[3,idx_plt]),xlabel="RV",ylabel="PC3",legend=:none)
  local plt4 = scatter(vec(scores[1,idx_plt]),vec(scores[2,idx_plt]),xlabel="PC1",ylabel="PC2",legend=:none)
  local plt5 = scatter(vec(scores[1,idx_plt]),vec(scores[3,idx_plt]),xlabel="PC1",ylabel="PC3",legend=:none)
  local plt6 = scatter(vec(scores[2,idx_plt]),vec(scores[3,idx_plt]),xlabel="PC2",ylabel="PC3",legend=:none)
  local plt7 = scatter(vec(scores[1,idx_plt]),vec(scores[4,idx_plt]),xlabel="PC2",ylabel="PC4",legend=:none)
  local plt8 = scatter(vec(scores[2,idx_plt]),vec(scores[4,idx_plt]),xlabel="PC2",ylabel="PC4",legend=:none)
  local plt9 = scatter(vec(scores[3,idx_plt]),vec(scores[4,idx_plt]),xlabel="PC3",ylabel="PC4",legend=:none)
  local pltall = plot(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8,plt9,layout=(3,3))
  pltall
end
