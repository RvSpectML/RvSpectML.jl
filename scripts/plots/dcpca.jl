"""
Convenience code for plotting result of DCPCA

Author: Eric Ford
Created: August 2020
"""

using RvSpectML
using Plots

include("spectra.jl")

function plot_basis_vectors(λ::VR1, f_mean::V1, deriv::V2, proj::A3;
                              num_basis::Integer = 4, idx_plt::AR = 1:length(f_mean),
                              label = 0 ) where {
                              R1<:AbstractRange, VR1<:AbstractVector{R1}, T1<:Real, V1<:AbstractVector{T1},
                              T2<:Real, V2<:AbstractVector{T2}, T3<:Real, A3<:AbstractArray{T3,2}, AR<:AbstractUnitRange }
  @assert num_basis >=1
  @assert sum(length.(λ)) == length(f_mean)
  @assert length(f_mean) == length(deriv)
  @assert size(proj,1) == length(f_mean)
  @assert minimum(idx_plt) >= 1
  @assert maximum(idx_plt) <= length(f_mean)
  if num_basis > 4   @warn "plot_basis_vectors hasn't been generalized for num_basis!=4."   end
  λ_plt = get_λs(λ,idx_plt)
  λc = mean(λ_plt)
  plt0 = plot(λ_plt, (f_mean[idx_plt].-1.0)./maximum(abs.(f_mean[idx_plt].-1.0)),linecolor=:black, label=:none) #"Std Mean")
  plt0 = plot!(plt0,λ_plt, deriv[idx_plt]./maximum(abs.(deriv[idx_plt])),linecolor=:red,label=:none) #"Std Deriv")
  plt0 = plot!(plt0,[λc,λc],[-1.1,1.0], linecolor=:black, label=:none)
  ylims!(-1.1,1.1)
  pltall = [plt0]
  if num_basis >= 1
      plt0 = plot!(plt0,λ_plt, proj[idx_plt,1]./maximum(abs.(proj[idx_plt,1])),linecolor=:blue,label=:none)
      plt1 = plot(λ_plt, proj[idx_plt,1]./maximum(abs.(proj[idx_plt,1])),linecolor=:blue,label=:none) #"PC 1")
      plt1 = plot!(plt1,[λc,λc],[-1.1,1.0], linecolor=:black, label=:none)
      ylabel!(plt1,"PC1")
      ylims!(-1.1,1.1)
      push!(pltall,plt1)
  end
  if num_basis >= 2
      plt0 = plot!(plt0,λ_plt, proj[idx_plt,2]./maximum(abs.(proj[idx_plt,2])),linecolor=:green,label=:none)
      plt2 = plot(λ_plt, proj[idx_plt,2]./maximum(abs.(proj[idx_plt,2])),linecolor=:green,label=:none) #"PC 2")
      plt2 = plot!(plt2,[λc,λc],[-1.1,1.0], linecolor=:black, label=:none)
      ylabel!(plt2,"PC2")
      ylims!(-1.1,1.1)
      push!(pltall,plt2)
  end
  if num_basis >= 3
      plt3 = plot(λ_plt, proj[idx_plt,3]./maximum(abs.(proj[idx_plt,3])),linecolor=:cyan,label=:none) #"PC 3")
      plt3 = plot!(plt3,[λc,λc],[-1.1,1.0], linecolor=:black, label=:none)
      ylabel!(plt3,"PC3")
      ylims!(-1.1,1.1)
      push!(pltall,plt3)
  end
  if num_basis >= 4
      plt4 = plot(λ_plt, proj[idx_plt,4]./maximum(abs.(proj[idx_plt,4])),linecolor=:magenta,label=:none) #"PC 4")
      plt4 = plot!(plt4,[λc,λc],[-1.1,1.0], linecolor=:black, label=:none)
      ylabel!(plt4,"PC4")
      ylims!(-1.1,1.1)
      push!(pltall,plt4)
  end
  if num_basis == 1
      xlabel!(plt1,"λ (Å)  line# "  * string(label))
      pltall = plot(plt0,plt1, layout = (2,1) )
  elseif num_basis == 2
      xlabel!(plt2,"λ (Å)  line# "  * string(label))
      pltall = plot(plt0,plt1,plt2, layout = (3,1) )
  elseif num_basis == 3
      xlabel!(plt3,"λ (Å)  line# "  * string(label))
      pltall = plot(plt0,plt1,plt2,plt3, layout = (4,1) )
  elseif num_basis == 4
      xlabel!(plt4,"λ (Å)  line# "  * string(label))
      pltall = plot(plt0,plt1,plt2,plt3,plt4, layout = (5,1) )
  elseif num_basis == 5
      xlabel!(plt5,"λ (Å)  line# "  * string(label))
      pltall = plot(plt0,plt1,plt2,plt3,plt4,plt5, layout = (6,1) )
  end


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
  if num_basis > 4   @warn "plot_basis_scores hasn't been generalized for num_basis>4."   end
  local plt0 = scatter(times[idx_plt],rvs[idx_plt],label=:none,color=:black)
  ylabel!(plt0,"RV (m/s)")
  if num_basis >= 1
      plt1 = scatter(times[idx_plt],scores[1,idx_plt],label=:none,color=:blue)
      pltall = plot(plt0,plt1, layout = (2,1) )
  end
  if num_basis >= 2
      plt2 = scatter(times[idx_plt],scores[2,idx_plt],label=:none,color=:green)
      pltall = plot(plt0,plt1,plt2, layout = (3,1) )
  end
  if num_basis >= 3
      plt3 = scatter(times[idx_plt],scores[3,idx_plt],label=:none,color=:cyan)
      pltall = plot(plt0,plt1,plt2,plt3, layout = (4,1) )
  end
  if num_basis >= 4
      plt4 = scatter(times[idx_plt],scores[4,idx_plt],label=:none,color=:magenta)
      pltall = plot(plt0,plt1,plt2,plt3,plt4, layout = (5,1) )
  end
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
  if num_basis > 4   @warn "plot_basis_scores_cor hasn't been generalized for num_basis>4."   end
  if num_basis >= 1
      plt1 = scatter(rvs[idx_plt],vec(scores[1,idx_plt]),xlabel="RV",ylabel="PC1",legend=:none)
      pltall = plt1
  end
  if num_basis >= 2
      plt2 = scatter(rvs[idx_plt],vec(scores[2,idx_plt]),xlabel="RV",ylabel="PC2",legend=:none)
      plt4 = scatter(vec(scores[1,idx_plt]),vec(scores[2,idx_plt]),xlabel="PC1",ylabel="PC2",legend=:none)
      pltall = plot(plt1,plt2,plt4,layout=3)
  end
  if num_basis >= 3
      plt3 = scatter(rvs[idx_plt],vec(scores[3,idx_plt]),xlabel="RV",ylabel="PC3",legend=:none)
      plt5 = scatter(vec(scores[1,idx_plt]),vec(scores[3,idx_plt]),xlabel="PC1",ylabel="PC3",legend=:none)
      plt6 = scatter(vec(scores[2,idx_plt]),vec(scores[3,idx_plt]),xlabel="PC2",ylabel="PC3",legend=:none)
      pltall = plot(plt1,plt2,plt3,plt4,plt5,plt6,layout=(2,3))
  end
  if num_basis >= 4
      plt7 = scatter(vec(scores[1,idx_plt]),vec(scores[4,idx_plt]),xlabel="PC2",ylabel="PC4",legend=:none)
      plt8 = scatter(vec(scores[2,idx_plt]),vec(scores[4,idx_plt]),xlabel="PC2",ylabel="PC4",legend=:none)
      plt9 = scatter(vec(scores[3,idx_plt]),vec(scores[4,idx_plt]),xlabel="PC3",ylabel="PC4",legend=:none)
      pltall = plot(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8,plt9,layout=(3,3))
  end
  pltall
end
