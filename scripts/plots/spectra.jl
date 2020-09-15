"""
Convenience code for plotting spectra, associated data types and results of analysis

Author: Eric Ford
Created: August 2020
"""

using RvSpectML
using Plots

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

function add_time_gap_lines(plt::AbstractPlot, times::AbstractVector{T} ; Δt_boundaries::AbstractVector{T2} = [0.5, 1.5, 7, 14, 42],
                                    lweights::AbstractVector{T3} = [ 1, 1, 1, 1, 2, 2, 3, 3  ],
                                    lstyles::AbstractVector{Symbol} = [ :dot, :dashdot, :dash, :solid, :dash, :solid, :dash, :solid ]
                                          ) where { T<:Real, T2<:Real, T3<:Real  }
   @assert 2 <= length(Δt_boundaries) <= length(lstyles)-1
   Δt = times[2:end] .- times[1:end-1]
   xvals = [xlims(plt)[1], xlims(plt)[2]]
   for i in 1:(length(Δt_boundaries)-1)
      idx = findall(x -> Δt_boundaries[i]<=x<Δt_boundaries[i+1], Δt)
      if length(idx) >= 1
         map( y->plot!(plt,xvals, (y+0.5).*[1, 1],label=:none, color=:black, linewidth=lweights[i], linestyle=lstyles[i]), idx )
      end
   end
   idx = findall(x -> x > last(Δt_boundaries), Δt)
   if length(idx) >= 1
      map(y->plot!(plt,xvals, (y+0.5).*[1, 1],label=:none, color=:black, linewidth=lweights[length(idx)], linestyle=lstyles[length(idx)]), idx )
   end
   return plt
end
