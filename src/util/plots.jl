using Plots
export plot_spectrum_chunks

function plot_spectrum_chunks(clt::CLT, chunks::Union{Integer,AbstractUnitRange};
    time_idx::Union{Integer,AbstractUnitRange}=1,
    plt::PT = Plots.plot(legend=:none)) where {CLT<:AbstractChunkListTimeseries, PT<:AbstractPlot}
    xmin = Inf # minimum(order_list_df.lambda_lo[chunk_idx])
    xmax = 0   # maximum(order_list_df.lambda_hi[chunk_idx])
    for c in chunks
        t = 1
        plot!(plt,clt.chunk_list[t].data[c].λ ,clt.chunk_list[t].data[c].flux)
        xmin = min(xmin,minimum(clt.chunk_list[t].data[c].λ))
        xmax = max(xmax,maximum(clt.chunk_list[t].data[c].λ))
    end
    xlabel!(plt,"λ (Å)")
    ylabel!(plt,"Flux)")
    #xlims!(xmin,xmax)
    plt
end
