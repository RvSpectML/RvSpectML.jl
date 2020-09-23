using RvSpectML

function extract_orders(data::AbstractVector{ST}, pipeline::PipelinePlan; recalc::Bool = false, verbose::Bool = false, orders_to_use = RvSpectML.orders_to_use_default(first(data).inst) ) where { ST<:AbstractSpectra }
   if need_to(pipeline,:extract_orders) || recalc
      if verbose println("# Extracting order list timeseries from spectra.")  end
      @assert !need_to(pipeline,:read_spectra)
      order_list_timeseries = make_order_list_timeseries(data, orders_to_use=orders_to_use)
      order_list_timeseries = filter_bad_chunks(order_list_timeseries,verbose=true)
      #RvSpectML.normalize_spectra!(order_list_timeseries,data);
      set_cache!(pipeline,:extract_orders,order_list_timeseries)
      dont_need_to!(pipeline,:extract_orders)
    end
    if has_cache(pipeline,:extract_orders) return read_cache(pipeline,:extract_orders)
    else   @error("Invalid pipeline state.")          end
end

#=
default_min_chunk_length = 75
default_min_telluric_model = 1.0

function make_chunk_from_existing_chunk(parent::RvSpectML.AbstractChunkOfSpectrum, loc::NamedTuple{(:pixels,:order),Tuple{UnitRange{Int64},Int64}}) # pixels::UnitRange, order::Integer
   λ = parent.λ.parent
   flux = parent.flux.parent
   var = parent.var.parent
   return RvSpectML.ChunkOfSpectrum(λ,flux,var,loc)
end

function find_chunks_avoiding_tellurics(order_list_timeseries::AbstractChunkListTimeseries, chid::Integer; min_chunk_length::Integer = default_min_chunk_length, min_telluric_model::Real = default_min_telluric_model )
   @assert 1 <= chid <= num_chunks(order_list_timeseries)
   @assert 3 <= min_chunk_length <= 128

   order = first(order_list_timeseries.chunk_list).data[chid].flux.indices[2]
   mask_affected_by_tellurics_at_any_time = mapreduce(obsid->(.!(view(order_list_timeseries.metadata[obsid][:tellurics], order_list_timeseries.chunk_list[obsid][chid].flux.indices[1], order_list_timeseries.chunk_list[obsid][chid].flux.indices[2])
                                          .>= min_telluric_model)), +, 1:length(order_list_timeseries) )
   order_length = length(order_list_timeseries.chunk_list[1].data[chid].flux)
   start_looking_at = 1
   chunk_list = NamedTuple{(:pixels,:order),Tuple{UnitRange{Int64},Int64}}[]
   while start_looking_at < order_length
      chunk_begin = findfirst(m->m==0,view(mask_affected_by_tellurics_at_any_time,start_looking_at:order_length) )
      if isnothing(chunk_begin)    # There's nothing left unaffected by tellurics in this order.
         break
      end
      chunk_begin += start_looking_at-1 # offset since didn't start at beginning
      chunk_end = findfirst(m->m!=0,view(mask_affected_by_tellurics_at_any_time,chunk_begin:order_length) )
      if isnothing(chunk_end)       # Don't need to break up into smaller chunks further, at least not on account of tellurics
         chunk_end = order_length
         if chunk_end-chunk_begin >= min_chunk_length  # Chunk big enough to be useful
            push!(chunk_list, ( pixels= chunk_begin:chunk_end, order=order ) )
         end
         break
      else
         chunk_end += chunk_begin-1    # offset since didn't start at beginning
         start_looking_at = chunk_end  # Start looking for next chunk at first pixel flagged as having a telluric
         chunk_end -= 1                # End chunk without tellurics one pixel earlier
         println("# start_looking_at = ", start_looking_at, "  chunk_begin = ", chunk_begin, "  chunk_end = ", chunk_end)
         println(" Found chunk w/ length=",chunk_end-chunk_begin-1)
         if chunk_end-chunk_begin >= min_chunk_length  # Chunk big enough to be useful
            push!(chunk_list, ( pixels= chunk_begin:chunk_end, order=order ))
         end
      end
   end
   #return chunk_list

   return map(obsid->RvSpectML.ChunkList(map(loc->make_chunk_from_existing_chunk(order_list_timeseries.chunk_list[obsid].data[chid], loc), chunk_list)),1:length(order_list_timeseries) )
end


function merge_chunk_lists_from_multiple_orders(vclt::Vector{Vector{ACLT}}) where { ACLT<:AbstractChunkList}
   @assert length(vclt) >= 1
   num_orders = length(vclt)
   num_times = length(first(vclt))
   for t_idx in 1:num_times
      for ord_idx in 2:num_orders
         append!(vclt[1][t_idx], vclt[ord_idx][t_idx])
      end
   end
   return vclt[1]
end



function find_wavelength_ranges_with_tellurics(order_list_timeseries::AbstractChunkListTimeseries, chid::Integer; min_chunk_length::Integer = default_min_chunk_length, min_telluric_model::Real = default_min_telluric_model )
   @assert 1 <= chid <= num_chunks(order_list_timeseries)
   @assert 3 <= min_chunk_length <= 128


   order = get_order_index(first(order_list_timeseries.chunk_list).data[chid]) # .flux.indices[2]

   mask_affected_by_tellurics_at_any_time = mapreduce(obsid->(.!(view(order_list_timeseries.metadata[obsid][:tellurics], get_pixels_range(order_list_timeseries.chunk_list[obsid][chid]), get_order_index(order_list_timeseries.chunk_list[obsid][chid]) )
                                          .>= min_telluric_model)), +, 1:length(order_list_timeseries) )
   order_length = length(order_list_timeseries.chunk_list[1].data[chid].flux)
   start_looking_at = 1
   chunk_list = NamedTuple{(:pixels,:order),Tuple{UnitRange{Int64},Int64}}[]
   while start_looking_at < order_length
      chunk_begin = findfirst(m->m==0,view(mask_affected_by_tellurics_at_any_time,start_looking_at:order_length) )
      if isnothing(chunk_begin)    # There's nothing left unaffected by tellurics in this order.
         break
      end
      chunk_begin += start_looking_at-1 # offset since didn't start at beginning
      chunk_end = findfirst(m->m!=0,view(mask_affected_by_tellurics_at_any_time,chunk_begin:order_length) )
      if isnothing(chunk_end)       # Don't need to break up into smaller chunks further, at least not on account of tellurics
         chunk_end = order_length
         if chunk_end-chunk_begin >= min_chunk_length  # Chunk big enough to be useful
            push!(chunk_list, ( pixels= chunk_begin:chunk_end, order=order ) )
         end
         break
      else
         chunk_end += chunk_begin-1    # offset since didn't start at beginning
         start_looking_at = chunk_end  # Start looking for next chunk at first pixel flagged as having a telluric
         chunk_end -= 1                # End chunk without tellurics one pixel earlier
         println("# start_looking_at = ", start_looking_at, "  chunk_begin = ", chunk_begin, "  chunk_end = ", chunk_end)
         println(" Found chunk w/ length=",chunk_end-chunk_begin-1)
         if chunk_end-chunk_begin >= min_chunk_length  # Chunk big enough to be useful
            push!(chunk_list, ( pixels= chunk_begin:chunk_end, order=order ))
         end
      end
   end
   #return chunk_list

   return map(obsid->RvSpectML.ChunkList(map(loc->make_chunk_from_existing_chunk(order_list_timeseries.chunk_list[obsid].data[chid], loc), chunk_list)),1:length(order_list_timeseries) )
end
=#
