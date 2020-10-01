using RvSpectML

function extract_orders(data::AbstractVector{ST}, pipeline::PipelinePlan; recalc::Bool = false, verbose::Bool = false ) where { ST<:AbstractSpectra }
   if need_to(pipeline,:extract_orders) || recalc
      if verbose println("# Extracting order list timeseries from spectra.")  end
      @assert !need_to(pipeline,:read_spectra)
      order_list_timeseries = make_order_list_timeseries(data)
      order_list_timeseries = filter_bad_chunks(order_list_timeseries,verbose=true)
      #RvSpectML.normalize_spectra!(order_list_timeseries,data);
      set_cache!(pipeline,:extract_orders,order_list_timeseries)
      dont_need_to!(pipeline,:extract_orders)
    end
    if has_cache(pipeline,:extract_orders) return read_cache(pipeline,:extract_orders)
    else   @error("Invalid pipeline state.")          end
end
