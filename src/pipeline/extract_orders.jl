
function extract_orders(data::AbstractVector{ST}, pipeline::PipelinePlan; recalc::Bool = false, verbose::Bool = false, orders_to_use = RvSpectMLBase.orders_to_use_default(first(data).inst) ) where { ST<:AbstractSpectra }
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
