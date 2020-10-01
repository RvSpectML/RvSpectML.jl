
function ccf_orders(order_list_timeseries::AbstractChunkListTimeseries,  line_list_df::DataFrame, pipeline::PipelinePlan; verbose::Bool = false, recalc::Bool = false,
   range_no_mask_change::Real=30e3, ccf_mid_velocity::Real=0.0, mask_scale_factor::Real=1, mask_type::Symbol = :tophat )

   if need_to(pipeline, :ccf_orders)  # Compute order CCF's & measure RVs
      if mask_type == :tophat
        mask_shape = TopHatCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      else
        @error("Requested mask shape (" * string(mask_type) * " not avaliable.")
      end

      line_list = BasicLineList(line_list_df.lambda, line_list_df.weight)
      ccf_plan = BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=range_no_mask_change)
      v_grid = calc_ccf_v_grid(ccf_plan)
      tstart = now()    # Compute CCFs for each order
      order_ccfs = EchelleCCFs.calc_order_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
      if verbose   println("# Order CCFs runtime: ", now()-tstart)  end

      set_cache!(pipeline,:ccf_orders, (order_ccfs = order_ccfs, v_grid=v_grid) )
      if save_data(pipeline, :ccf_orders)
         for (i, order) in orders_to_use
            if !(sum(order_ccfs[:,i,:]) > 0)   continue    end
            local t = Tables.table( order_ccfs[:,i,:]', header=Symbol.(v_grid) )
            CSV.write(joinpath(output_dir,target_subdir * "_ccf_order=" * string(order) * ".csv"),t)
         end
      end
      dont_need_to!(pipeline, :ccf_orders);
   end
   if has_cache(pipeline,:ccf_orders) return read_cache(pipeline,:ccf_orders)
   else   @error("Invalid pipeline state.")          end
end
