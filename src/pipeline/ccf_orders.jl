function ccf_orders(order_list_timeseries::AbstractChunkListTimeseries,  line_list_df::DataFrame, pipeline::PipelinePlan; verbose::Bool = false, calc_ccf_var::Bool = false, recalc::Bool = false,
   orders_to_use = RvSpectMLBase.orders_to_use_default(order_list_timeseries.inst),
   range_no_mask_change::Real=RvSpectMLBase.max_bc, ccf_mid_velocity::Real=0.0,
    v_step::Real=250.0, v_max=RvSpectMLBase.max_bc, mask_scale_factor::Real=1,
    mask_type::Symbol = :tophat, Δfwhm::AbstractVector{T} = zeros(0),
    allow_nans::Bool = false ) where { T<:Real }

   if need_to(pipeline, :ccf_orders)  || recalc # Compute order CCF's & measure RVs
      if mask_type == :tophat
         mask_shape = CCFs.TopHatCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      elseif mask_type == :gaussian
         mask_shape = CCFs.GaussianCCFMask(order_list_timeseries.inst, σ_scale_factor=mask_scale_factor)
      elseif mask_type == :supergaussian
         mask_shape = CCFs.SuperGaussianCCFMask(order_list_timeseries.inst, σ_scale_factor=mask_scale_factor)
      elseif mask_type == :halfcos
         mask_shape = CCFs.CosCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      else
         @error("Requested mask shape (" * string(mask_type) * " not avaliable.")
      end

      line_list = hasproperty(line_list_df, :order) ? CCFs.BasicLineList2D(line_list_df.lambda, line_list_df.weight, line_list_df.order) :
                  CCFs.BasicLineList(line_list_df.lambda, line_list_df.weight )
      ccf_plan = BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=range_no_mask_change, step=v_step, max=v_max, allow_nans=allow_nans)
      v_grid = calc_ccf_v_grid(ccf_plan)
      tstart = now()    # Compute CCFs for each order

      if calc_ccf_var
         (order_ccfs, order_ccf_vars ) = EchelleCCFs.calc_order_ccf_and_var_chunklist_timeseries(order_list_timeseries, ccf_plan, Δfwhm=Δfwhm)
         set_cache!(pipeline,:ccf_orders, (order_ccfs = order_ccfs, order_ccf_vars = order_ccf_vars, v_grid=v_grid) )
      else
         order_ccfs = EchelleCCFs.calc_order_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan, Δfwhm=Δfwhm)
         set_cache!(pipeline,:ccf_orders, (order_ccfs = order_ccfs, v_grid=v_grid) )
      end
      if verbose   println("# Order CCFs runtime: ", now()-tstart)  end

      if save_data(pipeline, :ccf_orders)
         for (i, order) in orders_to_use
            if !(sum(order_ccfs[:,i,:]) > 0)   continue    end
            local t = Tables.table( order_ccfs[:,i,:]', header=Symbol.(v_grid) )
            CSV.write(joinpath(output_dir,target_subdir * "_ccf_order=" * string(order) * ".csv"),t)
            if calc_ccf_var
               local t = Tables.table( order_ccf_vars[:,i,:]', header=Symbol.(v_grid) )
               CSV.write(joinpath(output_dir,target_subdir * "_ccf_var_order=" * string(order) * ".csv"),t)
            end
         end
      end
      dont_need_to!(pipeline, :ccf_orders);
   end
   if has_cache(pipeline,:ccf_orders) return read_cache(pipeline,:ccf_orders)
   else   @error("Invalid pipeline state.")          end
end
