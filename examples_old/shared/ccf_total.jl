using RvSpectML
using DataFrames, CSV

function ccf_total(order_list_timeseries::AbstractChunkListTimeseries, line_list_df::DataFrame, pipeline::PipelinePlan; recalc::Bool = false,
                  output_fn_suffix::String = "", range_no_mask_change::Real=30e3, ccf_mid_velocity::Real=0.0, v_step::Real=250.0,
                  mask_scale_factor::Real=1, mask_type::Symbol = :tophat, use_old::Bool = false, use_pixel_vars::Bool = false )
    if need_to(pipeline,:ccf_total) || recalc
      if verbose println("# Computing CCF.")  end
      @assert !need_to(pipeline,:extract_orders)
      @assert !need_to(pipeline,:clean_line_list_tellurics)
      if mask_type == :tophat
          mask_shape = CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      elseif mask_type == :gaussian
          mask_shape = CCF.GaussianCCFMask(order_list_timeseries.inst, σ_scale_factor=mask_scale_factor)
      elseif mask_type == :supergaussian
            mask_shape = CCF.SuperGaussianCCFMask(order_list_timeseries.inst, σ_scale_factor=mask_scale_factor)
      elseif mask_type == :halfcos
          mask_shape = CCF.CosCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      else
        @error("Requested mask shape (" * string(mask_type) * " not avaliable.")
      end

      line_list = CCF.BasicLineList(line_list_df.lambda, line_list_df.weight)
      ccf_plan = CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=range_no_mask_change, step=v_step)
      v_grid = CCF.calc_ccf_v_grid(ccf_plan)
      if use_old
        @time ccfs = CCF.calc_ccf_chunklist_timeseries_old(order_list_timeseries, ccf_plan)
      else
        @time ccfs = CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan, use_pixel_vars=use_pixel_vars)
      end
      if save_data(pipeline, :ccf_total)
         CSV.write(joinpath(output_dir,target_subdir * "_ccfs" * output_fn_suffix * ".csv"),Tables.table(ccfs',header=Symbol.(v_grid)))
         #CSV.write(joinpath(output_dir,target_subdir * "_ccfs_expr.csv"),Tables.table(ccfs_expr',header=Symbol.(v_grid)))
      end
      set_cache!(pipeline, :ccf_total, (ccfs=ccfs, v_grid=v_grid) )
      dont_need_to!(pipeline,:ccf_total)
    end

    if has_cache(pipeline,:ccf_total) return read_cache(pipeline,:ccf_total)
    else   @error("Invalid pipeline state.")          end
end
