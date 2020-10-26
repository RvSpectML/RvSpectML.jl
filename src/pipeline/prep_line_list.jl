
""" prepare_line_list( linelist_fn, spectra, pipeline; Δv_to_avoid_tellurics, v_center_to_avoid_tellurics )
linelist_fn is the full path to the line list file
"""
function prepare_line_list( linelist_fn::String, all_spectra::AbstractVector{SpecT}, pipeline::PipelinePlan; recalc::Bool = false,
         orders_to_use = RvSpectML.orders_to_use_default(first(all_spectra).inst),
         Δv_to_avoid_tellurics::Real = EchelleInstruments.default_Δv_to_avoid_tellurics, v_center_to_avoid_tellurics::Real = 0.0, verbose::Bool = false ) where { SpecT <: AbstractSpectra }
   @assert length(linelist_fn) >= 1
   @assert length(all_spectra) >= 1
   local line_list_df
   if need_to(pipeline,:read_line_list) || recalc
      lambda_range_with_good_data = get_λ_range(all_spectra)
      if verbose println("# Reading line list for CCF: ", linelist_fn, ".")  end
      espresso_df = EchelleCCFs.read_linelist_espresso(linelist_fn)
      inst_module = get_inst_module(first(all_spectra).inst)
      #if verbose   println("# extrema(λ) = ", extrema(espresso_df.lambda), "   pre-filter" )   end
      λmin_data = maximum(map(obsid->all_spectra[obsid].λ[1,first(orders_to_use)],1:length(all_spectra)))
      λmax_data = minimum(map(obsid->all_spectra[obsid].λ[end,last(orders_to_use)],1:length(all_spectra)))
      #if verbose   println("# λmin_data = ",λmin_data, "  λmax_data = ",λmax_data)   end
      line_list_df = inst_module.filter_line_list(espresso_df,first(all_spectra).inst, λmin=λmin_data, λmax=λmax_data)
      #if verbose   println("# extrema(λ) = ", extrema(line_list_df.lambda), "   post-filter" )   end
      if eltype(all_spectra) <: EchelleInstruments.AnyEXPRES
         RvSpectMLBase.discard_pixel_mask(all_spectra)
         RvSpectMLBase.discard_excalibur_mask(all_spectra)
      end
      set_cache!(pipeline,:read_line_list,deepcopy(line_list_df))
      dont_need_to!(pipeline,:read_line_list);
    end

    if need_to(pipeline,:clean_line_list_tellurics) || recalc
      if verbose println("# Removing lines with telluric contamination.")  end
      @assert !need_to(pipeline,:read_line_list)
      @assert !need_to(pipeline,:read_spectra)
      line_list_df = make_clean_line_list_from_tellurics(line_list_df, all_spectra,
               Δv_to_avoid_tellurics = Δv_to_avoid_tellurics, v_center_to_avoid_tellurics=v_center_to_avoid_tellurics)
      #if verbose   println("# extrema(λ) = ", extrema(line_list_df.lambda), "   post-clean-tellurics" )   end
      if typeof(first(all_spectra).inst) <: EchelleInstruments.AnyEXPRES
         # RvSpectML.discard_tellurics(all_spectra)  # Keep, since individual line fits use the tellurics info data later
      elseif typeof(first(all_spectra).inst) <: EchelleInstruments.AnyNEID
         # Nothing to do for NEID
      else
         @warn("Removing lines with telluric contamination currently only works with EXPRES data.")
      end

      set_cache!(pipeline,:clean_line_list_tellurics,deepcopy(line_list_df))
      dont_need_to!(pipeline,:clean_line_list_tellurics);
    end

    if need_to(pipeline,:calc_line_list_snr_weights) || recalc
      order_info = get_order_info(all_spectra, orders_to_use=orders_to_use)
      #println("# order_info contains: ", names(order_info))
      #println("# line_list_df contains: ", names(line_list_df))
      line_list_df = RvSpectML.assign_lines_to_orders(line_list_df, order_info)
      if verbose   println("# extrema(λ) = ", extrema(line_list_df.lambda)) end #, "   post-assign-to-orders" )   end
      RvSpectML.calc_snr_weights_for_lines!(line_list_df, all_spectra)
      set_cache!(pipeline,:calc_line_list_snr_weights,line_list_df)
      dont_need_to!(pipeline,:calc_line_list_snr_weights);
    end

  if     has_cache(pipeline,:calc_line_list_snr_weights) return read_cache(pipeline,:calc_line_list_snr_weights)
  elseif has_cache(pipeline,:clean_line_list_tellurics)  return read_cache(pipeline,:clean_line_list_tellurics)
  elseif has_cache(pipeline,:read_line_list)             return read_cache(pipeline,:read_line_list)
  else   @error("Invalid pipeline state.")               end
end
