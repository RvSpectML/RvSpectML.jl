
""" prepare_line_list( linelist_fn, spectra, pipeline; Δv_to_avoid_tellurics, v_center_to_avoid_tellurics )
"""
function prepare_line_list( linelist_fn::String, all_spectra::AbstractVector{SpecT}, pipeline::PipelinePlan; recalc::Bool = false,
         Δv_to_avoid_tellurics::Real = EchelleInstruments.default_Δv_to_avoid_tellurics, v_center_to_avoid_tellurics::Real = 0.0, verbose::Bool = false ) where { SpecT <: AbstractSpectra }
   @assert length(linelist_fn) >= 1
   @assert length(all_spectra) >= 1
   if need_to(pipeline,:read_line_list) || recalc
      lambda_range_with_good_data = get_λ_range(all_spectra)
      if verbose println("# Reading line list for CCF: ", linelist_fn, ".")  end
      espresso_filename = joinpath(pkgdir(EchelleCCFs),"data","masks",linelist_fn)
      espresso_df = EchelleCCFs.read_linelist_espresso(espresso_filename)
      inst_module = get_inst_module(first(all_spectra).inst)
      line_list_df = inst_module.filter_line_list(espresso_df,first(all_spectra).inst)
      #println(line_list_df)
      if eltype(all_spectra) <: EchelleInstruments.AnyEXPRES
         RvSpectMLBase.discard_pixel_mask(all_spectra)
         RvSpectMLBase.discard_excalibur_mask(all_spectra)
      end
      set_cache!(pipeline,:read_line_list,line_list_df)
      dont_need_to!(pipeline,:read_line_list);
    end

   if need_to(pipeline,:clean_line_list_tellurics) || recalc
      if verbose println("# Removing lines with telluric contamination.")  end
      @assert !need_to(pipeline,:read_line_list)
      @assert !need_to(pipeline,:read_spectra)
      line_list_no_tellurics_df = make_clean_line_list_from_tellurics(line_list_df, all_spectra,
               Δv_to_avoid_tellurics = Δv_to_avoid_tellurics, v_center_to_avoid_tellurics=v_center_to_avoid_tellurics)
      if typeof(first(all_spectra).inst) <: EchelleInstruments.AnyEXPRES
         # RvSpectML.discard_tellurics(all_spectra)  # Keep, since individual line fits use the tellurics info data later
      elseif typeof(first(all_spectra).inst) <: EchelleInstruments.AnyNEID
         # Nothing to do for NEID
      else
         @warn("Removing lines with telluric contamination currently only works with EXPRES data.")
      end
      set_cache!(pipeline,:clean_line_list_tellurics,line_list_no_tellurics_df)
      dont_need_to!(pipeline,:clean_line_list_tellurics);
    end

    if has_cache(pipeline,:clean_line_list_tellurics) return read_cache(pipeline,:clean_line_list_tellurics)
    elseif has_cache(pipeline,:read_line_list)        return read_cache(pipeline,:read_line_list)
    else   @error("Invalid pipeline state.")          end

end
