using RvSpectML

""" prepare_line_list_pass1( linelist_fn, spectra, pipeline; Δv_to_avoid_tellurics, v_center_to_avoid_tellurics )
"""
function prepare_line_list_pass1( linelist_fn::String, all_spectra::AbstractVector{SpecT}, pipeline::PipelinePlan; recalc::Bool = false,
         Δv_to_avoid_tellurics::Real = RvSpectMLBase.max_bc, v_center_to_avoid_tellurics::Real = 0.0 ) where { SpecT <: AbstractSpectra }
   @assert length(linelist_fn) >= 1
   @assert length(all_spectra) >= 1
   if need_to(pipeline,:read_line_list) || recalc
      lambda_range_with_good_data = get_λ_range(all_spectra)
      if verbose println("# Reading line list for CCF: ", linelist_for_ccf_filename, ".")  end
      espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks",linelist_for_ccf_filename)
      espresso_df = RvSpectML.read_linelist_espresso(espresso_filename)
      #inst_module = RvSpectML.get_inst_module(typeof(first(all_spectra).inst))
      #line_list_df = EXPRES.filter_line_list(espresso_df,first(all_spectra).inst)
      line_list_df = EXPRES.filter_line_list(espresso_df,first(all_spectra).inst)
      if eltype(all_spectra) <: AnyEXPRES
         RvSpectML.discard_pixel_mask(all_spectra)
         RvSpectML.discard_excalibur_mask(all_spectra)
      end
      set_cache!(pipeline,:read_line_list,line_list_df)
      dont_need_to!(pipeline,:read_line_list);
    end

   if need_to(pipeline,:clean_line_list_tellurics) || recalc
      if verbose println("# Removing lines with telluric contamination.")  end
      @assert !need_to(pipeline,:read_line_list)
      @assert !need_to(pipeline,:read_spectra)
      if typeof(first(all_spectra).inst) <: AnyEXPRES
         line_list_no_tellurics_df = make_clean_line_list_from_tellurics_expres(line_list_df, all_spectra, Δv_to_avoid_tellurics = Δv_to_avoid_tellurics, v_center_to_avoid_tellurics=v_center_to_avoid_tellurics)
               # RvSpectML.discard_tellurics(all_spectra)  # Keep, since individual line fits use the tellurics info data later
               set_cache!(pipeline,:clean_line_list_tellurics,line_list_no_tellurics_df)
         else
         @warn("Removing lines with telluric contamination currently only works with EXPRES data.")
      end
      dont_need_to!(pipeline,:clean_line_list_tellurics);
    end

    if has_cache(pipeline,:clean_line_list_tellurics) return read_cache(pipeline,:clean_line_list_tellurics)
    elseif has_cache(pipeline,:read_line_list)        return read_cache(pipeline,:read_line_list)
    else   @error("Invalid pipeline state.")          end
end
