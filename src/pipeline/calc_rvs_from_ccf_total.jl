
function calc_rvs_from_ccf_total(ccfs::AbstractArray{T1,2}, pipeline_plan::PipelinePlan; v_grid::AbstractVector{T2}, times::AbstractVector{T3},
                                 alg_fit_rv::AbstractMeasureRvFromCCF = MeasureRvFromCCFGaussian(),
                                 bin_nightly::Bool = true, bin_consecutive::Integer = 0 , recalc::Bool = false, verbose::Bool = true) where {T1<:Real, T2<:Real, T3<:Real }
   @assert length(v_grid) == size(ccfs,1)
   need_to!(pipeline_plan, :rvs_ccf_total)
   if need_to(pipeline_plan, :rvs_ccf_total)
      if verbose println("# Measuring RVs from CCF.")  end
      @assert !need_to(pipeline_plan,:ccf_total)
      #fit_gaussian_to_ccf = RVFromCCF.MeasureRvFromCCFGaussian()
      #fit_quadratic_to_ccf = RVFromCCF.MeasureRvFromCCFQuadratic()
      #rvs_ccf = RVFromCCF.measure_rvs_from_ccf(v_grid,ccfs,alg=fit_gaussian_to_ccf)
      #rvs_ccf = RVFromCCF.measure_rvs_from_ccf(v_grid,ccfs,alg=fit_quadratic_to_ccf)
      (rvs_ccf, σ_rvs_ccf) = measure_rvs_from_ccf(v_grid,ccfs,alg=alg_fit_rv)

      if bin_nightly
         rms_rv_nightly = bin_rvs_nightly(times=times,rvs=rvs_ccf.-mean(rvs_ccf))
         rms_rv_within_night = rms_rvs_within_night(times=times,rvs=rvs_ccf.-mean(rvs_ccf))
         if verbose   println("# RMS of RVs: ", std(rvs_ccf), "  nightly RVs: ", std(rms_rv_nightly), "   within night ", (rms_rv_within_night) )   end
      elseif 1 <= bin_consecutive <= floor(size(ccfs,2)//2)
         rms_rv_binned= bin_rvs_consecutive(rvs_ccf.-mean(rvs_ccf),bin_consecutive)
         if verbose   println("# RMS of RVs: ", std(rvs_ccf), "  binned RVs: ", std(rms_rv_binned) )   end
      else
         if verbose   println("# RMS of RVs: ", std(rvs_ccf))   end
      end

      if save_data(pipeline_plan, :rvs_ccf_total)
         CSV.write(joinpath(output_dir,target_subdir * "_rvs_ccf.csv"),DataFrame("Time [MJD]"=>times,"CCF RV [m/s]"=>rvs_ccf))
      end
      set_cache!(pipeline_plan, :rvs_ccf_total, rvs_ccf )
      dont_need_to!(pipeline_plan, :rvs_ccf_total)
   end
   if has_cache(pipeline_plan,:rvs_ccf_total) return read_cache(pipeline_plan,:rvs_ccf_total)
   else   @error("Invalid pipeline_plan state.")          end
end
