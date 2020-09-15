using RvSpectML
using CSV

function calc_rvs_from_ccf_total(ccfs::AbstractArray{T1,2}, pipeline::PipelinePlan; v_grid::AbstractVector{T2}, times::AbstractVector{T3},
                                 recalc::Bool = false, verbose::Bool = true) where {T1<:Real, T2<:Real, T3<:Real }
   @assert length(v_grid) == size(ccfs,1)
   need_to!(pipeline_plan, :rvs_ccf_total)
   if need_to(pipeline_plan, :rvs_ccf_total)
      if verbose println("# Measuring RVs from CCF.")  end
      @assert !need_to(pipeline_plan,:ccf_total)

      rvs_ccf = RVFromCCF.measure_rv_from_ccf(v_grid,ccfs,fit_type = :gaussian)
      #rvs_ccf = RVFromCCF.measure_rv_from_ccf(v_grid,ccfs,fit_type = "quadratic")

      rms_rv_nightly = bin_rvs_nightly(times=times,rvs=rvs_ccf.-mean(rvs_ccf))
      rms_rv_within_night = rms_rvs_within_night(times=times,rvs=rvs_ccf.-mean(rvs_ccf))
      if verbose   println("# RMS of RVs: ", std(rvs_ccf), "  nightly RVs: ", std(rms_rv_nightly), "   within night ", (rms_rv_within_night) )   end

      if save_data(pipeline_plan, :rvs_ccf_total)
         CSV.write(joinpath(output_dir,target_subdir * "_rvs_ccf.csv"),DataFrame("Time [MJD]"=>times,"CCF RV [m/s]"=>rvs_ccf))
      end
      set_cache!(pipeline, :rvs_ccf_total, rvs_ccf )
      dont_need_to!(pipeline_plan, :rvs_ccf_total)
   end
   if has_cache(pipeline,:rvs_ccf_total) return read_cache(pipeline,:rvs_ccf_total)
   else   @error("Invalid pipeline state.")          end
end
