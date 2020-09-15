using Pkg
 Pkg.activate(".")

verbose = true
 if verbose   println("# Loading RvSpecML")    end
 using RvSpectML
 include("shared/scripts.jl")
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

include("read_expres_data_101501.jl")

order_list_timeseries = extract_orders(all_spectra,pipeline_plan)

line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 21e3)

mask_scale_factors = [1.0, 2.0, 4.0, 6.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0 ]
 rms_rvs = zeros(length(mask_scale_factors))
 rms_nightly_rvs = zeros(length(mask_scale_factors))
 rms_within_night_rvs = zeros(length(mask_scale_factors))
 for (i,mask_scale_factor) in enumerate(mask_scale_factors)
   local (ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_expr = false)
   local rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true)
   rms_rvs[i] = std(rvs_ccf.-mean(rvs_ccf))
   rms_nightly_rvs[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf.-mean(rvs_ccf)))
   rms_within_night_rvs[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf.-mean(rvs_ccf))
 end
 using Plots

rms_rvs2 = zeros(length(mask_scale_factors))
 rms_nightly_rvs2 = zeros(length(mask_scale_factors))
 rms_within_night_rvs2 = zeros(length(mask_scale_factors))
 for (i,mask_scale_factor) in enumerate(mask_scale_factors)
   local (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_expr = true)
   local rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true)
   rms_rvs2[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
   rms_nightly_rvs2[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr)))
   rms_within_night_rvs2[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr))
 end

scatter(mask_scale_factors, rms_rvs)
 scatter!(mask_scale_factors, rms_nightly_rvs)
 scatter!(mask_scale_factors, rms_within_night_rvs)
 plot!(mask_scale_factors, rms_rvs2, color=1)
 plot!(mask_scale_factors, rms_nightly_rvs2, color=2)
 plot!(mask_scale_factors, rms_within_night_rvs2, color=3)
 scatter!(mask_scale_factors, rms_rvs2, color=1, label=:none)
 scatter!(mask_scale_factors, rms_nightly_rvs2, color=2, label=:none)
 scatter!(mask_scale_factors, rms_within_night_rvs2, color=3, label=:none)
