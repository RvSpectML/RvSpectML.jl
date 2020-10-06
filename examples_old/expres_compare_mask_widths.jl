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

line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc)

#mask_scale_factors = [  1, 2, 4.0, 8.0, 10, 12, 16, 18  ]
mask_scale_factors = [  1, 2, 4.0,8 , 12 ]
 rms_rvs = zeros(length(mask_scale_factors))
 rms_nightly_rvs = zeros(length(mask_scale_factors))
 rms_within_night_rvs = zeros(length(mask_scale_factors))
 println("Starting New Tophat CCFs...")
 for (i,mask_scale_factor) in enumerate(mask_scale_factors)
   println("# mask_scale_factor = ", mask_scale_factor)
   local (ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_scale_factor=mask_scale_factor, mask_type=:tophat,
   ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
   alg_fit_rv = RVFromCCF.MeasureRvFromCCFGaussian(frac_of_width_to_fit=5)
   local rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv=alg_fit_rv)
   rms_rvs[i] = std(rvs_ccf.-mean(rvs_ccf))
   rms_nightly_rvs[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf.-mean(rvs_ccf)))
   rms_within_night_rvs[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf.-mean(rvs_ccf))
 end

#=
mask_scale_factors2 = [  1, 2, 4.0, 8.0, 10, 12, 16, 18, 20, 22, 24, 26, 28, 30, 32 ]
 rms_rvs2 = zeros(length(mask_scale_factors2))
 rms_nightly_rvs2 = zeros(length(mask_scale_factors2))
 rms_within_night_rvs2 = zeros(length(mask_scale_factors2))
 println("Starting Old Tophat CCFs...")
 for (i,mask_scale_factor) in enumerate(mask_scale_factors2)
   println("# mask_scale_factor = ", mask_scale_factor)
   (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:tophat, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = true)
   rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
   rms_rvs2[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
   rms_nightly_rvs2[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr)))
   rms_within_night_rvs2[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr))
 end
=#

mask_scale_factors3 = [  4.0, 8.0, 10, 12, 16, 18, 20, 22, 24 , 26, 28, 30, 32]
 rms_rvs3 = zeros(length(mask_scale_factors3))
 rms_nightly_rvs3 = zeros(length(mask_scale_factors3))
 rms_within_night_rvs3 = zeros(length(mask_scale_factors3))
  println("Starting Half Cos...")
  for (i,mask_scale_factor) in enumerate(mask_scale_factors3)
    println("# mask_scale_factor = ", mask_scale_factor)
    (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:halfcos, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
    rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
    rms_rvs3[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
    rms_nightly_rvs3[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr)))
    rms_within_night_rvs3[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr))
  end

mask_scale_factors4 = [  4.0, 6.0, 8.0, 9, 10, 11, 12, 15, 18 ]
 rms_rvs4 = zeros(length(mask_scale_factors4))
 rms_nightly_rvs4 = zeros(length(mask_scale_factors4))
 rms_within_night_rvs4 = zeros(length(mask_scale_factors4))
   println("Starting Gaussian...")
    for (i,mask_scale_factor) in enumerate(mask_scale_factors4)
      println("# mask_scale_factor = ", mask_scale_factor)
      (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:gaussian, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
      rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
      rms_rvs4[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
      rms_nightly_rvs4[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr)))
      rms_within_night_rvs4[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr))

    end


mask_scale_factors5 = [ 6.0, 8.0, 9.0, 10.0]
 rms_rvs5 = zeros(length(mask_scale_factors5))
  rms_nightly_rvs5 = zeros(length(mask_scale_factors5))
  rms_within_night_rvs5 = zeros(length(mask_scale_factors5))

   println("Starting Super Gaussian...")
    for (i,mask_scale_factor) in enumerate(mask_scale_factors5)
      println("# mask_scale_factor = ", mask_scale_factor)
      (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:supergaussian, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
      rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
      rms_rvs5[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
      rms_nightly_rvs5[i] = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr)))
      rms_within_night_rvs5[i] = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_expr.-mean(rvs_ccf_expr))

    end


using Plots

plot(mask_scale_factors, rms_rvs, color=1, label="Tophat New")
    scatter!(mask_scale_factors, rms_nightly_rvs, color=1, label=:none)
    #plot!(mask_scale_factors2, rms_rvs2, color=2, label="Tophat Old")
    #scatter!(mask_scale_factors2, rms_nightly_rvs2, color=2, label=:none)
plot!(mask_scale_factors3, rms_rvs3, color=3, label="Half Cos")
    scatter!(mask_scale_factors3, rms_nightly_rvs3, color=3, label=:none)
    plot!(mask_scale_factors4, rms_rvs4, color=4, label="Gaussian")
    scatter!(mask_scale_factors4, rms_nightly_rvs4, color=4, label=:none)
plot!(mask_scale_factors5, rms_rvs5, color=5, label="Super-Gaussian")
    scatter!(mask_scale_factors5, rms_nightly_rvs5, color=5, label=:none)

plot(mask_scale_factors, rms_rvs, color=6, label="Tophat New")
        scatter!(mask_scale_factors, rms_nightly_rvs, color=6, label=:none)
