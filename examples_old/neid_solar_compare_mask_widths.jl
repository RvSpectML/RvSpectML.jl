using Pkg
 Pkg.activate(".")

verbose = true
 if verbose   println("# Loading RvSpecML")    end
 using RvSpectML
include("shared/scripts.jl")
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

include("read_neid_solar_data_20190918.jl")

order_list_timeseries = extract_orders(all_spectra,pipeline_plan)

line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 21e3)

mask_scale_factors = [ 1.0,  2, 4, 6, 8, 10, 12, 14, 16 ]
 rms_rvs = zeros(length(mask_scale_factors))
 rms_binned_rvs = zeros(length(mask_scale_factors))
 println("Starting old Tophat CCFs...")
 for (i,mask_scale_factor) in enumerate(mask_scale_factors)
   println("# mask_scale_factor = ", mask_scale_factor)
   local (ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_scale_factor=mask_scale_factor, mask_type=:tophat,ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
   alg_fit_rv = RVFromCCF.MeasureRvFromCCFGaussian(frac_of_width_to_fit=3)
   local rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv=alg_fit_rv)
   rms_rvs[i] = std(rvs_ccf.-mean(rvs_ccf))
   rms_binned_rvs[i] = std(RvSpectML.bin_rvs_consecutive(rvs_ccf.-mean(rvs_ccf), 4))
 end

 using Plots
   plt1 = plot(mask_scale_factors, rms_rvs, color=1, label="New Tophat")
   plt2 = plot(mask_scale_factors, rms_binned_rvs, color=1, label=:none)
   plot(plt1,plt2,layout=(2,1) )
   xlabel!("Mask width scale parameter")
   ylabel!("RMS (m/s)")

#=
# Commented out since redundant now that old and new match
mask_scale_factors2 = [ 1.0,  2, 4, 6, 8, 10, 12 ]
 rms_rvs2 = zeros(length(mask_scale_factors2))
 rms_binned_rvs2 = zeros(length(mask_scale_factors2))
 println("Starting old Tophat CCFs...")
 for (i,mask_scale_factor) in enumerate(mask_scale_factors2)
   println("# mask_scale_factor = ", mask_scale_factor)
   (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:tophat, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = true)
   rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
   rms_rvs2[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
   rms_binned_rvs2[i] = std(bin_rvs_consecutive(rvs_ccf_expr.-mean(rvs_ccf_expr), 4))
 end
 scatter!(plt1,mask_scale_factors2, rms_rvs2, color=2, label="Old Tophat")
 scatter!(plt2,mask_scale_factors2, rms_binned_rvs2, color=2, label=:none)
 plot(plt1,plt2,layout=(2,1) )
 #xlabel!("Mask width scale parameter")
 #ylabel!("RMS (m/s)")
=#

mask_scale_factors3 = [ 1.0,  2, 4, 6, 8, 10, 12 ]
 rms_rvs3 = zeros(length(mask_scale_factors3))
 rms_binned_rvs3 = zeros(length(mask_scale_factors3))
  println("Starting Half Cos...")
  for (i,mask_scale_factor) in enumerate(mask_scale_factors3)
    println("# mask_scale_factor = ", mask_scale_factor)
    (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:halfcos, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
    rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
    rms_rvs3[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
    rms_binned_rvs3[i] = std(bin_rvs_consecutive(rvs_ccf_expr.-mean(rvs_ccf_expr), 4))
  end
  plot!(plt1,mask_scale_factors3, rms_rvs3, color=3, label="Half Cos")
  plot!(plt2,mask_scale_factors3, rms_binned_rvs3, color=3, label=:none)
  plot(plt1,plt2,layout=(2,1) )

mask_scale_factors4 = [ 1.0,  2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0 ]
 rms_rvs4 = zeros(length(mask_scale_factors4))
   rms_binned_rvs4 = zeros(length(mask_scale_factors4))
   rms_rvs4_quad = zeros(length(mask_scale_factors4))
   rms_rvs4_cent = zeros(length(mask_scale_factors4))
   println("Starting Gaussian...")
    for (i,mask_scale_factor) in enumerate(mask_scale_factors4)
      println("# mask_scale_factor = ", mask_scale_factor)
      (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:gaussian, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
      rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
      rms_rvs4[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
      rms_binned_rvs4[i] = std(bin_rvs_consecutive(rvs_ccf_expr.-mean(rvs_ccf_expr), 4))
      rvs_ccf_quad = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv= RVFromCCF.MeasureRvFromCCFQuadratic())
      #rvs_ccf_cent = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv= RVFromCCF.MeasureRvFromCCFCentroid())
      rms_rvs4_quad[i] = std(rvs_ccf_quad.-mean(rvs_ccf_quad))
      #rms_rvs4_cent[i] = std(rvs_ccf_cent.-mean(rvs_ccf_cent))
    end

scatter!(plt1,mask_scale_factors4, rms_rvs4, color=4, label="Gaussian- Gaussian")
    scatter!(plt2,mask_scale_factors4, rms_binned_rvs4, color=4, label=:none)
    plot!(plt1,mask_scale_factors4, rms_rvs4_quad, color=5, label="Gaussian- Quadratic")
    #scatter!(plt1,mask_scale_factors4, rms_rvs4_cent, color=6, label="Gaussian- Centroid")
    plot(plt1,plt2,layout=(2,1) )

#=
# Commented out since this is slow, so people only run it when they really want to
mask_scale_factors5 = [ 4.0, 6, 8 ]
 rms_rvs5 = zeros(length(mask_scale_factors5))
   rms_binned_rvs5 = zeros(length(mask_scale_factors5))
   println("Starting Super Gaussian...")
    for (i,mask_scale_factor) in enumerate(mask_scale_factors5)
      println("# mask_scale_factor = ", mask_scale_factor)
      (ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:supergaussian, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
      rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
      rms_rvs5[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
      rms_binned_rvs5[i] = std(bin_rvs_consecutive(rvs_ccf_expr.-mean(rvs_ccf_expr), 4))
    end

plot!(plt1,mask_scale_factors5, rms_rvs5, color=5, label="Super-Gaussian")
    plot!(plt2,mask_scale_factors5, rms_binned_rvs5, color=5, label=:none)
    plot(plt1,plt2,layout=(2,1) )
=#

plt1 = plot(mask_scale_factors, rms_rvs, color=1, label="New Tophat")
 plt2 = plot(mask_scale_factors, rms_binned_rvs, color=1, label=:none)
 plot(plt1,plt2,layout=(2,1) )
 xlabel!("Mask width scale parameter")
 ylabel!("RMS (m/s)")
 #scatter!(plt1,mask_scale_factors2, rms_rvs2, color=2, label="Old Tophat")
 #scatter!(plt2,mask_scale_factors2, rms_binned_rvs2, color=2, label=:none)
 plot!(plt1,mask_scale_factors3, rms_rvs3, color=3, label="Half Cos")
 plot!(plt2,mask_scale_factors3, rms_binned_rvs3, color=3, label=:none)
 plot(plt1,plt2,layout=(2,1) )
 scatter!(plt1,mask_scale_factors4, rms_rvs4, color=4, label="Gaussian- Gaussian")
 scatter!(plt2,mask_scale_factors4, rms_binned_rvs4, color=4, label=:none)
 scatter!(plt1,mask_scale_factors4, rms_rvs4_quad, color=5, label="Gaussian- Quadratic")
 #scatter!(plt1,mask_scale_factors4, rms_rvs4_cent, color=6, label="Gaussian- Centroid")
 #plot!(plt1,mask_scale_factors5, rms_rvs5, color=5, label="Super-Gaussian")
 #plot!(plt2,mask_scale_factors5, rms_binned_rvs5, color=5, label=:none)

ylims!(plt1,0.8,1.0)
 ylims!(plt2,0,0.5)
 plot(plt1,plt2,layout=(2,1) )

savefig("neid_mask_width_comp.png")
