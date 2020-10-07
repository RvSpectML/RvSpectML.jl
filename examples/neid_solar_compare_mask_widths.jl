using Pkg
 Pkg.activate(".")

verbose = true
 if verbose   println("# Loading RvSpecML")    end
 using RvSpectML
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_neid_solar_data_20190918.jl"))

order_list_timeseries = extract_orders(all_spectra,pipeline_plan)

line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 21e3)

mask_scale_factors = [ 1.0,  2, 4, 6, 8, 10, 12, 14, 16 ]
 rms_rvs = zeros(length(mask_scale_factors))
 rms_binned_rvs = zeros(length(mask_scale_factors))
 println("Starting Tophat CCFs (w/o var)...")
 for (i,mask_scale_factor) in enumerate(mask_scale_factors)
   println("# mask_scale_factor = ", mask_scale_factor)
   local (ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_scale_factor=mask_scale_factor, mask_type=:tophat,ccf_mid_velocity=ccf_mid_velocity, recalc=true, use_old = false)
   alg_fit_rv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=3)
   local rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv=alg_fit_rv)
   rms_rvs[i] = std(rvs_ccf.-mean(rvs_ccf))
   rms_binned_rvs[i] = std(RvSpectML.bin_rvs_consecutive(rvs_ccf.-mean(rvs_ccf), 4))
 end

if make_plot(pipeline_plan,:ccf_total)
   using Plots
   plt1 = plot(mask_scale_factors, rms_rvs, color=1, label="Tophat w/o var")
   plt2 = plot(mask_scale_factors, rms_binned_rvs, color=1, label=:none)
   xlabel!("Mask width scale parameter")
   plot(plt1,plt2,layout=(2,1) )
   ylabel!("RMS (m/s)")
end


mask_scale_factors2 = [ 1.0,  2, 4, 6, 8, 10, 12 ]
 rms_rvs2 = zeros(length(mask_scale_factors2))
 rms_binned_rvs2 = zeros(length(mask_scale_factors2))
 println("Starting Tophat CCFs (w/ vars)...")
 for (i,mask_scale_factor) in enumerate(mask_scale_factors2)
   println("# mask_scale_factor = ", mask_scale_factor)
   ((ccfs_expr, ccf_vars_expr), v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_type=:tophat, mask_scale_factor=mask_scale_factor,  ccf_mid_velocity=ccf_mid_velocity, calc_ccf_var=true, recalc=true)
   rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, ccf_vars_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false)
   rms_rvs2[i] = std(rvs_ccf_expr.-mean(rvs_ccf_expr))
   rms_binned_rvs2[i] = std(bin_rvs_consecutive(rvs_ccf_expr.-mean(rvs_ccf_expr), 4))
 end
 if make_plot(pipeline_plan,:ccf_total)
   scatter!(plt1,mask_scale_factors2, rms_rvs2, color=2, label="Tophat w/ var")
   scatter!(plt2,mask_scale_factors2, rms_binned_rvs2, color=2, label=:none)
   plot(plt1,plt2,layout=(2,1) )
   xlabel!("Mask width scale parameter")
   ylabel!("RMS (m/s)")
 end

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
  if make_plot(pipeline_plan,:ccf_total)
    plot!(plt1,mask_scale_factors3, rms_rvs3, color=3, label="Half Cos")
    scatter!(plt2,mask_scale_factors3, rms_binned_rvs3, color=3, label=:none)
    plot(plt1,plt2,layout=(2,1))
  end

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
      rvs_ccf_quad = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv= EchelleCCFs.MeasureRvFromCCFQuadratic())
      #rvs_ccf_cent = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true, bin_consecutive=4, bin_nightly=false, alg_fit_rv= EchelleCCFs.MeasureRvFromCCFCentroid())
      rms_rvs4_quad[i] = std(rvs_ccf_quad.-mean(rvs_ccf_quad))
      #rms_rvs4_cent[i] = std(rvs_ccf_cent.-mean(rvs_ccf_cent))
    end

  if make_plot(pipeline_plan,:ccf_total)
    scatter!(plt1,mask_scale_factors4, rms_rvs4, color=4, label="Gaussian- Gaussian")
    scatter!(plt2,mask_scale_factors4, rms_binned_rvs4, color=4, label=:none)
    plot!(plt1,mask_scale_factors4, rms_rvs4_quad, color=5, label="Gaussian- Quadratic")
    #scatter!(plt1,mask_scale_factors4, rms_rvs4_cent, color=6, label="Gaussian- Centroid")
    plot(plt1,plt2,layout=(2,1) )
  end

#=
# Commented out since this is the EXPRES LSF model and slow, so people only run it for NEID when they really want to
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
    if make_plot(pipeline_plan,:ccf_total)
      plot!(plt1,mask_scale_factors5, rms_rvs5, color=5, label="Super-Gaussian")
      plot!(plt2,mask_scale_factors5, rms_binned_rvs5, color=5, label=:none)
      plot(plt1,plt2,layout=(2,1) )
    end
=#

if make_plot(pipeline_plan,:ccf_total)
 plt1 = plot(mask_scale_factors, rms_rvs, color=1, label="Tophat w/o vars")
 plt2 = plot(mask_scale_factors, rms_binned_rvs, color=1, label=:none)
 plot(plt1,plt2,layout=(2,1) )
 xlabel!("Mask width scale parameter")
 ylabel!("RMS (m/s)")
 scatter!(plt1,mask_scale_factors2, rms_rvs2, color=2, label="Tophat w/ vars")
 scatter!(plt2,mask_scale_factors2, rms_binned_rvs2, color=2, label=:none)
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
 if save_plot(pipeline_plan,:ccf_total)
   savefig("neid_mask_width_comp.png")
 end
end
