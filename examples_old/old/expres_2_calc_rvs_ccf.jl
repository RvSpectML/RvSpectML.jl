using RvSpectML
 using Statistics
 using Dates

make_plots = false
include("expres_1_read.jl")
order_list_timeseries = RvSpectML.make_order_list_timeseries(expres_data)
order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
lambda_range_with_data = (min = maximum(d->minimum(d.λ[expres_data[1].metadata[:excalibur_mask]]),expres_data), max = minimum(d->maximum(d.λ[expres_data[1].metadata[:excalibur_mask]]),expres_data) )
RvSpectML.normalize_spectra!(order_list_timeseries,expres_data);

espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks","G2.espresso.mas")
 espresso_df = RvSpectML.read_linelist_espresso(espresso_filename)
 line_list_df = espresso_df |>
     @filter(lambda_range_with_data.min <= _.lambda ) |>
    @filter( _.lambda < lambda_range_with_data.max) |>
#    @filter( _.lambda < 6000.0 ) |>                       # Avoid tellurics at redder wavelengths
#    @filter( _.lambda >6157 || _.lambda < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame

# Setup to run CCF
mask_shape = RvSpectML.CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=1.6)
 line_list = RvSpectML.CCF.BasicLineList(line_list_df.lambda, line_list_df.weight)
 ccf_plan = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=-5e3)
 v_grid = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan)

# Compute CCF's & measure RVs
tstart = now()
 @time ccfs = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
 println("# CCF runtime: ", now()-tstart)
 # Write CCFs to file
 using CSV
 CSV.write("101501_ccfs.csv",Tables.table(ccfs',header=Symbol.(v_grid)))

make_plots = true
# Plot CCFs (normalized by their maximimum)
if make_plots
  using Plots
 plot(v_grid,ccfs./maximum(ccfs,dims=1),label=:none)
 xlabel!("v (m/s)")
 ylabel!("CCF")
end

# Compute RVs
make_plots = true
 rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i],fit_type = "gaussian") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_quad  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i], fit_type = "quadratic") for i in 1:length(order_list_timeseries) ]
 println("RMS of RVs (ESPRESSO lines):  Gaussian: ",std(rvs_ccf_gauss), "   Quadratic: ", std(rvs_ccf_quad))
 # Write to RVs to CSV file
 CSV.write("101501_rvs_ccf_std.csv",DataFrame("Time [MJD]"=>order_list_timeseries.times,"CCF RV [m/s]"=>rvs_ccf_gauss))
 if make_plots
   t_idx = 1
   using Plots
   plot(v_grid,ccfs[:,t_idx],label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
 end


rvs_yale = CSV.read(joinpath(expres_data_path,target_subdir,"101501_activity.csv"), DataFrame)

plt1 = plot()
 rv_ref = mean(rvs_yale[!,"CBC RV [m/s]"])
 scatter!(plt1,rvs_yale[!,"Time [MJD]"],rvs_yale[!,"CBC RV [m/s]"].-rv_ref,yerr=rvs_yale[!,"CBC RV Error [m/s]"], color=1, label="Yale CBC")
 rv_ref = mean(rvs_yale[!,"CCF RV [m/s]"])
 scatter!(plt1,rvs_yale[!,"Time [MJD]"],rvs_yale[!,"CCF RV [m/s]"].-rv_ref,yerr=rvs_yale[!,"CCF RV Error [m/s]"], color=2, label="Yale CCF")
 times_plt = order_list_timeseries.times .- 40000
 rv_ref = mean(rvs_ccf_gauss)
 scatter!(plt1,times_plt,rvs_ccf_gauss.-mean(rvs_ccf_gauss),color=3,label="Gaussian")
 plt2 = plot( legend=:topleft)
 rv_ref = rvs_yale[!,"CCF RV [m/s]"].-mean(rvs_yale[!,"CCF RV [m/s]"]).+mean(rvs_yale[!,"CBC RV [m/s]"])
 scatter!(plt2,rvs_yale[!,"Time [MJD]"],rvs_yale[!,"CBC RV [m/s]"].-rv_ref,yerr=rvs_yale[!,"CBC RV Error [m/s]"], color=1, label="Yale CBC-Yale CCF")
 #rv_ref = rvs_yale[!,"CCF RV [m/s]"].-mean(rvs_yale[!,"CCF RV [m/s]"]).+mean(rvs_yale[!,"CCF RV [m/s]"])
 #scatter!(plt2,rvs_yale[!,"Time [MJD]"],rvs_yale[!,"CCF RV [m/s]"].-rv_ref,yerr=rvs_yale[!,"CCF RV Error [m/s]"], color=2, label="Yale CCF-Yale CCF")
 times_plt = order_list_timeseries.times .- 40000   # Why do I need to subtract this off?
 rv_ref = rvs_yale[!,"CCF RV [m/s]"].-mean(rvs_yale[!,"CCF RV [m/s]"]).+mean(rvs_ccf_gauss)
 scatter!(plt2,times_plt,rvs_ccf_gauss.-rv_ref,color=3,label="RvSpecML CCF - Yale CCF")
 plot(plt1,plt2,layout=(2, 1))

 plt1 = plot()# legend=:topleft)
  plt2 = plot()# legend=:topleft)
  rv_ref1 = mean(rvs_yale[!,"CCF RV [m/s]"])
  rv_ref2 = mean(rvs_ccf_gauss)
  scatter!(plt1,rvs_yale[!,"CCF RV [m/s]"].-rv_ref1,rvs_ccf_gauss.-rv_ref2,xerr=rvs_yale[!,"CCF RV Error [m/s]"],yerr=rvs_yale[!,"CCF RV Error [m/s]"],color=1, label=:none)
  xlabel!(plt1,"RV (m/s) from Yale's CCF")
  ylabel!(plt1,"RV (m/s) from RvSpectML's CCF")
  x_line = range(minimum(rvs_yale[!,"CCF RV [m/s]"]),stop=maximum(rvs_yale[!,"CCF RV [m/s]"]),length=20).-mean(rvs_yale[!,"CCF RV [m/s]"])
  plot!(plt1,x_line,x_line, color=6,label=:none)
  scatter!(plt2,rvs_yale[!,"CCF RV [m/s]"].-rv_ref1,rvs_ccf_gauss.-rv_ref2.-(rvs_yale[!,"CCF RV [m/s]"].-rv_ref1),xerr=rvs_yale[!,"CCF RV Error [m/s]"],yerr=rvs_yale[!,"CCF RV Error [m/s]"],color=1, label=:none)
  xlabel!(plt2,"RV (m/s) from Yale's CCF")
  ylabel!(plt2,"ΔRV (m/s)")
  plot(plt1,plt2,layout=(2, 1))


# Plot Nightly mean RVs
if make_plots
   nbin = 4
   plt_t = (order_list_timeseries.times .- minimum(order_list_timeseries.times) )
   plt = scatter(plt_t,rvs_ccf_gauss.-mean(rvs_ccf_gauss),label="Gaussian")
   scatter!(plt,plt_t,rvs_ccf_quad.-mean(rvs_ccf_quad),label="Quadratic")
   times_binned = RvSpectML.bin_times_nightly(plt_t)
   rvs_ccf_gauss_binned = RvSpectML.bin_rvs_nightly(times=plt_t,rvs=rvs_ccf_gauss)
   rvs_ccf_quad_binned = RvSpectML.bin_rvs_nightly(times=plt_t,rvs=rvs_ccf_quad)
   println("# RMS of nightly binned RVs: ", std(rvs_ccf_gauss_binned), "  ", std(rvs_ccf_quad_binned) ) #, "  ", std(rvs_ccf_cent_binned))
   scatter!(times_binned,rvs_ccf_gauss_binned.-mean(rvs_ccf_gauss),label="Binned Gauss")
   scatter!(times_binned,rvs_ccf_quad_binned.-mean(rvs_ccf_quad),label="Binned Quad")
   xlabel!(plt,"Time (days)")
   ylabel!(plt,"RV (m/s)")
   display(plt)
 end
