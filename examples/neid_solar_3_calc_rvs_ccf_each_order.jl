using RvSpectML
 using Statistics
 using Dates

# Read data
include("neid_solar_1_read.jl")
 order_list_timeseries = RvSpectML.make_order_list_timeseries(solar_data)
 order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
 RvSpectML.normalize_spectra!(order_list_timeseries,solar_data);

# Read & filter the line list file
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
 ccf_plan = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list)
 v_grid = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan)

tstart = now()    # Compute CCFs for each order
 @time order_ccfs = RvSpectML.CCF.calc_order_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan) # line_list, mask_shape=mask_shape, plan=ccf_plan)
 println("# Order CCFs runtime: ", now()-tstart)

  make_plots = true
 if make_plots
    using Plots
  end

if make_plots
   t_idx = 1
   order_idx = 30:40
   plot(v_grid,order_ccfs[:,order_idx,t_idx],label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
 end


rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,order_ccfs[:,j,i],fit_type = "gaussian")  for i in 1:length(order_list_timeseries), j in 1:size(order_ccfs,2) ]
 for i in 1:size(rvs_ccf_gauss,2)
   rvs_ccf_gauss[:,i] .-= mean(rvs_ccf_gauss[:,i])
 end

if make_plots
    using Plots
   nbin = 4
   plt_t = (order_list_timeseries.times .- minimum(order_list_timeseries.times) ) .* 24
   times_binned = RvSpectML.bin_times_consecutive(plt_t, nbin)
   plt = plot()
   plot!(plt,plt_t,rvs_ccf_gauss[:,1],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,10],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,20],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,30],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,30],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,40],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,50],label=:none)
   plot!(plt,plt_t,rvs_ccf_gauss[:,60],label=:none)
   xlabel!("Time (hours)")
   ylabel!("v (m/s)")
end

if make_plots
   order_idx = 5:10:70
   plt = plot()
   for i in order_idx
      rms = std(rvs_ccf_gauss[:,i])
      rvs_binned = RvSpectML.bin_rvs_consecutive(rvs_ccf_gauss[:,i], nbin)
      rms_binned = std(rvs_binned)
      println("  order = ", i, " RMS(order) = ", rms, "  RMS_binned = ", rms_binned)
      scatter!(plt,times_binned, rvs_binned,label=:none)
      plot!(plt,times_binned, rvs_binned,label=:none)
   end
   xlabel!("Time (hours)")
   ylabel!("v (m/s)")
   display(plt)
   end
end
