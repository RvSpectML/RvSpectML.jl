using RvSpectML
 using Statistics
 using Dates

make_plots = true
include("neid_solar_1_read.jl")
 order_list_timeseries = RvSpectML.make_order_list_timeseries(solar_data)
 order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
 RvSpectML.normalize_spectra!(order_list_timeseries,solar_data);

espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks","G2.espresso.mas")
 espresso_df = RvSpectML.read_mask_espresso(espresso_filename)
 lambda_range_with_data = (min = maximum(d->minimum(d.λ),solar_data), max = minimum(d->maximum(d.λ),solar_data) )
 line_list_df = espresso_df |>
     @filter(lambda_range_with_data.min <= _.lambda ) |>
     @filter( _.lambda < lambda_range_with_data.max) |>
 #    @filter( _.lambda < 6000.0 ) |>                       # Avoid tellurics at redder wavelengths
 #    @filter( _.lambda >6157 || _.lambda < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame

# Setup to run CCF
#mask_entry_doppler_factor = 0.8*calc_doppler_factor(mean(log.(order_list_timeseries.chunk_list[1].data[1].λ[2:end]./
                                                              order_list_timeseries.chunk_list[1].data[1].λ[1:end-1])))
 mask_shape =RvSpectML.CCF.TopHatCCFMask(mask_entry_doppler_factor)
 #line_list = RvSpectML.CCF.BasicLineList(line_list_df.lambda./mask_entry_doppler_factor,
    #                                         line_list_df.lambda.*mask_entry_doppler_factor, line_list_df.weight)
 line_list = RvSpectML.CCF.BasicLineList(line_list_df.lambda_lo, line_list_df.lambda_hi, line_list_df.weight)

 ccf_plan = RvSpectML.CCF.BasicCCFPlan()
 v_grid = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan)


# Compute CCF's & measure RVs
tstart = now()
 @time ccfs = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, line_list, mask_shape=mask_shape, plan=ccf_plan)
 println("# CCF runtime: ", now()-tstart)
 rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i],fit_type = "gaussian") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_quad  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i], fit_type = "quadratic") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_cent  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i], fit_type = "centroid") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_best  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i], fit_type = "bestfit") for i in 1:length(order_list_timeseries) ]
 println("RMS of RVs (ESPRESSO lines):  Gaussian: ",std(rvs_ccf_gauss), "   Quadratic: ", std(rvs_ccf_quad), "   Centroid: ", std(rvs_ccf_cent) )
 if make_plots
   using Plots
  plot(v_grid,ccfs,label=:none)
  xlabel!("v (m/s)")
  ylabel!("CCF")
 end
