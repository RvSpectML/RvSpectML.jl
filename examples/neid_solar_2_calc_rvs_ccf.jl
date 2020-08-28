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
     @filter(lambda_range_with_data.min <= _.lambda_lo ) |>
    @filter( _.lambda_hi < lambda_range_with_data.max) |>
#    @filter( _.lambda_hi < 6000.0 ) |>
 #   @filter( _.lambda_lo >6157 || _.lambda_hi < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame
 mask_entry_doppler_factor = 1+0.8*mean(log.(order_list_timeseries.chunk_list[1].data[1].λ[2:end]./order_list_timeseries.chunk_list[1].data[1].λ[1:end-1]))
 my_mask = hcat(line_list_df.lambda./mask_entry_doppler_factor,line_list_df.lambda.*mask_entry_doppler_factor, line_list_df.weight)
 Δv_grid = RvSpectML.CCFTophat.calc_ccf_Δv_grid()


tstart = now()
 @time ccfs = RvSpectML.CCFTophat.calc_ccf_chunklist_timeseries(order_list_timeseries, my_mask)
 println("# CCF runtime: ", now()-tstart)
 rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i],fit_type = "gaussian") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_quad  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "quadratic") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_cent  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "centroid") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_best  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "bestfit") for i in 1:length(order_list_timeseries) ]
 println("RMS of RVs (ESPRESSO lines):  Gaussian: ",std(rvs_ccf_gauss), "   Quadratic: ", std(rvs_ccf_quad), "   Centroid: ", std(rvs_ccf_cent) )
 if make_plots
   using Plots
  plot(Δv_grid,ccfs,label=:none)
  xlabel!("Δv (m/s)")
  ylabel!("CCF")
 end

 #make_plots = true
 if make_plots
   nbin = 4
   plt_t = (order_list_timeseries.times .- minimum(order_list_timeseries.times) ) .* 24
   plt = scatter(plt_t,rvs_ccf_gauss.-mean(rvs_ccf_gauss),label="Gaussian")
   scatter!(plt,plt_t,rvs_ccf_quad.-mean(rvs_ccf_quad),label="Quadratic")
   #scatter!(plt,plt_t,rvs_ccf_best.-mean(rvs_ccf_best),label="Best fit")
   #scatter!(plt,plt_t,rvs_ccf_cent.mean(rvs_ccf_cent),label="Centroid")
   times_binned = RvSpectML.bin_times(plt_t, nbin)
   rvs_ccf_gauss_binned = RvSpectML.bin_times(rvs_ccf_gauss, nbin)
   println("# RMS of binned (n=",nbin,") RVs: ", std(rvs_ccf_gauss_binned))
   scatter!(times_binned,rvs_ccf_gauss_binned.-mean(rvs_ccf_gauss),label="Binned")
   plot!(times_binned,rvs_ccf_gauss_binned.-mean(rvs_ccf_gauss),label=:none,color=3)
   xlabel!(plt,"Time (hr)")
   ylabel!(plt,"RV (m/s)")
   display(plt)
 end

0



# Try using Alex's clean lines
vald_filename = joinpath(ancilary_data_path,"VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas")
 vald_df = RvSpectML.read_mask_vald(vald_filename)
 ssn_out = RvSpectML.searchsortednearest(line_list_df.lambda,vald_df.lambda)
 line_list_clean = line_list_df[ssn_out,:]
 @assert maximum(line_list_clean.lambda .- vald_df.lambda) < 0.05
 clean_mask = hcat(line_list_clean.lambda./mask_entry_doppler_factor,line_list_clean.lambda.*mask_entry_doppler_factor, line_list_clean.weight)

 @time ccfs = RvSpectML.CCFTophat.calc_ccf_chunklist_timeseries(order_list_timeseries, my_mask)
  rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i],fit_type = "gaussian") for i in 1:length(order_list_timeseries) ]
  rvs_ccf_quad  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "quadratic") for i in 1:length(order_list_timeseries) ]
  rvs_ccf_cent  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "centroid") for i in 1:length(order_list_timeseries) ]
  rvs_ccf_best  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "bestfit") for i in 1:length(order_list_timeseries) ]
  println("RMS of RVs (Alex's lines):  Gaussian: ",std(rvs_ccf_gauss), "   Quadratic: ", std(rvs_ccf_quad), "   Centroid: ", std(rvs_ccf_cent) )
  if make_plots
    using Plots
   plot(Δv_grid,ccfs,label=:none)
   xlabel!("Δv (m/s)")
   ylabel!("CCF")
  end



# Try breaking CCFs up by order
@time order_ccfs = RvSpectML.CCFTophat.calc_order_ccf_chunklist_timeseries(order_list_timeseries, my_mask)

if make_plots
  order_idx = 8
 #plot(Δv_grid,mapreduce(i->order_ccfs[:,i],hcat, 1:20:size(order_ccfs,2)),label=:none)
 plot(Δv_grid,order_ccfs[:,order_idx,:],label=:none)
 xlabel!("Δv (m/s)")
 ylabel!("CCF")
end


order_ccfs_time_ave = reshape(mean(order_ccfs,dims=3),size(order_ccfs,1),size(order_ccfs,2))
for i in 1:size(order_ccfs_time_ave,2)
   order_ccfs_time_ave[:,i] .-= maximum(order_ccfs_time_ave[:,i])
end
if make_plots
 plot(Δv_grid,order_ccfs_time_ave,label=:none)
 xlabel!("Δv (m/s)")
 ylabel!("CCF")
end
