# Run code for previous steps with plotting turned off.
make_plots_orig_5 = isdefined(Main,:make_plots) ? make_plots : true
 make_plots = false
 include("neid_solar_3_calc_rvs_proj.jl")
 make_plots = make_plots_orig_5
 if make_plots
   using Plots
 end
using MultivariateStats
 #using Stheno, TemporalGPs

# Set parameters for plotting analysis
plt_order = 42
 plt_order_pix = 3301:3800


espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks","G2.espresso.mas")
espresso_df = RvSpectML.read_mask_espresso(espresso_filename)

lambda_range_with_data = (min = maximum(d->minimum(d.λ),solar_data), max = minimum(d->maximum(d.λ),solar_data) )
line_list_df = espresso_df |>
     @filter(lambda_range_with_data.min <= _.lambda_lo ) |>
    @filter( _.lambda_hi < lambda_range_with_data.max) |>
#    @filter( _.lambda_hi < 6000.0 ) |>
 #   @filter( _.lambda_lo >6157 || _.lambda_hi < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame

mask_entry_doppler_factor = 1+2*mean(log.(order_list_timeseries.chunk_list[1].data[1].λ[2:end]./order_list_timeseries.chunk_list[1].data[1].λ[1:end-1]))
my_mask = hcat(line_list_df.lambda./mask_entry_doppler_factor,line_list_df.lambda.*mask_entry_doppler_factor, line_list_df.weight)
Δv_grid = RvSpectML.CCFTophat.calc_ccf_Δv_grid()

function calc_ccf_chunk(chunk::AbstractChuckOfSpectrum, mask )
  res = RvSpectML.CCFTophat.ccf_1D(chunk.λ,chunk.flux, mask)
end

function calc_ccf_chunklist(chunk_list::AbstractChunkList, mask )
  mapreduce(chunk->calc_ccf_chunk(chunk, mask), +, chunk_list.data)
end

function calc_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries, mask )
  mapreduce(cl->calc_ccf_chunklist(cl, mask),hcat,clt.chunk_list)
end


@time ccfs = calc_ccf_chunklist_timeseries(order_list_timeseries, my_mask)

#maximum(abs2.(ccfs.-ccfs_comp))

if make_plots
  plot(Δv_grid,ccfs,label=:none)
  xlabel("Δv (m/s)")
  ylabel("CCF")
end

rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i],fit_type = "gaussian") for i in 1:length(order_list_timeseries) ]
rvs_ccf_quad  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "quadratic") for i in 1:length(order_list_timeseries) ]
rvs_ccf_cent  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "centroid") for i in 1:length(order_list_timeseries) ]
rvs_ccf_best  = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(Δv_grid,ccfs[:,i], fit_type = "bestfit") for i in 1:length(order_list_timeseries) ]

make_plots = true
if make_plots
  plt_t = (order_list_timeseries.times .- minimum(order_list_timeseries.times) ) .* 24
  plt = scatter(plt_t,rvs_ccf_gauss.-mean(rvs_ccf_gauss),label="Gaussian")
  scatter!(plt,plt_t,rvs_ccf_quad.-mean(rvs_ccf_quad),label="Quadratic")
  #scatter!(plt,plt_t,rvs_ccf_best.-mean(rvs_ccf_best),label="Best fit")
  #scatter!(plt,plt_t,rvs_ccf_cent.mean(rvs_ccf_cent),label="Centroid")
  xlabel!(plt,"Time (hr)")
  ylabel!(plt,"RV (m/s)")
  display(plt)
end
