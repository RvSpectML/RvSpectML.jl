# Run code for previous steps with plotting turned off.
make_plots_orig_3 = isdefined(Main,:make_plots) ? make_plots : true
 make_plots = false
 include("neid_solar_2_extract_chunks.jl")
 make_plots = make_plots_orig_3
 if make_plots
   using Plots
 end

# Set parameters for this analysis
oversample_fac_chunks = 2
 oversample_fac_orders = 2
 num_spectra_to_bin = 1   # 1 results in no binning in time
 idx_chunk_min = 20
 idx_chunk_max = 21
 plt_order = 42
 plt_order_pix = 3301:3800

chunk_grids = map(c->RvSpectML.make_grid_for_chunk(chunk_list_timeseries,c,oversample_factor=oversample_fac_chunks), 1:num_chunks(chunk_list_timeseries) )
 @time spectra_all = RvSpectML.pack_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_grids)
#=
RvSpectML.GPs.reset_ncalls()
 @time spectra_all_gp = RvSpectML.pack_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_grids[1:1], alg=:GP)
 plt1 = scatter(spectra_all.λ,spectra_all.flux[:,1],markersize=1.5,label="Obs")
 plot!(plt1,spectra_all_gp.λ[spectra_all_gp.chunk_map[1]],spectra_all_gp.flux[spectra_all_gp.chunk_map[1],1],markersize=1.5,label="GP")
 #plt2 = scatter(spectra_all.λ,spectra_all.flux[:,1]-mean(spectra_all.flux,dims=2),markersize=1.5,label="Obs")
 plt2 = plot(spectra_all_gp.λ[spectra_all_gp.chunk_map[1]],spectra_all_gp.flux[spectra_all_gp.chunk_map[1],1].-mean(spectra_all.flux[:,1],dims=2),markersize=1.5,label="GP")
 plot(plt1,plt2,layout=(2,1) )
=#

f_mean = calc_mean_spectrum(spectra_all.flux,spectra_all.var)
 deriv = calc_mean_deriv(spectra_all.flux,spectra_all.var,spectra_all.λ,spectra_all.chunk_map)
 plt_times = (chunk_list_timeseries.times .-minimum(chunk_list_timeseries.times)).*24
 spectra_matrix = spectra_all

if 2 <= num_spectra_to_bin <=20
  spectra_binned = RvSpectML.bin_consecutive_spectra(spectra_all,num_spectra_to_bin)
  f_mean = calc_mean_spectrum(spectra_binned.flux,spectra_binned.var)
  deriv = calc_mean_deriv(spectra_binned.flux,spectra_binned.var,spectra_binned.λ,spectra_binned.chunk_map)
  times_binned = RvSpectML.bin_times(spectra_binned,chunk_list_timeseries.times,num_spectra_to_bin)
  plt_times = (times_binned.-minimum(chunk_list_timeseries.times)).*24
  spectra_matrix = spectra_binned
end

if make_plots
   #idx_chunk_min = 20
   #idx_chunk_max = 21
   global idx_plt = first(spectra_matrix.chunk_map[idx_chunk_min]):last(spectra_matrix.chunk_map[idx_chunk_max])
   mean_in_plt = mean(f_mean[idx_plt])
   std_in_plt = stdm(f_mean[idx_plt],mean_in_plt)
   plt1 = plot()
   plot!(plt1,spectra_matrix.λ[idx_plt],(spectra_matrix.flux[idx_plt,:].-mean_in_plt)./std_in_plt,linecolor=:black, label=:none)
   plot!(plt1,spectra_matrix.λ[idx_plt],(f_mean[idx_plt].-mean_in_plt)./std_in_plt,linecolor=:red, label="Standardized mean spectrum")
   plot!(plt1,spectra_matrix.λ[idx_plt], deriv[idx_plt]./std(deriv[idx_plt]), label="Standardized deriv", linecolor=:green)
   display(plt1)
end


println("Computing RVs using dflux/dlnlambda from chunks.")
 (rvs_1, σ_rvs_1) = RvSpectML.calc_rvs_from_taylor_expansion(spectra_matrix,mean=f_mean,deriv=deriv)
 rms_rvs_1 = std(rvs_1)
 mean_σ_rvs_1 = mean(σ_rvs_1)
 println("# rms(RVs)  = ", rms_rvs_1, "  σ_RVs = ", mean_σ_rvs_1, "   N_obs = ", length(rvs_1) )
 flush(stdout)
 if make_plots
   plt2 = scatter(plt_times, rvs_1, yerr=σ_rvs_1, label="RV chunks", color=:blue)
 end
 plt2



# Analyze spectra as a set of (large sections of) orders
order_grids = map(c->RvSpectML.make_grid_for_chunk(order_list_timeseries,c,oversample_factor=oversample_fac_orders), 1:num_chunks(order_list_timeseries) )
  spectral_orders_matrix = RvSpectML.pack_chunk_list_timeseries_to_matrix(order_list_timeseries,order_grids)
  if 2 <= num_spectra_to_bin <=20
    spectral_orders_matrix = RvSpectML.bin_consecutive_spectra(spectral_orders_matrix,num_spectra_to_bin)
  end
  f_mean_orders = calc_mean_spectrum(spectral_orders_matrix.flux,spectral_orders_matrix.var)
  deriv_orders = calc_mean_deriv(spectral_orders_matrix.flux,spectral_orders_matrix.var,spectral_orders_matrix.λ,spectral_orders_matrix.chunk_map)

if make_plots
  plt_order = 40
  plt_order_pix = 3501:3800
  idx_plt = spectral_orders_matrix.chunk_map[plt_order][plt_order_pix]
  plt_λ = RvSpectML.get_λs(order_grids, idx_plt)
     mean_in_plt = mean(f_mean_orders[idx_plt])
     std_in_plt = stdm(f_mean_orders[idx_plt],mean_in_plt)
     local plt = plot()
     plot!(plt,plt_λ,(f_mean_orders[idx_plt].-mean_in_plt)./std_in_plt,label=:none,linecolor=:black)
     plot!(plt,plt_λ, deriv_orders[idx_plt]./std(deriv_orders[idx_plt]),label=:none, linecolor=:green)
     scatter!(plt,plt_λ,(spectral_orders_matrix.flux[idx_plt,:].-mean_in_plt)./std_in_plt,markersize=1, label=:none)
     xlabel!(plt,"λ (Å)")
     ylabel!(plt,"Mean & Deriv")
     local plt3 = scatter(plt_λ,spectral_orders_matrix.flux[idx_plt,:].-f_mean_orders[idx_plt],markersize=1, label=:none)
     ylabel!(plt3,"Residuals")
     local pltall = plot(plt,plt3,layout=(2,1))
     display(pltall)
end

println("Computing RVs using dflux/dlnlambda from orders.")
  order_rvs = RvSpectML.calc_chunk_rvs_from_taylor_expansion(spectral_orders_matrix,mean=f_mean_orders,deriv=deriv_orders)
  ave_order_rvs = vec( sum(mapreduce(c ->order_rvs[c].rv./order_rvs[c].σ_rv.^2, hcat, 1:length(order_rvs)),dims=2) ./
                     sum(mapreduce(c ->(1.0./order_rvs[c].σ_rv.^2), hcat, 1:length(order_rvs)),dims=2) )
  rms_order_rvs = map(order->std(order.rv), order_rvs)
  mean_order_σrvs = map(order->mean(order.σ_rv), order_rvs)
  println("# rms(ave_order_rvs)/√N = ", std(ave_order_rvs)/sqrt(length(ave_order_rvs)), "  <RMS RVs> = ", mean(rms_order_rvs)#=/sqrt(length(ave_order_rvs))=#, "  <σ_RVs> = ", mean(mean_order_σrvs) )
  if make_plots
    plt4 = histogram(abs.(ave_order_rvs),bins=40,label="|Ave Order RVs|", alpha=0.75)
    histogram!(plt4,rms_order_rvs,bins=40,label="RMS Order RVs", alpha=0.75)
    histogram!(plt4,mean_order_σrvs,bins=40,label="σ Order RVs", alpha=0.75)
    xlabel!("(m/s)")
    ylabel!("Counts")
  end

chunk_rms_cut_off = quantile(rms_order_rvs,0.9)
  idx_good_chunks = findall(x-> x<= chunk_rms_cut_off, rms_order_rvs)
  @assert length(idx_good_chunks)>=1
  ave_good_chunks_rvs = vec( sum(mapreduce(c ->order_rvs[c].rv./order_rvs[c].σ_rv.^2, hcat, idx_good_chunks),dims=2) ./
                      sum(mapreduce(c ->(1.0./order_rvs[c].σ_rv.^2), hcat, idx_good_chunks),dims=2) )
  sigma_good_chunks_rvs = sqrt.(vec( sum(mapreduce(c ->1.0./order_rvs[c].σ_rv.^2, hcat, idx_good_chunks),dims=2) ./
                      sum(mapreduce(c ->(1.0./order_rvs[c].σ_rv.^4), hcat, idx_good_chunks),dims=2) ))

 rms_rvs_ave_good_orders = std(ave_good_chunks_rvs)
 mean_sigma_good_chunks = mean(sigma_good_chunks_rvs)
 println("# rms(RVs_good_orders) = ", rms_rvs_ave_good_orders,  "  <σ_RV good orders> = ", mean_sigma_good_chunks, "   N_obs = ", length(ave_good_chunks_rvs) )
 #map(o->plot!(plt,chunk_list_timeseries.times,order_rvs[o].rv,yerr=order_rvs[o].σ_rv,markersize=1,label=:none),idx_good_chunks )

if make_plots
   plt = plot()
   scatter!(plt,plt_times, rvs_1, yerr=σ_rvs_1, label="Equal weighted chunks", markersize=3, color=:blue, legend=:topleft)
   scatter!(plt,plt_times, ave_order_rvs,label="Ave Order RVs", markersize=3, color=:green)
   scatter!(plt,plt_times, ave_good_chunks_rvs, yerr=sigma_good_chunks_rvs, label="Ave good chunks", markersize=3, color=:red)
   xlabel!("Time (hours)")
   ylabel!("RV (m/s)")
   display(plt)
end
