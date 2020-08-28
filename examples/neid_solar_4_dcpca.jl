# Run code for previous steps with plotting turned off.
make_plots_orig_4 = isdefined(Main,:make_plots) ? make_plots : true
 make_plots = false
 include("neid_solar_3_calc_rvs_proj.jl")
 make_plots = make_plots_orig_4
 if make_plots
   using Plots
 end
using MultivariateStats

# Set parameters for plotting analysis
plt_order = 42
 plt_order_pix = 3301:3800

fm_perp = RvSpectML.compute_spectra_perp_doppler_shift(spectral_orders_matrix.flux,deriv_orders, ave_good_chunks_rvs )
  idx_good_obs = 1:length(ave_good_chunks_rvs)
  M = fit(PCA, fm_perp[:,idx_good_obs]; maxoutdim=10)
  pca_out = MultivariateStats.transform(M,fm_perp[:,idx_good_obs])
  frac_var_explained = 1.0.-cumsum(principalvars(M))./tvar(M)
  if make_plots
    plt = scatter(frac_var_explained, xlabel="Number of PCs", ylabel="Frac Variance Unexplained")
    display(plt)
  end
  frac_var_explained

if make_plots
  # plt_order = 42
  # plt_order_pix = 3301:3800
  RvSpectML.plot_basis_vectors(order_grids, f_mean_orders, deriv_orders, M.proj, idx_plt = spectral_orders_matrix.chunk_map[plt_order][plt_order_pix] )
end

if make_plots
  RvSpectML.plot_basis_scores(plt_times, ave_good_chunks_rvs, pca_out )
end

if make_plots
  RvSpectML.plot_basis_scores_cor( ave_good_chunks_rvs, pca_out)
end
