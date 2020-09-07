using Pkg
 Pkg.activate(".")

verbose = true
if verbose   println("# Loading RvSpecML")    end
 using RvSpectML
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics

# USER:  The default paths that specify where datafiles can be entered here or overridden in examples/data_paths.jl
target_subdir = "101501"   # USER: Replace with directory of your choice
 fits_target_str = "101501"
 output_dir = "examples/output"
 default_paths_to_search = [pwd(),"examples",joinpath(pkgdir(RvSpectML),"examples"),"/gpfs/group/ebf11/default/ebf11/expres/inputs"]

 make_plots = true
 write_ccf_to_csv = false
 write_rvs_to_csv = false
 write_template_to_csv = false
 write_spectral_grid_to_jld2 = false
 write_dcpca_to_csv = false

has_loaded_data = false
 has_computed_ccfs = false
 has_computed_rvs = false
 has_computed_tempalte = false
 has_computed_dcpca = false

if !has_loaded_data
  df_files = make_manifest(target_subdir, EXPRES )
  eval(code_to_include_param_jl())

  if verbose println("# Reading in FITS files.")  end
  @time expres_data = map(EXPRES.read_data,eachrow(df_files_use))
  has_loaded_data = true

  if verbose println("# Extracting order list timeseries from spectra.")  end
    order_list_timeseries = RvSpectML.make_order_list_timeseries(expres_data)
    order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
    RvSpectML.normalize_spectra!(order_list_timeseries,expres_data);
  end

end


if !has_computed_ccfs
 if verbose println("# Reading line list for CCF: ", linelist_for_ccf_filename, ".")  end
 lambda_range_with_good_data = get_λ_range(expres_data)
 espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks",linelist_for_ccf_filename)
 espresso_df = RvSpectML.read_linelist_espresso(espresso_filename)
 line_list_df = EXPRES.filter_line_list(espresso_df,first(expres_data).inst)

    # Compute CCF's & measure RVs
 if verbose println("# Computing CCF.")  end
 mask_shape = RvSpectML.CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=tophap_ccf_mask_scale_factor)
 line_list = RvSpectML.CCF.BasicLineList(line_list_df.lambda, line_list_df.weight)
 ccf_plan = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity)
 v_grid = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan)
 @time ccfs = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
 # Write CCFs to file
 if write_ccf_to_csv
    using CSV
    CSV.write(target_subdir * "_ccfs.csv",Tables.table(ccfs',header=Symbol.(v_grid)))
 has_computed_ccfs = true
 end

if make_plots
   using Plots
   t_idx = 1
   using Plots
   plot(v_grid,ccfs[:,t_idx],label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
 end


if !has_computed_rvs
 if verbose println("# Measuring RVs from CCF.")  end
 #rvs_ccf_gauss = [ RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs[:,i],fit_type = "gaussian") for i in 1:length(order_list_timeseries) ]
 rvs_ccf_gauss = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs,fit_type = "gaussian")
 # Store estimated RVs in metadata
 map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf_gauss[i]-mean(rvs_ccf_gauss), 1:length(order_list_timeseries) )
 if write_rvs_to_csv
    using CSV
    CSV.write(target_subdir * "_rvs_ccf.csv",DataFrame("Time [MJD]"=>order_list_timeseries.times,"CCF RV [m/s]"=>rvs_ccf_gauss))
 end
 has_computed_rvs = true
end

if !has_computed_tempalte
   if verbose println("# Making template spectra.")  end
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries)

   if write_template_to_csv
      using CSV
      CSV.write(target_subdir * "_template.csv",DataFrame("λ"=>spectral_orders_matrix.λ,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
   if write_spectral_grid_to_jld2
      using JLD2, FileIO
      save(target_subdir * "_matrix.jld2", Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
end

if !has_computed_dcpca
   if verbose println("# Performing Doppler constrained PCA analysis.")  end
   using MultivariateStats
   dcpca_out, M = RvSpectML.DCPCA.doppler_constrained_pca(spectral_orders_matrix.flux, deriv, rvs_ccf_gauss)
   frac_var_explained = 1.0.-cumsum(principalvars(M))./tvar(M)
   println("# Fraction of variance explained = ", frac_var_explained[1:min(5,length(frac_var_explained))])
   if write_dcpca_to_csv
      using CSV
      CSV.write(target_subdir * "_dcpca_basis.csv",  Tables.table(M.proj) )
      CSV.write(target_subdir * "_dcpca_scores.csv", Tables.table(dcpca_out) )
   end
end

#=
using MultivariateStats
function doppler_constrained_pca(flux::AbstractArray{T1,2}, deriv::AbstractVector{T2}, rvs::AbstractVector{T3} )  where { T1<:Real, T2<:Real, T3<:Real  }#, times::AbstractVector, rvs::AbstractVector)
  fm_perp = RvSpectML.compute_spectra_perp_doppler_shift(flux,deriv, rvs )
  idx_good_obs = 1:length(rvs)
  M = fit(PCA, fm_perp[:,idx_good_obs]; maxoutdim=10)
  pca_out = MultivariateStats.transform(M,fm_perp[:,idx_good_obs])
  return pca_out, M
end

=#
if make_plots
  using Plots
  # Set parameters for plotting analysis
  plt_order = 42
  plt_order_pix = 4500:5000
  plt = scatter(frac_var_explained, xlabel="Number of PCs", ylabel="Frac Variance Unexplained")
  display(plt)
end

if make_plots
  plt_order = 13
  RvSpectML.plot_basis_vectors(order_grids, f_mean, deriv, M.proj, idx_plt = spectral_orders_matrix.chunk_map[plt_order], num_basis=3 )
  #xlims!(5761.5,5766)
end

if make_plots
  RvSpectML.plot_basis_scores(order_list_timeseries.times, rvs_ccf_gauss, pca_out, num_basis=3 )
end

if make_plots
  RvSpectML.plot_basis_scores_cor( rvs_ccf_gauss, pca_out, num_basis=3)
end

has_loaded_data = false
has_computed_ccfs = false
has_computed_rvs = false
has_computed_tempalte = false
has_computed_dcpca = false
