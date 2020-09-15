using Pkg|
 Pkg.activate(".")

verbose = true
if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
 include("shared/scripts.jl")
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

include("read_neid_solar_data_20190918.jl")

if need_to(pipeline,:extract_orders)
   if verbose println("# Extracting order list timeseries from spectra.")  end
   @assert !need_to(pipeline,:read_spectra)
   order_list_timeseries = RvSpectML.make_order_list_timeseries(all_spectra)
   order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
   #RvSpectML.normalize_spectra!(order_list_timeseries,all_spectra);
   dont_need_to!(pipeline,:extract_orders);
 end
 order_list_timeseries

if need_to(pipeline,:read_line_list)
   if verbose println("# Reading line list for CCF: ", linelist_for_ccf_filename, ".")  end
   lambda_range_with_good_data = get_λ_range(all_spectra)
   espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks",linelist_for_ccf_filename)
   espresso_df = RvSpectML.read_linelist_espresso(espresso_filename)
   line_list_df = NEID.filter_line_list(espresso_df,first(all_spectra).inst)
   dont_need_to!(pipeline,:read_line_list);
 end
 #line_list_df

line_list_no_tellurics_df = line_list_df
dont_need_to!(pipeline,:clean_line_list_tellurics)
# TODO: Filter line list for tellurics when using non-EXPRES data
#=
#need_to!(pipeline,:clean_line_list_tellurics)
if need_to(pipeline,:clean_line_list_tellurics)
   if verbose println("# Removing lines with telluric contamination.")  end    # Currently only works for EXPRES data
   @assert !need_to(pipeline,:read_line_list)
   @assert !need_to(pipeline,:read_spectra)
   line_list_no_tellurics_df  = make_clean_line_list_from_tellurics_expres(line_list_df, all_spectra, Δv_to_avoid_tellurics = 30.0e3) #14000.0)
   dont_need_to!(pipeline,:clean_line_list_tellurics);
 end
 #line_list_no_tellurics_df
=#

#need_to!(pipeline,:ccf_total)
if need_to(pipeline,:ccf_total)
   if verbose println("# Computing CCF.")  end
   @assert !need_to(pipeline,:extract_orders)
   @assert !need_to(pipeline,:clean_line_list_tellurics)
   mask_shape = RvSpectML.CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=tophap_ccf_mask_scale_factor)
   #line_list = RvSpectML.CCF.BasicLineList(line_list_df.lambda, line_list_df.weight)
   line_list = RvSpectML.CCF.BasicLineList(line_list_no_tellurics_df.lambda, line_list_no_tellurics_df.weight)
   ccf_plan = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=30e3)
   v_grid = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan)
   @time ccfs = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
   mask_shape_expr = RvSpectML.CCF.GaussianCCFMask(order_list_timeseries.inst, scale_factor=9)
   # Warning:  CCF with a shape other than a tophat is still experimental
   ccf_plan_expr = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape_expr, line_list=line_list, midpoint=ccf_mid_velocity)
   @time ccfs_expr = RvSpectML.CCF.calc_ccf_chunklist_timeseries_expr(order_list_timeseries, ccf_plan_expr) #, verbose=true)
   println("# Ratio of max(ccfs_expr)/max(ccfs) = ", mean(maximum(ccfs_expr,dims=1)./maximum(ccfs,dims=1)) )
   #=
   mask_shape_expr2 = RvSpectML.CCF.CosCCFMask(order_list_timeseries.inst, scale_factor=18)   # Why does such a large value help so much?
   ccf_plan_expr2 = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape_expr2, line_list=line_list, midpoint=ccf_mid_velocity)
   @time ccfs_expr2 = RvSpectML.CCF.calc_ccf_chunklist_timeseries_expr(order_list_timeseries, ccf_plan_expr2) #, verbose=true)
   =#
   if save_data(pipeline, :ccf_total)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_ccfs.csv"),Tables.table(ccfs',header=Symbol.(v_grid)))
      CSV.write(joinpath(output_dir,target_subdir * "_ccfs_expr.csv"),Tables.table(ccfs_expr',header=Symbol.(v_grid)))
   end
   dont_need_to!(pipeline,:ccf_total)
 end


if make_plot(pipeline, :ccf_total)
   using Plots
   t_idx = 20
   plt = plot(v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   scatter!(plt,v_grid,ccfs_expr[:,t_idx]./maximum(ccfs_expr[:,t_idx],dims=1),markersize=1.2,label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
   if save_plot(pipeline,:ccf_total)   savefig(plt,joinpath(output_dir, target_subdir * "_ccf_sum.png"))   end
   display(plt)
end


if make_plot(pipeline, :ccf_total)
   using Plots
   zvals = ccfs./maximum(ccfs,dims=1).-mean(ccfs./maximum(ccfs,dims=1),dims=2)
   colorscale = cgrad(:balance)
   plt = heatmap(v_grid,collect(1:size(ccfs,2)),zvals', c=colorscale, clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
   xlabel!("v (m/s)")
   ylabel!("Observation #")
   title!("CCF(v,t)-<CCF>(v) vs time")
   if save_plot(pipeline,:ccf_total)   savefig(plt,joinpath(output_dir, target_subdir * "_ccf_sum_vs_time_heatmap.png"))   end
   display(plt)
end

#need_to!(pipeline, :rvs_ccf_total)
if need_to(pipeline, :rvs_ccf_total)
   if verbose println("# Measuring RVs from CCF.")  end
   @assert !need_to(pipeline,:ccf_total)
   rvs_ccf_gauss = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs,fit_type = :gaussian)
   rvs_ccf_gauss2 = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs_expr,fit_type = :gaussian)
   #rvs_ccf_gauss3 = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid,ccfs_expr2,fit_type = :gaussian)
   println("# RMS of RVs:          Tophat ", std(rvs_ccf_gauss), "   Gaussian ", std(rvs_ccf_gauss))
   rms_rv_binned = bin_rvs_consecutive(rvs_ccf_gauss,4)
   rms_rv_expr_binned = bin_rvs_consecutive(rvs_ccf_gauss2,4)
   println("# RMS of binned RVs:  Tophat ", std(rms_rv_binned), "   Gaussian ", std(rms_rv_expr_binned))
   #rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_gauss)
   #rms_rv_expr_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_gauss2)
   #println("# RMS within night RVs: Tophat ",rms_rv_within_night,  "   Gaussian ",rms_rv_expr_within_night)
   # Store estimated RVs in metadata
   map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf_gauss[i]-mean(rvs_ccf_gauss), 1:length(order_list_timeseries) )
   if save_data(pipeline, :rvs_ccf_total)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_rvs_ccf.csv"),DataFrame("Time [MJD]"=>order_list_timeseries.times,"CCF RV [m/s]"=>rvs_ccf_gauss))
   end
   dont_need_to!(pipeline, :rvs_ccf_total)
end


if make_plot(pipeline, :rvs_ccf_total)
   using Plots
   rvs_ccf_gauss .-= mean(rvs_ccf_gauss)
   rvs_ccf_gauss2 .-= mean(rvs_ccf_gauss2)
   plt = scatter(rvs_ccf_gauss,markersize=3,label="RVs CCF Tophat", legend=:bottomright)
   scatter!(plt,rvs_ccf_gauss2,markersize=3,label="RVs CCF Gaussian")
   ylabel!("v (m/s)")
   xlabel!("Time (#)")
   if save_plot(pipeline,:rvs_ccf_total)   savefig(plt,joinpath(output_dir, target_subdir * "_rvs_ccf_sum.png"))   end
   display(plt)

   #=
   using CSV
   df_yale_resutls = CSV.read(joinpath(homedir(),"Data/EXPRES/inputs/101501/101501_activity.csv"),DataFrame)
   rvs_yale_ccf = df_yale_resutls["CCF RV [m/s]"]
   rvs_yale_ccf .-= mean(rvs_yale_ccf)
   rvs_yale_cbc = df_yale_resutls["CBC RV [m/s]"]
   rvs_yale_cbc .-= mean(rvs_yale_cbc)
   scatter!(rvs_yale_cbc,markersize=2,label="Yale CBC")
   scatter!(rvs_yale_ccf,markersize=2,label="Yale CCF")
   =#
   #diff = rvs_yale_ccf.-rvs_yale_cbc
   #diff = rvs_ccf_gauss.-rvs_yale_cbc
   #diff = rvs_ccf_gauss2.-rvs_yale_cbc
   #diff = rvs_ccf_gauss.-rvs_yale_ccf
   #diff = rvs_ccf_gauss2.-rvs_yale_ccf
end

if make_plot(pipeline, :rvs_ccf_total)
   diff = rvs_ccf_gauss.-rvs_ccf_gauss2
   println(std(diff))
   plt = scatter(order_list_timeseries.times,diff,markersize=4,label="Delta RV")
   ylabel!("Δv (m/s) (Two mask shapes)")
   xlabel!("Time (d)")
   if save_plot(pipeline,:rvs_ccf_total)   savefig(plt,joinpath(output_dir, target_subdir * "_rvs_ccf_sum.png"))   end
   display(plt)
end

#need_to!(pipeline, :ccf_orders)
 if need_to(pipeline, :ccf_orders)  # Compute order CCF's & measure RVs
   tstart = now()    # Compute CCFs for each order
   order_ccfs = RvSpectML.CCF.calc_order_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
   println("# Order CCFs runtime: ", now()-tstart)

   if save_data(pipeline, :ccf_orders)
      using CSV
      inst = EXPRES.EXPRES2D()
      for (i, order) in orders_to_use
         if !(sum(order_ccfs[:,i,:]) > 0)   continue    end
         local t = Tables.table( order_ccfs[:,i,:]', header=Symbol.(v_grid) )
         CSV.write(joinpath(output_dir,target_subdir * "_ccf_order=" * string(order) * ".csv"),t)
      end
   end
   dont_need_to!(pipeline, :ccf_orders);
end

if make_plot(pipeline, :ccf_orders)
   # order ccfs averaged over time
   obs = 1:size(order_ccfs,3)
   #plot(v_grid,reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1), size(order_ccfs,1),size(order_ccfs[:,ord,:],2)), label=:none)
   # resudiuals of order ccfs to time and chunk averaged ccf
   #plt = plot(v_grid,reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,ord,:],2)).-
   #            reshape(sum(order_ccfs[:,ord,:],dims=(2,3))./maximum(sum(order_ccfs[:,ord,:],dims=(2,3))),size(order_ccfs,1) ), label=:none)
   order_labels = map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2))
   zvals = reshape(sum(order_ccfs[:,:,obs],dims=3)./maximum(sum(order_ccfs[:,:,obs],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,:,obs],2)).-
   reshape(sum(order_ccfs[:,:,obs],dims=(2,3))./maximum(sum(order_ccfs[:,:,obs],dims=(2,3))),size(order_ccfs,1) )
   plt = heatmap(v_grid,order_labels, zvals',c = cgrad(:balance), clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )

   xlabel!("v (m/s)")
   ylabel!("Order ID")
   title!("CCF-<CCF> for obs ID=" * string(obs))
   if save_plot(pipeline,:ccf_orders)   savefig(plt,joinpath(output_dir, target_subdir * "_ccf_orders.png"))   end
   display(plt)
end

if make_plot(pipeline, :ccf_orders)
   ord = 3 # TODO:  I think this order illustrates the issue Alex is going to fix with lines dropping off at the boundary of orders.  Check his fix solves the issue.
   #plot(v_grid,reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1), size(order_ccfs,1),size(order_ccfs[:,ord,:],2)), label=:none)
   # resudiuals of order ccfs to time and chunk averaged ccf
   #plt = plot(v_grid,reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,ord,:],2)).-
   #            reshape(sum(order_ccfs[:,ord,:],dims=(2,3))./maximum(sum(order_ccfs[:,ord,:],dims=(2,3))),size(order_ccfs,1) ), label=:none)
   zvals = reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,ord,:],2)).-
   reshape(sum(order_ccfs[:,ord,:],dims=(2,3))./maximum(sum(order_ccfs[:,ord,:],dims=(2,3))),size(order_ccfs,1) )
   plt = heatmap(v_grid,1:size(order_ccfs,3), zvals',c = cgrad(:balance), clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
   xlabel!("v (m/s)")
   ylabel!("Observation ID")
   title!("CCF-<CCF> for order=" * string(ord))
   if save_plot(pipeline,:ccf_orders)   savefig(plt,joinpath(output_dir, target_subdir * "_ccf_orders.png"))   end
   display(plt)
end

RvSpectML.Pipeline.make_plot!(pipeline,:scalpels)
if make_plot(pipeline, :scalpels)
   @assert !need_to(pipeline, :rvs_ccf_total)
   @assert !need_to(pipeline, :ccf_total)
   plt = Scalpels.make_plots_scalpels(rvs_ccf_gauss, ccfs, max_num_basis=4, v_grid=v_grid, times=order_list_timeseries.times, output_path="examples/output/figures")
   display(plt)
end


#need_to!(pipeline, :template)
if need_to(pipeline, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectra.")  end
   @assert !need_to(pipeline,:extract_orders)
   @assert !need_to(pipeline,:rvs_ccf_total)
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries)
   if save_data(pipeline, :template)
      #using CSV
      #CSV.write(joinpath(output_dir,target_subdir * "_template.csv"),DataFrame("λ"=>spectral_orders_matrix.λ,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
      using JLD2, FileIO
      save(joinpath(output_dir,target_subdir * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
   dont_need_to!(pipeline, :template);
end

if make_plot(pipeline,:template)
   chunkid = 10
   idx = spectral_orders_matrix.chunk_map[chunkid]
   plt = scatter(spectral_orders_matrix.λ[idx],(f_mean[idx].-1.0)./maximum(abs.((f_mean[idx].-1.0))),markersize=1.0,label=:none)
   scatter!(plt,spectral_orders_matrix.λ[idx],deriv[idx]./maximum(abs.(deriv[idx])),markersize=1.0,label=:none)
   scatter!(plt,spectral_orders_matrix.λ[idx],deriv2[idx]./maximum(abs.(deriv2[idx])),markersize=1.0,label=:none)
   xlabel!("λ(Å)")
   ylabel!("f(λ), f'(λ), f''(λ), all standardized")
   title!("Template spectrum for chunk " * string(chunkid) )
end

#need_to!(pipeline, :dcpca)
if need_to(pipeline, :dcpca)
   if verbose println("# Performing Doppler constrained PCA analysis.")  end
   @assert !need_to(pipeline,:rvs_ccf_total)
   @assert !need_to(pipeline,:template)
   using MultivariateStats
   Δrvs_dcpca = calc_rvs_from_taylor_expansion(spectral_orders_matrix; mean=f_mean, deriv=deriv).rv
   dcpca_out, M = RvSpectML.DCPCA.doppler_constrained_pca(spectral_orders_matrix.flux, deriv, Δrvs_dcpca)
   frac_var_unexplained = 1.0.-cumsum(principalvars(M))./tvar(M)
   num_basis = length(frac_var_unexplained)
   println("# Fraction of variance unexplained = ", frac_var_unexplained[1:min(5,num_basis)])
   if save_data(pipeline, :dcpca)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_basis.csv"),  Tables.table(M.proj) )
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_scores.csv"), Tables.table(dcpca_out) )
   end
   dont_need_to!(pipeline, :dcpca);
end


if make_plot(pipeline, :dcpca)  # Ploting results from DCPCA
  using Plots
  # Set parameters for plotting analysis
  plt = scatter(frac_var_unexplained, xlabel="Number of PCs", ylabel="Frac Variance Unexplained")
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca_frac_var.png"))   end
  display(plt)
end

if make_plot(pipeline, :dcpca)
  plt_order = 10
  plt = RvSpectML.plot_basis_vectors(order_grids, f_mean, deriv, M.proj, idx_plt = spectral_orders_matrix.chunk_map[plt_order], num_basis=2 )
  #xlims!(5761.5,5766)
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca_basis.png"))   end
  display(plt)
end

if make_plot(pipeline, :dcpca)
  plt = RvSpectML.plot_basis_scores(order_list_timeseries.times, rvs_ccf_gauss, dcpca_out, num_basis=min(5,num_basis) )
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca_scores.png"))   end
  display(plt)
end

if make_plot(pipeline, :dcpca)
  plt = RvSpectML.plot_basis_scores_cor( rvs_ccf_gauss, dcpca_out, num_basis=min(5,num_basis))
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca_scores_cor.png"))   end
  display(plt)
end

#need_to!(pipeline,:fit_lines)
if need_to(pipeline,:fit_lines)
   if verbose println("# Performing fresh search for lines in template spectra.")  end
   cl = ChunkList(map(grid->ChuckOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map))
   #= # We're done with the spectral_orders_matrix, so we can release the memory now
   spectral_orders_matrix = nothing
   GC.gc()
   need_to!(pipeline,:template)
   =#
   lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(min_deriv2=3))  # TODO: Automate threshold for finding a line

   if verbose println("# Finding above lines in all spectra.")  end
   @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_in_template )

   if save_data(pipeline,:fit_lines)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_linefinder_lines.csv"), lines_in_template )
      CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits.csv"), fits_to_lines )
      #CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits_clean.csv"), lines_to_try )
   end
   dont_need_to!(pipeline,:fit_lines);
end


#need_to!(pipeline,:clean_line_list_blends)
if need_to(pipeline,:clean_line_list_blends)
   # Exploratory data analysis of the distribution of line properties over time to figure out how to select "clean" lines
   @assert !need_to(pipeline,:fit_lines)
   fit_distrib = fits_to_lines |> @groupby(_.line_id) |>
            @map( { median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
                   std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc), line_id=_.line_id, frac_converged=mean(_.fit_converged)  } ) |>
            @filter(_.frac_converged == 1.0 ) |> DataFrame
                  # , min_telluric_model_all_obs=minimum(_.min_telluric_model_this_obs),
   good_lines = fit_distrib |>
            #@filter(_.min_telluric_model_all_obs == 1.0 ) |> # Already done
            @filter( _.median_depth > 0.1 ) |>
            @filter( _.median_σ² < 0.015) |>
            @filter( _.std_σ² < 0.004) |>
            @filter(_.std_b < 0.06) |>
            @filter( _.std_depth < 0.006 ) |>
            #@filter(_.std_a < 0.025) |>             # Unncessary?
            #@filter( -0.25 < _.median_b < 0.25) |>  # Unncessary?
            DataFrame
      println("# Found ", size(lines_in_template,1), " lines, including ",size(good_lines,1), " good lines, rejected ", size(lines_in_template,1)-size(good_lines,1), " lines.")
      bad_lines = fit_distrib |> @filter(_.std_λc >= 0.002 ) |> DataFrame
      good_lines_high_scatter = good_lines |> @filter(_.std_λc >= 0.002 ) |> DataFrame
      if size(bad_lines,1) >= 1
         println("# ", size(bad_lines,1), " lines have large λc scatter, including ", size(good_lines_high_scatter,1), " good lines.")
      end
      lines_to_try = lines_in_template[first.(good_lines[!,:line_id]),:]
      dont_need_to!(pipeline,:clean_line_list_blends);
 end
 fit_distrib

#=if make_plot(pipeline,:fit_lines) || true # Need to add key to pipeline
scatter(good_lines.median_depth,good_lines.std_λc,markersize=2,label=:none)
scatter(good_lines.std_depth,good_lines.std_λc,markersize=2,label=:none)
scatter(good_lines.std_σ²,good_lines.std_λc,markersize=2,label=:none)
scatter(good_lines.median_σ²,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.median_b,good_lines.std_λc,markersize=2,label=:none)
scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.median_b,good_lines.std_λc,markersize=2,label=:none)
scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)
end
=#

# Reset steps of pipeline to rerun with new linelist.
need_to!(pipeline, :ccf_total);
 need_to!(pipeline, :rvs_ccf_total);
 need_to!(pipeline, :ccf_orders);
 need_to!(pipeline, :rvs_ccf_orders);
 need_to!(pipeline, :template);
 need_to!(pipeline, :dcpca);

#need_to!(pipeline,:ccf_total)
if need_to(pipeline,:ccf_total)
   if verbose println("# Computing CCFs with new line list.")  end
   @assert !need_to(pipeline,:extract_orders)
   @assert !need_to(pipeline,:clean_line_list_blends)
   #mask_shape = RvSpectML.CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=tophap_ccf_mask_scale_factor)
   perm = sortperm(lines_to_try.fit_λc)
   line_list2 = RvSpectML.CCF.BasicLineList(lines_to_try.fit_λc[perm], lines_to_try.fit_depth[perm] )
   #mask_shape = RvSpectML.CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=tophap_ccf_mask_scale_factor*3)
   mask_shape = RvSpectML.CCF.GaussianCCFMask(order_list_timeseries.inst, scale_factor=4)
   ccf_plan2 = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape, step=100.0, line_list=line_list2, midpoint=0.0)
   v_grid2 = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan2)
   @time ccfs2 = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan2)
   # Write CCFs to file
   if save_data(pipeline,:ccf_total)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_ccfs2.csv"),Tables.table(ccfs2',header=Symbol.(v_grid)))
   end
   dont_need_to!(pipeline,:ccf_total)
end

if make_plot(pipeline, :ccf_total)
   using Plots
   t_idx = 20
   using Plots
   plt = plot(v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   scatter!(plt,v_grid2,ccfs2[:,t_idx]./maximum(ccfs2[:,t_idx],dims=1),markersize=1.2,label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
   if save_plot(pipeline,:ccf_total)   savefig(plt,joinpath(output_dir, target_subdir * "_ccf2_sum.png"))   end
   display(plt)
end

#need_to!(pipeline,:rvs_ccf_total)
if need_to(pipeline, :rvs_ccf_total)
   if verbose println("# Measuring RVs from CCF.")  end
   @assert !need_to(pipeline,:ccf_total)
   rvs_ccf_gauss2 = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid2,ccfs2,fit_type = :gaussian)
   # Store estimated RVs in metadata
   println("# RMS of RVs:           Orig lines ", std(rvs_ccf_gauss),  "    Cleaned lines ", std(rvs_ccf_gauss2))
   rms_rv_binned = bin_rvs_consecutive(rvs_ccf_gauss,4)
   rms_rv_binned2 = bin_rvs_consecutive(rvs_ccf_gauss2,4)
   println("# RMS of binned RVs:  Orig lines ", std(rms_rv_binned2), "   Cleanded lines ", std(rms_rv_expr_binned2))
   #rms_rv_nightly2 = bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_gauss2)
   #println("# RMS of nightly RVs:   Orig lines ", std(rms_rv_nightly), "    Cleaned lines ", std(rms_rv_nightly2))
   #rms_rv_within_night2 = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_gauss2)
   #println("# RMS within night RVs: Orig lines ",rms_rv_within_night,  "   Cleaned liens ",rms_rv_within_night2)
   map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf_gauss2[i]-mean(rvs_ccf_gauss2), 1:length(order_list_timeseries) )
   if save_data(pipeline, :rvs_ccf_total)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_rvs_ccf2.csv"),DataFrame("Time [MJD]"=>order_list_timeseries.times,"CCF RV [m/s]"=>rvs_ccf_gauss))
   end
   dont_need_to!(pipeline, :rvs_ccf_total);
end


#need_to!(pipeline, :template)
if need_to(pipeline, :template)
   @assert !need_to(pipeline,:rvs_ccf_total)
   @assert !need_to(pipeline,:clean_line_list_blends)
   # Revert to velocities before cleaning for now
   map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf_gauss[i]-mean(rvs_ccf_gauss), 1:length(order_list_timeseries) )
   #chunk_list_df2 = lines_to_try |> @select(:fit_min_λ,:fit_max_λ) |> @rename(:fit_min_λ=>:lambda_lo, :fit_max_λ=>:lambda_hi) |> DataFrame
   expand_chunk_factor = 4
   chunk_list_df2 = lines_to_try |> @select(:fit_λc,:fit_min_λ,:fit_max_λ) |> @map({:lambda_lo=>_.fit_λc-expand_chunk_factor*(_.fit_λc-_.fit_min_λ), :lambda_hi=>_.fit_λc+expand_chunk_factor*(_.fit_max_λ-_.fit_λc)}) |> DataFrame
   chunk_list_timeseries2 = RvSpectML.make_chunk_list_timeseries(all_spectra,chunk_list_df2)

   # Check that no NaN's included
   #(chunk_list_timeseries2, chunk_list_df2) = RvSpectML.filter_bad_chunks(chunk_list_timeseries2,chunk_list_df2)
   chunk_list_timeseries2 = RvSpectML.filter_bad_chunks(chunk_list_timeseries2)
   println(size(chunk_list_df2), " vs ", num_chunks(chunk_list_timeseries2) )

   if verbose println("# Making template spectra.")  end
   @time ( spectral_orders_matrix2, f_mean2, var_mean2, deriv_2, deriv2_2, order_grids2 )  = RvSpectML.make_template_spectra(chunk_list_timeseries2)

   if save_data(pipeline, :template)
      #using CSV
      #CSV.write(joinpath(output_dir,target_subdir * "_template2.csv"),DataFrame("λ"=>spectral_orders_matrix2.λ,"flux_template"=>f_mean2,"var"=>var_mean2, "dfluxdlnλ_template"=>deriv_2,"d²fluxdlnλ²_template"=>deriv2_2))
      using JLD2, FileIO
      save(joinpath(output_dir,target_subdir * "_matrix2.jld2"), Dict("λ"=>spectral_orders_matrix2.λ,"spectra"=>spectral_orders_matrix2.flux,"var_spectra"=>spectral_orders_matrix2.var,"flux_template"=>f_mean2,"var"=>var_mean2, "dfluxdlnλ_template"=>deriv_2,"d²fluxdlnλ²_template"=>deriv2_2))
   end
   dont_need_to!(pipeline, :template)
end

#rvs_dcpca = calc_rvs_from_taylor_expansion(spectral_orders_matrix; mean=f_mean, deriv=deriv)
#std(rvs_dcpca.rv)

#std(rvs_dcpca2.rv)


#need_to!(pipeline, :dcpca)
if need_to(pipeline, :dcpca)
   @assert !need_to(pipeline,:rvs_ccf_total)
   @assert !need_to(pipeline,:template)
   if verbose println("# Performing Doppler constrained PCA analysis.")  end
   using MultivariateStats
   Δrvs_dcpca2 = calc_rvs_from_taylor_expansion(spectral_orders_matrix2; mean=f_mean2, deriv=deriv_2).rv
   dcpca2_out, M2 = RvSpectML.DCPCA.doppler_constrained_pca(spectral_orders_matrix2.flux, deriv_2, Δrvs_dcpca2)
   frac_var_unexplained2 = 1.0.-cumsum(principalvars(M2))./tvar(M2)
   num_basis = length(frac_var_unexplained2)
   println("# Fraction of variance unexplained (orig) = ", frac_var_unexplained[1:min(5,num_basis)])
   println("# Fraction of variance unexplained (clean) = ", frac_var_unexplained2[1:min(5,num_basis)])
   if save_data(pipeline, :dcpca)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_basis2.csv"),  Tables.table(M2.proj) )
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_scores2.csv"), Tables.table(dcpca2_out) )
   end
   dont_need_to!(pipeline, :dcpca)
end


if make_plot(pipeline, :dcpca)  # Ploting results from DCPCA
  using Plots
  plt = scatter(frac_var_unexplained2, xlabel="Number of PCs", ylabel="Frac Variance Unexplained")
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca2_frac_var.png"))   end
  display(plt)
end

if make_plot(pipeline, :dcpca)
   plt_line = 10
   plt = RvSpectML.plot_basis_vectors(order_grids2, f_mean2, deriv_2, M2.proj, idx_plt = spectral_orders_matrix2.chunk_map[plt_line], num_basis=min(4,num_basis), label=plt_line)
   if save_plot(pipeline,:dcpca)
     savefig(plt,joinpath(output_dir, target_subdir * "_dcpca2_basis.png"))
   end
   suspicous_lines = Int[]
   anim = @animate for plt_line ∈ 1:length(order_grids2)
      local idx = spectral_orders_matrix2.chunk_map[plt_line]
      argmin_pixel = argmin(f_mean2[idx])
      if 0.45*length(idx) <= argmin_pixel <= 0.55*length(idx)
        RvSpectML.plot_basis_vectors(order_grids2, f_mean2, deriv_2, M2.proj, idx_plt = spectral_orders_matrix2.chunk_map[plt_line], num_basis=min(4,num_basis), label=plt_line  )
     else
        push!(suspicous_lines,plt_line)
     end
   end
   gif(anim, joinpath(output_dir, target_subdir * "_dcpca2_basis.gif"), fps = 25)
  display(plt)
end

if make_plot(pipeline, :dcpca)
  plt = RvSpectML.plot_basis_scores(order_list_timeseries.times, rvs_ccf_gauss2, dcpca2_out, num_basis=min(4,num_basis) )
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca2_scores.png"))   end
  display(plt)

end

if make_plot(pipeline, :dcpca)
  plt = RvSpectML.plot_basis_scores_cor( rvs_ccf_gauss2, dcpca2_out, num_basis=min(4,num_basis))
  if save_plot(pipeline,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca2_score_cor.png"))   end
  display(plt)
end
