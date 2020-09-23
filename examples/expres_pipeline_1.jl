using Pkg
 Pkg.activate(".")

verbose = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

include("../scripts/scripts.jl")

include("read_expres_101501_data.jl")
order_list_timeseries = extract_orders(all_spectra,pipeline_plan)
line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 21e3)

(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_scale_factor=1.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true)
#(ccfs_expr, v_grid_expr) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan, mask_scale_factor=12.0,  ccf_mid_velocity=ccf_mid_velocity, recalc=true)
#println("# Ratio of max(ccfs_expr)/max(ccfs) = ", mean(maximum(ccfs_expr,dims=1)./maximum(ccfs,dims=1)) )

if make_plot(pipeline_plan, :ccf_total)
   using Plots
   t_idx = 20
   #plt = plot(v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   #scatter!(plt,v_grid,ccfs_expr[:,t_idx]./maximum(ccfs_expr[:,t_idx],dims=1),markersize=1.2,label=:none)
   plt = plot(v_grid,ccfs[:,t_idx],label=:none)
   scatter!(plt,v_grid_expr,ccfs_expr[:,t_idx],markersize=1.2,label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
   if save_plot(pipeline_plan,:ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf_sum.png"))   end
   display(plt)
end


if make_plot(pipeline_plan, :ccf_total)
   include("../scripts/plots/spectra.jl")
   zvals = ccfs./maximum(ccfs,dims=1).-mean(ccfs./maximum(ccfs,dims=1),dims=2)
   colorscale = cgrad(:balance)
   plt = heatmap(v_grid,collect(1:size(ccfs,2)),zvals', c=colorscale, clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
   add_time_gap_lines(plt,order_list_timeseries.times)
   xlabel!("v (m/s)")
   ylabel!("Observation #")
   title!("CCF(v,t)-<CCF>(v) vs time")
   if save_plot(pipeline_plan,:ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf_sum_vs_time_heatmap.png"))   end
   display(plt)
end


rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true)
# Store estimated RVs in metadata for use when making template
map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf[i]-mean(rvs_ccf), 1:length(rvs_ccf) )
#rvs_ccf_expr = calc_rvs_from_ccf_total(ccfs_expr, pipeline_plan, v_grid=v_grid_expr, times = order_list_timeseries.times, recalc=true)

if make_plot(pipeline_plan, :rvs_ccf_total)
   using Plots
   rvs_ccf .-= mean(rvs_ccf)
   #rvs_ccf_expr .-= mean(rvs_ccf_expr)
   plt = scatter(rvs_ccf,markersize=3,label="RVs CCF Std", legend=:bottomright)
   #scatter!(plt,rvs_ccf_expr,markersize=3,label="RVs CCF Expr")
   ylabel!("v (m/s)")
   xlabel!("Time (#)")
   if save_plot(pipeline_plan,:rvs_ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_rvs_ccf_sum.png"))   end
   display(plt)
   #=
   df_yale_resutls = CSV.read(joinpath(homedir(),"Data/EXPRES/inputs/101501/101501_activity.csv"))
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
if make_plot(pipeline_plan, :rvs_ccf_total)
   diff = rvs_ccf.-rvs_ccf_expr
   println(std(diff))
   plt = scatter(order_list_timeseries.times,diff,markersize=4,label="Delta RV")
   ylabel!("Δv (m/s) (Two mask shapes)")
   xlabel!("Time (d)")
   if save_plot(pipeline_plan,:rvs_ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_rvs_ccf_sum.png"))   end
   display(plt)
end

#need_to!(pipeline_plan, :ccf_orders)
order_ccfs = ccf_orders(order_list_timeseries, line_list_df, pipeline_plan)

if make_plot(pipeline_plan, :ccf_orders)
   # order ccfs averaged over observations at multiple times
   obs = 1:length(order_list_timeseries.times)
   order_labels = map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2))
   orders_to_plot = findall(c->sum(order_ccfs[:,c,obs])>0, 1:size(order_ccfs,2))
   zvals =  reshape(sum(order_ccfs[:,orders_to_plot,obs],dims=3)./maximum(sum(order_ccfs[:,orders_to_plot,obs],dims=3),dims=1)   ,size(order_ccfs,1),size(order_ccfs[:,orders_to_plot,obs],2)) .-
            reshape(sum(order_ccfs[:,orders_to_plot,obs],dims=(2,3))./maximum(sum(order_ccfs[:,orders_to_plot,obs],dims=(2,3))),size(order_ccfs,1) )
   plt = heatmap(v_grid,order_labels[orders_to_plot], zvals',c = cgrad(:balance), clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )

   xlabel!("v (m/s)")
   ylabel!("Order ID")
   title!("CCF-<CCF> for obs ID=" * string(obs))
   if save_plot(pipeline_plan,:ccf_orders)   savefig(plt,joinpath(output_dir,target_subdir * " _ccf_orders.png"))   end
   display(plt)
end

#RvSpectML.make_plot!(pipeline_plan, :movie)
if make_plot(pipeline_plan, :ccf_orders) && make_plot(pipeline_plan, :movie)
   # order ccfs averaged over observations at multiple times
   anim = @animate for obs ∈  1:length(order_list_timeseries.times)
      local order_labels = map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2))
      local orders_to_plot = findall(c->sum(order_ccfs[:,c,obs])>0, 1:size(order_ccfs,2))
      local zvals =  reshape(sum(order_ccfs[:,orders_to_plot,obs],dims=3)./maximum(sum(order_ccfs[:,orders_to_plot,obs],dims=3),dims=1)   ,size(order_ccfs,1),size(order_ccfs[:,orders_to_plot,obs],2)) .-
               reshape(sum(order_ccfs[:,orders_to_plot,obs],dims=(2,3))./maximum(sum(order_ccfs[:,orders_to_plot,obs],dims=(2,3))),size(order_ccfs,1) )
      local plt = heatmap(v_grid,order_labels[orders_to_plot], zvals',c = cgrad(:balance), clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )

      xlabel!("v (m/s)")
      ylabel!("Order ID")
      title!("CCF-<CCF> for obs ID=" * string(obs))
   end
   gif(anim, joinpath(output_dir,target_subdir * "_ccf_order_movie_obs.gif"), fps = 5)
end

if make_plot(pipeline_plan, :ccf_orders)
   ord = 3
   zvals = reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,ord,:],2)).-
   reshape(sum(order_ccfs[:,ord,:],dims=(2,3))./maximum(sum(order_ccfs[:,ord,:],dims=(2,3))),size(order_ccfs,1) )
   plt = heatmap(v_grid,1:size(order_ccfs,3), zvals',c = cgrad(:balance), clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
   add_time_gap_lines(plt,order_list_timeseries.times)
   xlabel!("v (m/s)")
   ylabel!("Observation ID")
   title!("CCF-<CCF> for order=" * string(ord))
   if save_plot(pipeline_plan,:ccf_orders)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf_obs_order=" * string(ord) * ".png"))   end
   display(plt)
end

if make_plot(pipeline_plan, :ccf_orders) && make_plot(pipeline_plan, :movie)
   orders_to_plot = findall(c->sum(order_ccfs[:,c,obs])>0, 1:size(order_ccfs,2))
   min_z_val = Inf
   max_z_val = -Inf
   for ord ∈ orders_to_plot
      local zvals = reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,ord,:],2)) .-
                    reshape(sum(order_ccfs[:,ord,:],dims=(2,3))./maximum(sum(order_ccfs[:,ord,:],dims=(2,3))),size(order_ccfs,1) )
      lo, hi = extrema(zvals)
      #println(" ord = ", ord, " lo = ", lo, " hi = ", hi)
      global min_z_val = min(min_z_val,lo)
      global max_z_val = max(max_z_val,hi)
   end
   #println(min_z_val, " - ", max_z_val)
   if -min_z_val > max_z_val
         max_z_val = -min_z_val
   else
         min_z_val = -max_z_val
   end
   #println(min_z_val, " - ", max_z_val)
   anim = @animate for ord ∈ orders_to_plot
      local zvals = reshape(sum(order_ccfs[:,ord,:],dims=3)./maximum(sum(order_ccfs[:,ord,:],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,ord,:],2)).-
      reshape(sum(order_ccfs[:,ord,:],dims=(2,3))./maximum(sum(order_ccfs[:,ord,:],dims=(2,3))),size(order_ccfs,1) )
      local plt = heatmap(v_grid,1:size(order_ccfs,3), zvals',c = cgrad(:balance), clims=(min_z_val,max_z_val) )
      add_time_gap_lines(plt,order_list_timeseries.times)
      xlabel!("v (m/s)")
      ylabel!("Observation ID")
      title!("CCF-<CCF> for order=" * string(order_labels[ord]))
   end
   gif(anim, joinpath(output_dir,target_subdir * "_ccf_order_movie=.gif"), fps = 5)
end

if need_to(pipeline_plan,:scalpels)
   rvs_scalpels = map(n->Scalpels.clean_rvs_scalpels(rvs_ccf_gauss, ccfs, num_basis=n), 1:5)
   println("RMS RVs cleaned by Scalpels: ",std.(rvs_scalpels) )
   dont_need_to!(pipeline_plan,:scalpels)
end

make_plot!(pipeline_plan,:scalpels)
if make_plot(pipeline_plan, :scalpels)
   @assert !need_to(pipeline_plan, :rvs_ccf_total)
   @assert !need_to(pipeline_plan, :ccf_total)
   include("../scripts/plots/scalpels.jl")
   plt = Scalpels.make_plots_scalpels(rvs_ccf_gauss, ccfs, max_num_basis=2, v_grid=v_grid, times=order_list_timeseries.times, output_path="examples/output/solar/figures")
   display(plt)
end


#need_to!(pipeline_plan, :template)
if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectra.")  end
   @assert !need_to(pipeline_plan,:extract_orders)
   @assert !need_to(pipeline_plan,:rvs_ccf_total)
   GC.gc()   # run garbage collector for deallocated memory
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = make_template_spectra(order_list_timeseries)
   if save_data(pipeline_plan, :template)
      #using CSV
      #CSV.write(joinpath(output_dir,target_subdir * "_template.csv"),DataFrame("λ"=>spectral_orders_matrix.λ,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
      using JLD2, FileIO
      save(joinpath(output_dir,target_subdir * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
   dont_need_to!(pipeline_plan, :template);
end

if make_plot(pipeline_plan,:template)
   chunkid = 10
   idx = spectral_orders_matrix.chunk_map[chunkid]
   plt = scatter(spectral_orders_matrix.λ[idx],(f_mean[idx].-1.0)./maximum(abs.((f_mean[idx].-1.0))),markersize=1.0,label=:none)
   scatter!(plt,spectral_orders_matrix.λ[idx],deriv[idx]./maximum(abs.(deriv[idx])),markersize=1.0,label=:none)
   scatter!(plt,spectral_orders_matrix.λ[idx],deriv2[idx]./maximum(abs.(deriv2[idx])),markersize=1.0,label=:none)
   xlabel!("λ (Å)")
   ylabel!("f(λ), f'(λ), f''(λ), all standardized")
   title!("Template spectrum for chunk " * string(chunkid) )
end

#need_to!(pipeline_plan, :dcpca)
if need_to(pipeline_plan, :dcpca)
   if verbose println("# Performing Doppler constrained PCA analysis.")  end
   @assert !need_to(pipeline_plan,:rvs_ccf_total)
   @assert !need_to(pipeline_plan,:template)
   using MultivariateStats
   Δrvs_dcpca = calc_rvs_from_taylor_expansion(spectral_orders_matrix; mean=f_mean, deriv=deriv).rv
   dcpca_out, M = DCPCA.doppler_constrained_pca(spectral_orders_matrix.flux, deriv, Δrvs_dcpca)
   frac_var_unexplained = 1.0.-cumsum(principalvars(M))./tvar(M)
   num_basis = length(frac_var_unexplained)
   println("# Fraction of variance unexplained = ", frac_var_unexplained[1:min(5,num_basis)])
   if save_data(pipeline_plan, :dcpca)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_basis.csv"),  Tables.table(M.proj) )
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_scores.csv"), Tables.table(dcpca_out) )
   end
   dont_need_to!(pipeline_plan, :dcpca);
end


if make_plot(pipeline_plan, :dcpca)  # Ploting results from DCPCA
  using Plots
  # Set parameters for plotting analysis
  plt = scatter(frac_var_unexplained, xlabel="Number of PCs", ylabel="Frac Variance Unexplained")
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca_frac_var.png"))   end
  display(plt)
end

if make_plot(pipeline_plan, :dcpca)
  plt_order = 10
  plt = plot_basis_vectors(order_grids, f_mean, deriv, M.proj, idx_plt = spectral_orders_matrix.chunk_map[plt_order], num_basis=2 )
  #xlims!(5761.5,5766)
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca_basis.png"))   end
  display(plt)
end

if make_plot(pipeline_plan, :dcpca)
  plt = plot_basis_scores(order_list_timeseries.times, rvs_ccf_gauss, dcpca_out, num_basis=min(5,num_basis) )
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca_scores.png"))   end
  display(plt)
end

if make_plot(pipeline_plan, :dcpca)
  plt = plot_basis_scores_cor( rvs_ccf_gauss, dcpca_out, num_basis=min(5,num_basis))
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca_scores_cor.png"))   end
  display(plt)
end

#need_to!(pipeline_plan,:fit_lines)
if need_to(pipeline_plan,:fit_lines)
   if verbose println("# Performing fresh search for lines in template spectra.")  end
   cl = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map))
   #= # We're done with the spectral_orders_matrix, so we can release the memory now
   spectral_orders_matrix = nothing
   GC.gc()
   need_to!(pipeline_plan,:template)
   =#
   lines_in_template = LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(min_deriv2=3))  # TODO: Automate threshold for finding a line

   if verbose println("# Finding above lines in all spectra.")  end
   @time fits_to_lines = LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_in_template )

   if save_data(pipeline_plan,:fit_lines)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_linefinder_lines.csv"), lines_in_template )
      CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits.csv"), fits_to_lines )
      #CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits_clean.csv"), lines_to_try )
   end
   dont_need_to!(pipeline_plan,:fit_lines);
end


#need_to!(pipeline_plan,:clean_line_list_blends)
if need_to(pipeline_plan,:clean_line_list_blends)
   # Exploratory data analysis of the distribution of line properties over time to figure out how to select "clean" lines
   @assert !need_to(pipeline_plan,:fit_lines)
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
      dont_need_to!(pipeline_plan,:clean_line_list_blends);
 end
 fit_distrib

#=if make_plot(pipeline_plan,:fit_lines) || true # Need to add key to pipeline_plan
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

# Reset steps of pipeline_plan to rerun with new linelist.
need_to!(pipeline_plan, :ccf_total);
 need_to!(pipeline_plan, :rvs_ccf_total);
 need_to!(pipeline_plan, :ccf_orders);
 need_to!(pipeline_plan, :rvs_ccf_orders);
 need_to!(pipeline_plan, :template);
 need_to!(pipeline_plan, :dcpca);

#need_to!(pipeline_plan,:ccf_total)
if need_to(pipeline_plan,:ccf_total)
   if verbose println("# Computing CCFs with new line list.")  end
   @assert !need_to(pipeline_plan,:extract_orders)
   @assert !need_to(pipeline_plan,:clean_line_list_blends)
   #mask_shape = CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=tophap_ccf_mask_scale_factor)
   perm = sortperm(lines_to_try.fit_λc)
   line_list2 = CCF.BasicLineList(lines_to_try.fit_λc[perm], lines_to_try.fit_depth[perm] )
   #mask_shape = CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=tophap_ccf_mask_scale_factor*3)
   mask_shape = CCF.GaussianCCFMask(order_list_timeseries.inst, scale_factor=4)
   ccf_plan2 = CCF.BasicCCFPlan(mask_shape = mask_shape, step=100.0, line_list=line_list2, midpoint=0.0)
   v_grid2 = CCF.calc_ccf_v_grid(ccf_plan2)
   @time ccfs2 = CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan2)
   # Write CCFs to file
   if save_data(pipeline_plan,:ccf_total)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_ccfs2.csv"),Tables.table(ccfs2',header=Symbol.(v_grid)))
   end
   dont_need_to!(pipeline_plan,:ccf_total)
end

if make_plot(pipeline_plan, :ccf_total)
   using Plots
   t_idx = 20
   using Plots
   plt = plot(v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   scatter!(plt,v_grid2,ccfs2[:,t_idx]./maximum(ccfs2[:,t_idx],dims=1),markersize=1.2,label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
   if save_plot(pipeline_plan,:ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf2_sum.png"))   end
   display(plt)
end

#need_to!(pipeline_plan,:rvs_ccf_total)
if need_to(pipeline_plan, :rvs_ccf_total)
   if verbose println("# Measuring RVs from CCF.")  end
   @assert !need_to(pipeline_plan,:ccf_total)
   fit_gaussian_to_ccf = RVFromCCF.MeasureRvFromCCFGaussian()
   rvs_ccf_gauss2 = RVFromCCF.measure_rv_from_ccf(v_grid2,ccfs2,alg=fit_gaussian_to_ccf)
   # Store estimated RVs in metadata
   println("# RMS of RVs:           Orig lines ", std(rvs_ccf_gauss),  "    Cleaned lines ", std(rvs_ccf_gauss2))
   rms_rv_nightly2 = bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_gauss2)
   println("# RMS of nightly RVs:   Orig lines ", std(rms_rv_nightly), "    Cleaned lines ", std(rms_rv_nightly2))
   rms_rv_within_night2 = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_gauss2)
   println("# RMS within night RVs: Orig lines ",rms_rv_within_night,  "   Cleaned liens ",rms_rv_within_night2)
   map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf_gauss2[i]-mean(rvs_ccf_gauss2), 1:length(order_list_timeseries) )
   if save_data(pipeline_plan, :rvs_ccf_total)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_rvs_ccf2.csv"),DataFrame("Time [MJD]"=>order_list_timeseries.times,"CCF RV [m/s]"=>rvs_ccf_gauss))
   end
   dont_need_to!(pipeline_plan, :rvs_ccf_total);
end


#need_to!(pipeline_plan, :template)
if need_to(pipeline_plan, :template)
   @assert !need_to(pipeline_plan,:rvs_ccf_total)
   @assert !need_to(pipeline_plan,:clean_line_list_blends)
   # Revert to velocities before cleaning for now
   map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf_gauss[i]-mean(rvs_ccf_gauss), 1:length(order_list_timeseries) )
   #chunk_list_df2 = lines_to_try |> @select(:fit_min_λ,:fit_max_λ) |> @rename(:fit_min_λ=>:lambda_lo, :fit_max_λ=>:lambda_hi) |> DataFrame
   expand_chunk_factor = 4
   chunk_list_df2 = lines_to_try |> @select(:fit_λc,:fit_min_λ,:fit_max_λ) |> @map({:lambda_lo=>_.fit_λc-expand_chunk_factor*(_.fit_λc-_.fit_min_λ), :lambda_hi=>_.fit_λc+expand_chunk_factor*(_.fit_max_λ-_.fit_λc)}) |> DataFrame
   chunk_list_timeseries2 = make_chunk_list_timeseries(spectra,chunk_list_df2)

   # Check that no NaN's included
   #(chunk_list_timeseries2, chunk_list_df2) = filter_bad_chunks(chunk_list_timeseries2,chunk_list_df2)
   chunk_list_timeseries2 = filter_bad_chunks(chunk_list_timeseries2)
   println(size(chunk_list_df2), " vs ", num_chunks(chunk_list_timeseries2) )

   GC.gc()   # run garbage collector for deallocated memory
   if verbose println("# Making template spectra.")  end
   @time ( spectral_orders_matrix2, f_mean2, var_mean2, deriv_2, deriv2_2, order_grids2 )  = make_template_spectra(chunk_list_timeseries2)

   if save_data(pipeline_plan, :template)
      #using CSV
      #CSV.write(joinpath(output_dir,target_subdir * "_template2.csv"),DataFrame("λ"=>spectral_orders_matrix2.λ,"flux_template"=>f_mean2,"var"=>var_mean2, "dfluxdlnλ_template"=>deriv_2,"d²fluxdlnλ²_template"=>deriv2_2))
      using JLD2, FileIO
      save(joinpath(output_dir,target_subdir * "_matrix2.jld2"), Dict("λ"=>spectral_orders_matrix2.λ,"spectra"=>spectral_orders_matrix2.flux,"var_spectra"=>spectral_orders_matrix2.var,"flux_template"=>f_mean2,"var"=>var_mean2, "dfluxdlnλ_template"=>deriv_2,"d²fluxdlnλ²_template"=>deriv2_2))
   end
   dont_need_to!(pipeline_plan, :template)
end

#rvs_dcpca = calc_rvs_from_taylor_expansion(spectral_orders_matrix; mean=f_mean, deriv=deriv)
#std(rvs_dcpca.rv)

#std(rvs_dcpca2.rv)


need_to!(pipeline_plan, :dcpca)
 if need_to(pipeline_plan, :dcpca)
   @assert !need_to(pipeline_plan,:rvs_ccf_total)
   @assert !need_to(pipeline_plan,:template)
   if verbose println("# Performing Doppler constrained PCA analysis.")  end
   using MultivariateStats
   # Second EXPRES observation of 101501 is very weird. Let's try leaving it out of the DCPCA analysis.
   obs_incl_for_dcpca = vcat([1], collect(3:size(spectral_orders_matrix2.flux,2)) )
   spectral_orders_matrix2 = make_spectral_time_series_common_wavelengths_with_selected_times(spectral_orders_matrix2, obs_incl_for_dcpca)
   Δrvs_dcpca2 = calc_rvs_from_taylor_expansion(spectral_orders_matrix2; mean=f_mean2, deriv=deriv_2).rv
   dcpca2_out, M2 = DCPCA.doppler_constrained_pca(spectral_orders_matrix2.flux, deriv_2, Δrvs_dcpca2)
   frac_var_unexplained2 = 1.0.-cumsum(principalvars(M2))./tvar(M2)
   num_basis = length(frac_var_unexplained2)
   println("# Fraction of variance unexplained (orig) = ", frac_var_unexplained[1:min(5,length(frac_var_unexplained))])
   println("# Fraction of variance unexplained (clean) = ", frac_var_unexplained2[1:min(5,num_basis)])
   if save_data(pipeline_plan, :dcpca)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_basis2.csv"),  Tables.table(M2.proj) )
      CSV.write(joinpath(output_dir,target_subdir * "_dcpca_scores2.csv"), Tables.table(dcpca2_out) )
   end
   dont_need_to!(pipeline_plan, :dcpca)
end


if make_plot(pipeline_plan, :dcpca)  # Ploting results from DCPCA
  using Plots
  plt = scatter(frac_var_unexplained2, xlabel="Number of PCs", ylabel="Frac Variance Unexplained")
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca2_frac_var.png"))   end
  display(plt)
end

if make_plot(pipeline_plan, :dcpca)
   plt_line = 10
   plt = plot_basis_vectors(order_grids2, f_mean2, deriv_2, M2.proj, idx_plt = spectral_orders_matrix2.chunk_map[plt_line], num_basis=min(4,num_basis), label=plt_line)
   if save_plot(pipeline_plan,:dcpca)
     savefig(plt,joinpath(output_dir,target_subdir * "_dcpca2_basis.png"))
   end
   # Making animations is slow.  # TODO: Add animations key to pipeline_plan?
   #= suspicous_lines = Int[]
   anim = @animate for plt_line ∈ 1:length(order_grids2)
      local idx = spectral_orders_matrix2.chunk_map[plt_line]
      argmin_pixel = argmin(f_mean2[idx])
      if 0.45*length(idx) <= argmin_pixel <= 0.55*length(idx)
        plot_basis_vectors(order_grids2, f_mean2, deriv_2, M2.proj, idx_plt = spectral_orders_matrix2.chunk_map[plt_line], num_basis=min(4,num_basis), label=plt_line  )
     else
        push!(suspicous_lines,plt_line)
     end
   end
   gif(anim, joinpath(output_dir,"dcpca2_basis.gif"), fps = 25)
   =#

  display(plt)
end

if make_plot(pipeline_plan, :dcpca)
  plt = plot_basis_scores(order_list_timeseries.times[obs_incl_for_dcpca], rvs_ccf_gauss2[obs_incl_for_dcpca], dcpca2_out, num_basis=min(4,num_basis) )
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca2_scores.png"))   end
  display(plt)

end

if make_plot(pipeline_plan, :dcpca)
  plt = plot_basis_scores_cor( rvs_ccf_gauss2[obs_incl_for_dcpca], dcpca2_out, num_basis=min(4,num_basis))
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca2_score_cor.png"))   end
  display(plt)
end
