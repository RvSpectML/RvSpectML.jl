using Pkg
 Pkg.activate(".")

verbose = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
 include("shared/scripts.jl")
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

include("shared/scripts.jl")
include("read_expres_data_101501.jl")

order_list_timeseries = extract_orders(all_spectra,pipeline_plan)

line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,
         v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 21e3)

#need_to!(pipeline_plan, :template)
if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectra.")  end
   @assert !need_to(pipeline_plan,:extract_orders)
#   @assert !need_to(pipeline_plan,:rvs_ccf_total)
   GC.gc()   # run garbage collector for deallocated memory
   map(i->order_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(order_list_timeseries) )
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
   using Plots
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
   #@assert !need_to(pipeline_plan,:rvs_ccf_total)
   @assert !need_to(pipeline_plan,:template)
   using MultivariateStats
   Δrvs_dcpca = calc_rvs_from_taylor_expansion(spectral_orders_matrix; mean=f_mean, deriv=deriv).rv
   #map(i->order_list_timeseries.metadata[i][:rv_est] = Δrvs_dcpca[i], 1:length(order_list_timeseries) )
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
  include("../scripts/plots/dcpca.jl")
  #using Plots
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
  plt = plot_basis_scores(order_list_timeseries.times, Δrvs_dcpca, dcpca_out, num_basis=min(5,num_basis) )
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca_scores.png"))   end
  display(plt)
end

if make_plot(pipeline_plan, :dcpca)
  plt = plot_basis_scores_cor( Δrvs_dcpca, dcpca_out, num_basis=min(5,num_basis))
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

# Reset steps of pipeline_plan to rerun with new linelist.
need_to!(pipeline_plan, :ccf_total);
 need_to!(pipeline_plan, :rvs_ccf_total);
 need_to!(pipeline_plan, :ccf_orders);
 need_to!(pipeline_plan, :rvs_ccf_orders);
 need_to!(pipeline_plan, :template);
 need_to!(pipeline_plan, :dcpca);

RvSpectML.dont_save_data!(pipeline_plan,:ccf_total)
idx = sortperm(lines_to_try.fit_λc)
line_list_clean_df = DataFrame(:lambda=>lines_to_try.fit_λc[idx], :weight=>lines_to_try.fit_depth[idx] )
(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_scale_factor=1.0,
                           ccf_mid_velocity=ccf_mid_velocity, recalc=true)
(ccfs2, v_grid2) = ccf_total(order_list_timeseries, line_list_clean_df, pipeline_plan,  mask_scale_factor=1.0,
                           ccf_mid_velocity=0, recalc=true)

if make_plot(pipeline_plan, :ccf_total)
   using Plots
   t_idx = 20
   using Plots
   plt = plot(v_grid.-ccf_mid_velocity,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   scatter!(plt,v_grid2,ccfs2[:,t_idx]./maximum(ccfs2[:,t_idx],dims=1),markersize=1.2,label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
   if save_plot(pipeline_plan,:ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf2_sum.png"))   end
   display(plt)
end

#need_to!(pipeline_plan,:rvs_ccf_total)
# TODO: Replace using ccf's here with calc_rvs_from_taylor_expansion using new template?
# Δrvs_dcpca = calc_rvs_from_taylor_expansion(spectral_orders_matrix; mean=f_mean, deriv=deriv).rv

rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times=order_list_timeseries.times,recalc=true)
rvs_ccf2 = calc_rvs_from_ccf_total(ccfs2, pipeline_plan, v_grid=v_grid2, times=order_list_timeseries.times,recalc=true)

#need_to!(pipeline_plan, :template)
if need_to(pipeline_plan, :template)
   @assert !need_to(pipeline_plan,:rvs_ccf_total)
   @assert !need_to(pipeline_plan,:clean_line_list_blends)
   # Revert to velocities before cleaning for now
   map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf2[i]-mean(rvs_ccf2), 1:length(order_list_timeseries) )
   #chunk_list_df2 = lines_to_try |> @select(:fit_min_λ,:fit_max_λ) |> @rename(:fit_min_λ=>:lambda_lo, :fit_max_λ=>:lambda_hi) |> DataFrame
   expand_chunk_factor = 4
   chunk_list_df2 = lines_to_try |> @select(:fit_λc,:fit_min_λ,:fit_max_λ) |> @map({:lambda_lo=>_.fit_λc-expand_chunk_factor*(_.fit_λc-_.fit_min_λ), :lambda_hi=>_.fit_λc+expand_chunk_factor*(_.fit_max_λ-_.fit_λc)}) |> DataFrame
   chunk_list_timeseries2 = make_chunk_list_timeseries(all_spectra,chunk_list_df2)

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
  plt = plot_basis_scores(order_list_timeseries.times[obs_incl_for_dcpca], rvs_ccf2[obs_incl_for_dcpca], dcpca2_out, num_basis=min(4,num_basis) )
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir, target_subdir * "_dcpca2_scores.png"))   end
  display(plt)

end

if make_plot(pipeline_plan, :dcpca)
  plt = plot_basis_scores_cor( rvs_ccf2[obs_incl_for_dcpca], dcpca2_out, num_basis=min(4,num_basis))
  if save_plot(pipeline_plan,:dcpca)   savefig(plt,joinpath(output_dir,target_subdir * "_dcpca2_score_cor.png"))   end
  display(plt)
end
