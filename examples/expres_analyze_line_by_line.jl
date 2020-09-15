using Pkg
 Pkg.activate(".")

verbose = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
 include("shared/scripts.jl")
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

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

#need_to!(pipeline_plan,:fit_lines)
if need_to(pipeline_plan,:fit_lines)
   if verbose println("# Performing fresh search for lines in template spectra.")  end
   cl = ChunkList(map(grid->ChuckOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map))
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

#scatter(good_lines.median_depth,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.std_depth,good_lines.std_λc,markersize=2,label=:none)
#scatter(good_lines.std_σ²,good_lines.std_λc,markersize=2,label=:none)
#scatter(good_lines.median_σ²,good_lines.std_λc,markersize=2,label=:none)
#scatter(good_lines.median_b,good_lines.std_λc,markersize=2,label=:none)
#scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)
#scatter(good_lines.median_b,good_lines.std_λc,markersize=2,label=:none)
#scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)
