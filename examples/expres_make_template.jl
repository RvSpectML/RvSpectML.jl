cd("RvSpectML")
 using Pkg
 Pkg.activate("examples")

verbose = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
 if verbose   println("# Loading other packages")    end
 using DataFrames, Query, Statistics, Dates

all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_expres_data_101501.jl"))
#all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_neid_solar_data_20190918.jl"))
if typeof(get_inst(all_spectra)) <: AnyEXPRES   continuum_normalize_spectra!(all_spectra)   end

#order_list_timeseries = extract_orders(all_spectra,pipeline_plan,recalc=true)
# Something buggy if include order 1.  Order 86 essentially ignored due to tellurics. First and last orders have NaN issues
# order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=1:85, recalc=true )
order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=12:83, recalc=true )

obsid = 1
chid = 30
extrema(order_list_timeseries[obsid][chid].λ)
using Plots
plot(order_list_timeseries[obsid][chid].λ, order_list_timeseries[obsid][chid].flux)


chunk_λ_grid = RvSpectML.make_grid_for_chunk(order_list_timeseries,chid,oversample_factor=4.0, remove_rv_est=false)
length(chunk_λ_grid)
extrema(chunk_λ_grid)
length(order_list_timeseries[obsid][chid].flux)

res = interp_chunk_to_grid_gp_temporal( order_list_timeseries[obsid][chid], chunk_λ_grid ;   use_logy=false, smooth_factor=16 )
res2 = interp_chunk_to_grid_gp_temporal( order_list_timeseries[obsid][chid], chunk_λ_grid ;   use_logy=false, smooth_factor=1 )
plot(chunk_λ_grid, res.flux)
plot!(chunk_λ_grid, res2.flux)
scatter!(order_list_timeseries[obsid][chid].λ, order_list_timeseries[obsid][chid].flux, markersize=1.3)
xlims!(5080,5085)



xlims!(5080,5085)

# First, make line list using default orders, where laser comb calibration is available
line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=true)
(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true)
line_width_50 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.5)
line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)
line_list_excalibur_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=43:72, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05, recalc=true, verbose=true)
 ((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_excalibur_df, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf = calc_rvs_from_ccf_total(ccfs, ccf_vars, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)

if make_plot(pipeline_plan, :ccf_total)
    using Plots
    t_idx = 20
    #plt = plot(v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
    #scatter!(plt,v_grid,ccfs_expr[:,t_idx]./maximum(ccfs_expr[:,t_idx],dims=1),markersize=1.2,label=:none)
    plt = plot(v_grid,ccfs[:,t_idx],label=:none)
    xlabel!("v (m/s)")
    ylabel!("CCF")
    if save_plot(pipeline_plan,:ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf_sum.png"))   end
    display(plt)
end


# Now repeat, but using custom orders
# Something buggy if include order 1.  Order 86 essentially ignored due to tellurics.  Above 75 tellurics become major issue.  37 is also particularly problematic
#line_list_blue_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=vcat( 12:36, 38:42), v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05, recalc=true, verbose=true)
line_list_blue_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=vcat(32:35,38:39,41:43), v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05, recalc=true, verbose=true)
 ((ccfs_blue, ccf_vars_blue), v_grid) = ccf_total(order_list_timeseries, line_list_blue_df, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
 rvs_ccf_blue = calc_rvs_from_ccf_total(ccfs_blue, ccf_vars_blue, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)

# Compare RVs from red (where LFC is) and selected orders in blue (without LFC)
using Plots
scatter(rvs_ccf.-mean(rvs_ccf), label="red")
 scatter!(rvs_ccf_blue.-mean(rvs_ccf_blue), label="blue")
 scatter!(rvs_ccf_blue.-rvs_ccf.-mean(rvs_ccf_blue.-rvs_ccf), label="diff")

scatter(rvs_ccf.-mean(rvs_ccf), rvs_ccf_blue.-rvs_ccf.-mean(rvs_ccf_blue.-rvs_ccf),label="none")

# If want to check RVs from specific orders
#=
for ord in 12:35
   println("# order = ", ord)
 line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=ord:ord, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05, recalc=true, verbose=true)
 ((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=3)
 rvs_ccf = calc_rvs_from_ccf_total(ccfs, ccf_vars, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
end
=#

need_to!(pipeline_plan, :template)
if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectra.")  end
   @assert !need_to(pipeline_plan,:extract_orders)
   GC.gc()   # run garbage collector for deallocated memory
   map(i->order_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(order_list_timeseries) )
   # Smothing is broken with GP interpolation.  Need to fix.  In mean time, here's a Sinc interpolation workaround
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, smooth_factor=2.0)
   #@time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, alg=:Sinc)
   if save_data(pipeline_plan, :template)
      using JLD2, FileIO
      save(joinpath(output_dir,target_subdir * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
   dont_need_to!(pipeline_plan, :template);
end

if make_plot(pipeline_plan,:template)
   using Plots
   chunkid = 11
   idx = spectral_orders_matrix.chunk_map[chunkid]
   plt = plot(spectral_orders_matrix.λ[idx],f_mean[idx],markersize=1.0,label="Template")
   plot!(plt,spectral_orders_matrix.λ[idx],deriv[idx]./mean(abs.(deriv[idx])),markersize=1.1,label="Deriv")
   plot!(plt,spectral_orders_matrix.λ[idx],deriv2[idx]./mean(abs.(deriv2[idx])),markersize=1.1,label="Deriv 2")
   xlabel!("λ (Å)")
   ylabel!("f(λ)")
   title!("Template spectrum for chunk " * string(chunkid) )
   xlims!(spectral_orders_matrix.λ[idx[floor(Int,1+0.28*length(idx))]],spectral_orders_matrix.λ[idx[floor(Int,1+0.32*length(idx))]])
end

if make_plot(pipeline_plan,:template)
   using Plots
   chunkid = 11
   idx = spectral_orders_matrix.chunk_map[chunkid]
   plt = plot(spectral_orders_matrix.λ[idx],(f_mean[idx].-1.0)./mean(abs.((f_mean[idx].-1.0))),markersize=1.0,label="Template")
   plot!(plt,spectral_orders_matrix.λ[idx],deriv[idx]./mean(abs.(deriv[idx])),markersize=1.1,label="Deriv")
   plot!(plt,spectral_orders_matrix.λ[idx],deriv2[idx]./mean(abs.(deriv2[idx])),markersize=1.1,label="Deriv 2")
   xlabel!("λ (Å)")
   ylabel!("f(λ), f'(λ), f''(λ), all standardized")
   title!("Template spectrum for chunk " * string(chunkid) )
   xlims!(spectral_orders_matrix.λ[idx[floor(Int,1+0.28*length(idx))]],spectral_orders_matrix.λ[idx[floor(Int,1+0.32*length(idx))]])
end

need_to!(pipeline_plan,:fit_lines)
if need_to(pipeline_plan,:fit_lines)
   if verbose println("# Performing fresh search for lines in template spectra.")  end
   cl = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
   # We're done with the spectral_orders_matrix, so we can release the memory now
   #spectral_orders_matrix = nothing
   #GC.gc()
   need_to!(pipeline_plan,:template)
   lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=true)  # TODO: Automate threshold for finding a line

   if verbose println("# Finding above lines in all spectra.")  end
   @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_in_template )

   if save_data(pipeline_plan,:fit_lines)
      using CSV
      CSV.write(joinpath(output_dir,target_subdir * "_linefinder_lines.csv"), lines_in_template )
      CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits.csv"), fits_to_lines )
      #CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits_clean.csv"), lines_to_try )
   end
   dont_need_to!(pipeline_plan,:fit_lines);
end

#=
# TODO: Reevaluate min_deriv2 threshold now that using normalized linear flux rather than unnormalized log flux
lines_in_template_vs_threshold = map(thresh->RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(min_deriv2=thresh)),[0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0])
  map(x->size(x,1),lines_in_template_vs_threshold)
lines_in_template = lines_in_template_vs_threshold[4]
 @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_in_template )
=#

using Printf
function select_line_fits_with_good_std_v(line_fits_df::DataFrame, quantile_threshold::Real; verbose::Bool = false, write_csv::Bool = false)
   fit_distrib = line_fits_df |> @groupby(_.line_id) |>
            @map( { median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
                   std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc),
                   line_id=first(_.line_id),  frac_converged=mean(_.fit_converged)  } ) |>
            @filter(_.frac_converged == 1.0 ) |> DataFrame

   std_v_threshold =  quantile(fit_distrib.std_λc./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps ,quantile_threshold)
   median_σ_width_treshold = 2000 # quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,1-quantile_threshold)
   good_lines_std_v = fit_distrib |>
            @filter( _.std_λc./_.median_λc.*RvSpectML.speed_of_light_mps <= std_v_threshold ) |>
            @filter( sqrt.(_.median_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps >= median_σ_width_treshold ) |>
                  DataFrame
   if verbose
      println("# Found ", size(good_lines_std_v,1), " good lines (std_v), rejected ", size(fit_distrib,1)-size(good_lines_std_v,1), " lines.")
   end
   if write_csv
      val_str = Printf.@sprintf("%1.2f",quantile_threshold)
      CSV.write(joinpath(output_dir,target_subdir * "_good_lines_stdv_quant=" * val_str * ".csv"), good_lines_std_v )
   end
   return good_lines_std_v
end

function select_line_fits_with_good_depth_width_slope(line_fits_df::DataFrame, quantile_threshold::Real; verbose::Bool = false, write_csv::Bool = false )
   fit_distrib = line_fits_df |> @groupby(_.line_id) |>
            @map( { median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
                   std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc),
                   line_id=first(_.line_id),  frac_converged=mean(_.fit_converged)  } ) |>
            @filter(_.frac_converged == 1.0 ) |> DataFrame

   std_depth_treshold = quantile(fit_distrib.std_depth,quantile_threshold)
   median_σ_width_treshold = 2000 # quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,1-quantile_threshold)
   std_σ_width_treshold = quantile(sqrt.(fit_distrib.std_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,quantile_threshold)
   std_b_treshold = quantile(fit_distrib.std_b,quantile_threshold)
   std_a_treshold = quantile(fit_distrib.std_a,quantile_threshold)
   good_lines_alt = fit_distrib |>
      @filter( 0.05 <= _.median_depth <= 0.95 ) |>
      #@filter( sqrt.(_.median_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps >= median_σ_width_treshold ) |>
      @filter( sqrt.(_.std_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps <= std_σ_width_treshold ) |>
      @filter( _.std_b < std_b_treshold) |>

      DataFrame
   if verbose
      println("# Found ", size(good_lines_alt,1), " good lines (std_depth_width_slope), rejected ", size(fit_distrib,1)-size(good_lines_alt,1), " lines.")
   end
   if write_csv
      val_str = Printf.@sprintf("%1.2f",quantile_threshold)
      CSV.write(joinpath(output_dir,target_subdir * "_good_lines_fit_quant=" * val_str * ".csv"), good_lines_alt )
   end
   return good_lines_alt
end


function calc_rvs_from_line_by_line_fits_stdv(fits_to_lines::DataFrame, threshold_stdv::Real)
   good_lines_stdv = select_line_fits_with_good_std_v(fits_to_lines,threshold_stdv,write_csv=false)
   df = fits_to_lines |> @join(good_lines_stdv, _.line_id, _.line_id, {obs_id=_.obs_idx, line_id=_.line_id, Δv=(_.fit_λc.-__.median_λc).*3e8./__.median_λc, weight=(__.median_λc/(3e8*__.std_λc))^2 }) |>
         @groupby(_.obs_id) |> @map({obs_id=first(_.obs_id), Δv= sum(_.Δv.*_.weight)./sum(_.weight) }) |> DataFrame
   rms_stdv = std(df.Δv)
end
function calc_rvs_from_line_by_line_fits_alt(fits_to_lines::DataFrame, threshold_alt::Real)
   good_lines_alt = select_line_fits_with_good_depth_width_slope(fits_to_lines,threshold_alt)
   df2 = fits_to_lines |> @join(good_lines_alt, _.line_id, _.line_id, {obs_id=_.obs_idx, line_id=_.line_id, Δv=(_.fit_λc.-__.median_λc).*3e8./__.median_λc, weight=(__.median_λc/(3e8*__.std_λc))^2 }) |>
            @groupby(_.obs_id) |> @map({obs_id=first(_.obs_id), Δv= sum(_.Δv.*_.weight)./sum(_.weight) }) |> DataFrame
   rms_alt = std(df2.Δv)
end

thresholds = 0.5:0.1:1 #[0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
 plt = plot()
 scatter!(plt,thresholds,map(t->calc_rvs_from_line_by_line_fits_stdv(fits_to_lines,t),thresholds),color=1,label="σ_λ")
 plot!(plt,thresholds,map(t->calc_rvs_from_line_by_line_fits_stdv(fits_to_lines,t),thresholds),color=1,label=:none)
 scatter!(plt,thresholds,map(t->calc_rvs_from_line_by_line_fits_alt(fits_to_lines,t),thresholds),color=2,label="σ_other")
 plot!(plt,thresholds,map(t->calc_rvs_from_line_by_line_fits_alt(fits_to_lines,t),thresholds),color=2,label=:none)
 xlabel!("Threshold for accepting lines")
 ylabel!("RMS v (m/s)")
 savefig("rms_rv_vs_line_accept_threshold.png")
 display(plt)

#= Code below is leftovers from tinkering

map(c->(typeof(c)),fits_to_lines[1,:])
map(c->(typeof(c)),good_lines_stdv[1,:])
good_lines_stdv

println(names(df))
for row in eachrow(df)
   obs_id = row.obs_id
   good_lines_stdv.line_id .== row.
   mean_λc = row.mean_λ
   #line_weights = 1.0 ./ row.std_λ.^2
   Δvs = map(λ->(λ.-mean_λc).*3e8./mean_λc,row.all_λ)
   ave_Δv = sum(Δvs.*lines_weights) / sum(line_weights)
   println(obs_id,ave_Δv)
end

@map( { obs_id=first(_.obs_id),
        median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
               std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc),
               line_id=first(_.line_id),  frac_converged=mean(_.fit_converged)  } ) |>
        @filter(_.frac_converged == 1.0 ) |> DataFrame

if make_plot(pipeline_plan,:template)
   using Plots
   chunkid = 26
   idx = spectral_orders_matrix.chunk_map[chunkid]
   plt = [plot() for i in 1:3]
   for i in 1:3
      local line_plt_idx = findall(x->minimum(spectral_orders_matrix.λ[idx])<=x<=maximum(spectral_orders_matrix.λ[idx]) , good_lines_alt.median_λc)
      map(x->plot!(plt[i],[x,x],[-0.5,0], color=:black,label=:none), good_lines_alt.median_λc[line_plt_idx])
      line_plt_idx = findall(x->minimum(spectral_orders_matrix.λ[idx])<=x<=maximum(spectral_orders_matrix.λ[idx]) , good_lines_stdv.median_λc)
      map(x->plot!(plt[i],[x,x],[ 0,0.5], color=:black,label=:none), good_lines_stdv.median_λc[line_plt_idx])
      scatter!(plt[i],spectral_orders_matrix.λ[idx],(f_mean[idx].-1.0)./maximum(abs.((f_mean[idx].-1.0))),markersize=1.0,color=2,label=:none)
      scatter!(plt[i],spectral_orders_matrix.λ[idx],deriv[idx]./maximum(abs.(deriv[idx])),markersize=1.0,color=3,label=:none)
      scatter!(plt[i],spectral_orders_matrix.λ[idx],deriv2[idx]./maximum(abs.(deriv2[idx])),markersize=1.0,color=4,label=:none)
      ylims!(plt[i],-0.5,0.5)
   end
   λ_div1 = spectral_orders_matrix.λ[idx][floor(Int,length(idx)/3)]
   λ_div2 = spectral_orders_matrix.λ[idx][floor(Int,length(idx)*2/3)]
   xlims!(plt[1],minimum(spectral_orders_matrix.λ[idx]),λ_div1)
   xlims!(plt[2],λ_div1,λ_div2)
   xlims!(plt[3],λ_div2,maximum(spectral_orders_matrix.λ[idx]))
   xlabel!(plt[3],"λ (Å)")
   ylabel!(plt[2],"f(λ), f'(λ), f''(λ), all standardized")
   title!(plt[1],"Template spectrum for chunk " * string(chunkid) )
   plt_combo = plot(plt[1], plt[2], plt[3], layout = (3,1))
   display(plt_combo)
end

# Reset steps of pipeline_plan to rerun with new linelist.
need_to!(pipeline_plan, :ccf_total);
 need_to!(pipeline_plan, :rvs_ccf_total);
 need_to!(pipeline_plan, :ccf_orders);
 need_to!(pipeline_plan, :rvs_ccf_orders);
 need_to!(pipeline_plan, :template);
 need_to!(pipeline_plan, :dcpca);

lines_to_try = good_lines_stdv |> @map({ lambda=_.median_λc, weight=_.median_depth} ) |> DataFrame
 need_to!(pipeline_plan, :ccf_total);
 need_to!(pipeline_plan, :rvs_ccf_total);
 sort!(lines_to_try,:lambda)
 (ccfs2, v_grid2) = ccf_total(order_list_timeseries, lines_to_try, pipeline_plan,  mask_scale_factor=1.0, range_no_mask_change=line_width, ccf_mid_velocity=ccf_mid_velocity, v_step=100, recalc=true)
 alg_fit_rv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=3)
 rvs_ccf2 = calc_rvs_from_ccf_total(ccfs2, pipeline_plan, v_grid=v_grid2, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)

lines_to_try = good_lines_alt |> @map({ lambda=_.median_λc, weight=_.median_depth} ) |> DataFrame
  need_to!(pipeline_plan, :ccf_total);
  need_to!(pipeline_plan, :rvs_ccf_total);
  sort!(lines_to_try,:lambda)
  (ccfs2, v_grid2) = ccf_total(order_list_timeseries, lines_to_try, pipeline_plan,  mask_scale_factor=1.0, range_no_mask_change=line_width, ccf_mid_velocity=ccf_mid_velocity, v_step=100, recalc=true)
  alg_fit_rv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=3)
  rvs_ccf2 = calc_rvs_from_ccf_total(ccfs2, pipeline_plan, v_grid=v_grid2, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)


names(good_lines_stdv)

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



# Random plots that may be useful
histogram(log10.(fit_distrib.std_λc./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps ),nbins=100, alpha=0.5, label="All converged fits")
 good_lines = select_line_fits_with_good_depth_width_slope(fits_to_lines,0.8)
 histogram!(log10.(good_lines.std_λc./good_lines.median_λc.*RvSpectML.speed_of_light_mps),nbins=50,alpha=0.5,label="Consistent Fits")
 good_lines = select_line_fits_with_good_std_v(fits_to_lines,0.7)
 histogram!(log10.(good_lines.std_λc./good_lines.median_λc.*RvSpectML.speed_of_light_mps),nbins=50,alpha=0.5, label="Low σ v Fits")


scatter(log10.(fit_distrib.median_depth),(fit_distrib.median_σ²),markersize=2.0)
 scatter!(log10.(bad_lines.median_depth),(bad_lines.median_σ²),markersize=3.0)

quantile(fit_distrib.std_λc./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps , 0.8)

histogram(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,nbins=200)
 histogram!(sqrt.(good_lines.median_σ²)./good_lines.median_λc.*RvSpectML.speed_of_light_mps,nbins=200,alpha=0.5)

quantile(sqrt.(good_lines.median_σ²)./good_lines.median_λc.*RvSpectML.speed_of_light_mps,0.05)


scatter( sqrt.(fit_distrib.std_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,good_lines.std_λc,markersize=2,label=:none)
 scatter!( sqrt.(good_lines.std_σ²)./good_lines.median_λc.*RvSpectML.speed_of_light_mps,good_lines.std_λc,markersize=3,label=:none)

std_depth_treshold = quantile(fit_distrib.std_depth,0.90)
scatter(good_lines.std_depth,good_lines.std_λc,markersize=2,label=:none)

histogram(good_lines.median_depth,nbins=100)

sqrt.(good_lines.median_σ²)./good_lines.median_λc.*RvSpectML.speed_of_light_mps
#scatter(good_lines.median_depth,good_lines.std_λc,markersize=2,label=:none)




scatter(sqrt.(good_lines.median_σ²)./good_lines.median_λc.*RvSpectML.speed_of_light_mps,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.median_depth,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.median_a,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.median_b,good_lines.std_λc,markersize=2,label=:none)


scatter(sqrt.(good_lines.std_σ²)./good_lines.median_λc.*RvSpectML.speed_of_light_mps,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.std_depth,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.std_a,good_lines.std_λc,markersize=2,label=:none)

scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)


scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)


scatter(good_lines.std_b,good_lines.std_λc,markersize=2,label=:none)




names(lines_in_template)


println(names(lines_in_template_vs_threshold[1]))
=#
