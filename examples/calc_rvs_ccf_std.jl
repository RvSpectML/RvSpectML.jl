verbose = true
  if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
  using RvSpectML
  if verbose && !isdefined(Main,:RvSpectMLPlots)  println("# Loading RvSpecMLPlots")    end
  using RvSpectMLPlots
  if verbose   println("# Loading other packages")    end
  using Statistics

all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_expres_data_101501.jl"))
#all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_neid_solar_data_20190918.jl"))

order_list_timeseries = extract_orders(all_spectra,pipeline_plan, recalc=true )

line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 30e3, recalc=true)

(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:tophat, mask_scale_factor=2.0,
  ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=false)

line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)

line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05)

((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:tophat, mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=true)
# Code working towards propagating uncertainties through CCF
#((ccfs2, ccf_vars2), v_grid2) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:gaussian, mask_scale_factor=10.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=true)

rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true)

if make_plot(pipeline_plan, :ccf_total)
   using Plots
   t_idx = 1:length(all_spectra)
   plt = plot()
   plot!(plt,v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   xlabel!("v (m/s)")
   ylabel!("CCF")
   if save_plot(pipeline_plan,:ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_ccf_sum.png"))   end
   display(plt)
end

if make_plot(pipeline_plan, :ccf_total)
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

(order_ccfs, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_df, pipeline_plan)

if need_to(pipeline_plan,:scalpels)
   rvs_scalpels = map(n->Scalpels.clean_rvs_scalpels(rvs_ccf, ccfs, num_basis=n), 1:5)
   println("RMS RVs cleaned by Scalpels: ",std.(rvs_scalpels) )
   dont_need_to!(pipeline_plan,:scalpels)
end

make_plot!(pipeline_plan,:scalpels)
if make_plot(pipeline_plan, :scalpels)
   @assert !need_to(pipeline_plan, :rvs_ccf_total)
   @assert !need_to(pipeline_plan, :ccf_total)
   plt = make_plots_scalpels(rvs_ccf, ccfs, max_num_basis=2, v_grid=v_grid, times=order_list_timeseries.times, output_path="examples/output/figures")
   display(plt)
end
