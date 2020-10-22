verbose = true
  make_plots = false
  if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
  using RvSpectML
  if verbose && make_plots && !isdefined(Main,:RvSpectMLPlots)  println("# Loading RvSpecMLPlots")    end
  if make_plots    using RvSpectMLPlots    end
  if verbose   println("# Loading other packages")    end
  using Statistics

# Load data
all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_expres_data_101501.jl"))
#all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_neid_solar_data_20190918.jl"))
if typeof(get_inst(all_spectra)) <: AnyEXPRES   continuum_normalize_spectra!(all_spectra)   end
order_list_timeseries = extract_orders(all_spectra, pipeline_plan, orders_to_use=orders_to_use_default(first(all_spectra).inst), recalc=true )

# Preliminary CCF calculation
line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=orders_to_use_default(first(all_spectra).inst), v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true)
(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:tophat, mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=false)
rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true)
if abs(ccf_mid_velocity-mean(rvs_ccf)) > 1e3     ccf_mid_velocity=mean(rvs_ccf)   end

# Measure how much space needed to ensure no telluric contamination from from one high-SNR spectra and recompute line list
ref_obs_idx = RvSpectMLBase.choose_obs_idx_for_init_guess(df_files_use,first(all_spectra).inst)
line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,ref_obs_idx),frac_depth=0.05)
line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=orders_to_use_default(first(all_spectra).inst), v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics=line_width_05, recalc=true)

# Recompute CCF using wider mask (check results from expres_compare_mask_widths for mask suggestions)
((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:halfcos, mask_scale_factor=10.0, ccf_mid_velocity=ccf_mid_velocity, calc_ccf_var=true, recalc=true )
alg_fit_rv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=1.0, measure_width_at_frac_depth=0.5, init_guess_ccf_σ=ccf_mid_velocity)
 rvs_ccf = calc_rvs_from_ccf_total(ccfs, ccf_vars, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, alg_fit_rv=alg_fit_rv, recalc=true)

if !make_plots     make_no_plots!(pipeline_plan)    end
if make_plot(pipeline_plan, :ccf_total)
   using Plots
   using RvSpectMLPlots
   plt = RvSpectMLPlots.make_plot_ccf_vs_time(ccfs,v_grid = v_grid, output_path=output_dir, save_fig = save_plot(pipeline_plan,:ccf_total) )
   display(plt)
end

if make_plot(pipeline_plan, :ccf_total)
   plt = RvSpectMLPlots.make_heatmap_ccf_vs_time(ccfs,v_grid = v_grid, times=order_list_timeseries.times, output_path=output_dir, save_fig = save_plot(pipeline_plan,:ccf_total) )
   display(plt)
end


need_to!(pipeline_plan,:scalpels)
if need_to(pipeline_plan,:scalpels)
   rvs_scalpels = map(n->Scalpels.clean_rvs_scalpels(rvs_ccf, ccfs, num_basis=n), 1:5)
   println("RMS RVs cleaned by Scalpels: ",std.(rvs_scalpels) )
   dont_need_to!(pipeline_plan,:scalpels)
end

make_plot!(pipeline_plan,:scalpels)
if make_plot(pipeline_plan, :scalpels)
   @assert !need_to(pipeline_plan, :rvs_ccf_total)
   @assert !need_to(pipeline_plan, :ccf_total)
   plt = make_plots_scalpels(rvs_ccf, ccfs, max_num_basis=2, v_grid=v_grid, times=order_list_timeseries.times, output_path=output_dir, save_fig = save_plot(pipeline_plan, :scalpels) )
   display(plt)
end

# Continue, if want to compute and plot order CCFs (e.g., to figure out which to include/reject)
#(order_ccfs_novars, v_grid_order_ccfs_novars) = ccf_orders(order_list_timeseries, line_list_df, pipeline_plan,  ccf_mid_velocity=ccf_mid_velocity, recalc=true)
mask_shape = EchelleCCFs.TopHatCCFMask(order_list_timeseries.inst, scale_factor=2)
line_list = EchelleCCFs.BasicLineList2D(line_list_df)
v_step = 250
ccf_plan = EchelleCCFs.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=RvSpectMLBase.max_bc, #= step=v_step,=#  max=RvSpectMLBase.max_bc  )
v_grid_orders_ccfs = EchelleCCFs.calc_ccf_v_grid(ccf_plan)
(order_ccfs, order_ccf_vars ) = EchelleCCFs.calc_order_ccf_and_var_chunklist_timeseries(order_list_timeseries, ccf_plan)
#all(order_ccfs.-order_ccfs_novars .== 0)

plt = plot()
 obsid = 1
 for chid in 1:size(order_ccfs,2)
   plot!(plt,v_grid_orders_ccfs,order_ccfs[:,chid,1],label=:none)
 end
 display(plt)

 plt = plot()
  obsid = 1
  for chid in 1:size(order_ccfs,2)
    plot!(plt,v_grid_orders_ccfs,order_ccfs[:,chid,1]./maximum(order_ccfs[:,chid,1]),label=:none)
  end
  display(plt)

plt = plot()
  chid = 1
  for obsid in 1:size(order_ccfs,3)
    plot!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),label=:none)
    if obsid ==1
    #  scatter!(plt,v_grid_orders_ccfs,(order_ccfs_novars[:,chid,obsid])./maximum(order_ccfs_novars[:,chid,obsid]),markersize=1.3,label=:none)
    end
  end
  display(plt)


plt = plot()
  chid = 10
  for obsid in 1:size(order_ccfs,3)
    plot!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]) .-
                   mean(order_ccfs[:,chid,:]./maximum(order_ccfs[:,chid,:]),dims=2) ,label=:none)
  end
  begin
    obsid = 1
    scatter!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]) .- mean(order_ccfs[:,chid,:]./maximum(order_ccfs[:,chid,:]),dims=2), markersize=1.3, label=:none)
  end
  display(plt)

if make_plot(pipeline_plan, :ccf_orders)
   plt = make_heatmap_ccf_vs_order(order_ccfs, v_grid=v_grid_orders_ccfs, order_labels=map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2)),
                      output_path=output_dir * target_subdir, save_fig=save_plot(pipeline_plan,:ccf_orders)   )
   display(plt)
end
