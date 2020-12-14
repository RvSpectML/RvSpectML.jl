using Pkg
 Pkg.activate(".")

cd("RvSpectML")
Pkg.activate("examples")
verbose = true
  make_plots = true
  if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
  using RvSpectML
  if verbose && make_plots && !isdefined(Main,:RvSpectMLPlots)  println("# Loading RvSpecMLPlots")    end
  if make_plots
     using RvSpectMLPlots
  end
  if verbose   println("# Loading other packages")    end
  using Statistics

all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_expres_data_101501.jl"))
#all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_neid_solar_data_20190918.jl"))
if typeof(get_inst(all_spectra)) <: AnyEXPRES   continuum_normalize_spectra!(all_spectra)   end

order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=12:85, recalc=true )

function calc_depth_and_expected_rv_precission(spectrum::ST, pixels::AR, order::Integer) where { ST<:AbstractSpectra2D, AR<:AbstractRange{Int64}, AA1<:AbstractArray{Int64,1} }
  logλ = log.(view(spectrum.λ,pixels,order))
  flux = view(spectrum.flux,pixels,order)
  var = view(spectrum.var,pixels,order)
  (f_mean, f_deriv) = RvSpectML.TemporalGPInterpolation.predict_mean_and_deriv(logλ, flux, logλ;sigmasq_obs=var, use_logx=false, use_logy=false, smooth_factor=2)
  depth = 1.0 .- minimum(f_mean)/maximum(f_mean)
  exp_sigma_rv = RvSpectMLBase.speed_of_light_mps / sqrt( sum(f_deriv.^2 ./ var) )
  return (depth=depth, exp_σ_rv=exp_sigma_rv)
end

#order_info = get_order_info(all_spectra, orders_to_use=1:85)

#if !make_plots     make_no_plots!(pipeline_plan)    end
line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  orders_to_use=38:76, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true)
#line_list_orig = deepcopy(line_list_df)

#line_list_merged = RvSpectML.assign_lines_to_orders(line_list_df, order_info)
#line_list_reweighted = RvSpectML.calc_snr_weights_for_lines!(deepcopy(line_list_merged), all_spectra)
#result = calc_weights_for_lines(copy(line_list_merged), all_spectra )

# range_no_mask_change=RvSpectMLBase.max_bc, step=v_step,  max=RvSpectMLBase.max_bc

(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:tophat, mask_scale_factor=7.0,
  v_step=400, ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=false)
rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true).-ccf_mid_velocity

#alg_fit_rv = EchelleCCFs.MeasureRvFromCCFQuadratic()
#alg_fit_rv = EchelleCCFs.MeasureRvFromCCFTemplate(v_grid=v_grid2, template=vec(mean(ccfs2,dims=2)), frac_of_width_to_fit=1.0, measure_width_at_frac_depth=0.5)
alg_fit_rv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=1.0, measure_width_at_frac_depth=0.5, init_guess_ccf_σ=ccf_mid_velocity)
((ccfs2, ccf_vars2), v_grid2) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:tophat, mask_scale_factor=7.0,
      ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=true)
rvs_ccf2 = calc_rvs_from_ccf_total(ccfs2, ccf_vars2, pipeline_plan, v_grid=v_grid2, times = order_list_timeseries.times, alg_fit_rv=alg_fit_rv, recalc=true)

using Plots
plot(v_grid,ccfs./maximum(ccfs,dims=1),color=:blue,label=:none)
 plot!(v_grid2,ccfs2./maximum(ccfs2,dims=1),markersize=1.2, color=:red,label=:none)





x_plt = order_list_timeseries.times # mod.(order_list_timeseries.times,365.25)
scatter(x_plt,rvs_ccf.-mean(rvs_ccf), color=:blue); scatter!(x_plt,rvs_ccf2.-mean(rvs_ccf2),color=:red)

(ccfs2, v_grid2) = ccf_total(order_list_timeseries, line_list_reweighted, pipeline_plan,  mask_type=:tophat, mask_scale_factor=2.0,
  ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=false)

((ccfs2, ccf_vars2), v_grid2) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:halfcos, mask_scale_factor=11.0,
    ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=true)


rvs_ccf2 = calc_rvs_from_ccf_total(ccfs2, ccf_vars2, pipeline_plan, v_grid=v_grid2, times = order_list_timeseries.times, alg_fit_rv=alg_fit_rv, recalc=true)


#scatter(line_list_reweighted[:depth]./line_list_reweighted[:weight_snr],line_list_reweighted[:weight], markersize=1.2)
scatter(line_list_reweighted[!,:depth],line_list_reweighted[!,:weight], markersize=1.2)
scatter(1.0 ./line_list_reweighted[!,:weight_snr],line_list_reweighted[!,:weight], markersize=1.2)
scatter(line_list_reweighted[!,:depth].*line_list_reweighted[!,:weight_snr],line_list_reweighted[!,:weight], markersize=1.2)
scatter(line_list_reweighted[!,:weight_snr] ./line_list_reweighted[!,:exp_sigma_rv].^2,line_list_reweighted[!,:weight], markersize=1.2)



#snr_of_spectra = vec(mean(result[:,:depth],dims=2))
#plot((result./snr_of_spectra)', label=:none)

plot(v_grid2,ccfs2./maximum(ccfs2,dims=2).-mean(ccfs2./maximum(ccfs2,dims=2),dims=2),markersize=1.2, color=:blue,label=:none)

ref_obs_idx = RvSpectMLBase.choose_obs_idx_for_init_guess(df_files_use,first(all_spectra).inst)

line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,ref_obs_idx),frac_depth=0.05)

line_list_df = prepare_line_list(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05)

((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:tophat, mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true,
      calc_ccf_var=true)
# Code working towards propagating uncertainties through CCF
#((ccfs2, ccf_vars2), v_grid2) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=:gaussian, mask_scale_factor=10.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=true)

rvs_ccf = calc_rvs_from_ccf_total(ccfs, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true)

if make_plot(pipeline_plan, :ccf_total)
   using Plots
   plt = RvSpectMLPlots.make_plot_ccf_vs_time(ccfs,v_grid = v_grid, output_path=output_dir, save_fig = save_plot(pipeline_plan,:ccf_total) )
   #plt = RvSpectMLPlots.make_plot_ccf_vs_time(ccfs2,v_grid = v_grid2, output_path=output_dir, save_fig = save_plot(pipeline_plan,:ccf_total) )
   display(plt)
end

if make_plot(pipeline_plan, :ccf_total)
   plt = RvSpectMLPlots.make_heatmap_ccf_vs_time(ccfs,v_grid = v_grid, times=order_list_timeseries.times, output_path=output_dir, save_fig = save_plot(pipeline_plan,:ccf_total) )
   #plt = RvSpectMLPlots.make_heatmap_ccf_vs_time(ccfs2,v_grid = v_grid2, times=order_list_timeseries.times, output_path=output_dir, save_fig = save_plot(pipeline_plan,:ccf_total) )
   display(plt)
end



(order_ccfs, v_grid_order_ccfs_novars) = ccf_orders(order_list_timeseries, line_list_df, pipeline_plan, orders_to_use=1:85, ccf_mid_velocity=ccf_mid_velocity, recalc=true)

plot(order_ccfs_novars[:,:,1],label=:none)
scatter!(order_ccfs_novars[:,5,1],markersize=1.2,label=:none)



mask_shape = EchelleCCFs.TopHatCCFMask(order_list_timeseries.inst) #, scale_factor=8)
#line_list = EchelleCCFs.BasicLineList(line_list_df.lambda, line_list_df.weight)
line_list = EchelleCCFs.BasicLineList2D(line_list_df)
v_step = 400
ccf_plan = EchelleCCFs.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=RvSpectMLBase.max_bc, step=v_step,  max=RvSpectMLBase.max_bc  )
v_grid_orders_ccfs = EchelleCCFs.calc_ccf_v_grid(ccf_plan)
#(order_ccfs, order_ccf_vars )
 order_ccfs = EchelleCCFs.calc_order_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)

ccf_plan2 = EchelleCCFs.BasicCCFPlan(mask_shape = mask_shape, line_list=EchelleCCFs.BasicLineList2D(line_list_df), midpoint=ccf_mid_velocity, range_no_mask_change=RvSpectMLBase.max_bc, step=v_step,  max=RvSpectMLBase.max_bc )
(order_ccfs2, order_ccf_vars2 ) = EchelleCCFs.calc_order_ccf_and_var_chunklist_timeseries(order_list_timeseries, ccf_plan2)


plt = plot()
 obsid = 9
 order_ccfs = order_ccfs_novars
 for chid in 1:size(order_ccfs,2)
   plot!(plt,v_grid_orders_ccfs,order_ccfs[:,chid,1]./maximum(order_ccfs[:,chid,1]),label=:none)
   #scatter!(plt,v_grid_orders_ccfs,order_ccfs2[:,chid,1]./maximum(order_ccfs2[:,chid,1]),markersize=1.2,label=:none)
   #plot!(plt,v_grid_orders_ccfs,sqrt.(order_ccf_vars[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),label=:none)

 end
 display(plt)


chid = 84; scatter!(plt,v_grid_orders_ccfs,order_ccfs[:,chid,1]./maximum(order_ccfs[:,chid,1]),markersize=1.3, label=:none)

plt = plot()
  chid = 52
  mean_for_order = vec(mean(order_ccfs[:,chid,:]./maximum(order_ccfs[:,chid,:],dims=1),dims=2))
  for obsid in 1:size(order_ccfs,3)
    plot!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]).-mean_for_order .+ 0.003*obsid,label=:none)
    #if obsid ==1
  #    scatter!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),yerr=sqrt.(order_ccf_vars[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),label=:none)
  #  end
  end
  display(plt)
  size(mean_for_order)

ccf_comp  = reshape(sum(order_ccfs,dims=2),size(order_ccfs,1),size(order_ccfs,3))

ccfs .-ccf_comp

scatter(map(chid->minimum(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),1:85))

plt = plot()
    chid = 50
    for obsid in 1:size(order_ccfs,3)
      plot!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]) .-
                     mean(order_ccfs[:,chid,:]./maximum(order_ccfs[:,chid,:]),dims=2) ,label=:none)
      #if obsid ==1
      #  scatter!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),yerr=sqrt.(order_ccf_vars[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]),label=:none)
      #end
    end
    display(plt)
obsid = 6; scatter!(plt,v_grid_orders_ccfs,(order_ccfs[:,chid,obsid])./maximum(order_ccfs[:,chid,obsid]) .-
               mean(order_ccfs[:,chid,:]./maximum(order_ccfs[:,chid,:]),dims=2) , markersize=1.2, label=string(obsid))

if make_plot(pipeline_plan, :ccf_orders)
   plt = make_heatmap_ccf_vs_order(order_ccfs, v_grid=v_grid_orders_ccfs, order_labels=map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2)),
                      output_path=output_dir * target_subdir, save_fig=save_plot(pipeline_plan,:ccf_orders)   )
   display(plt)
end

if make_plot(pipeline_plan, :ccf_orders)
   plt = make_heatmap_ccf_vs_order(order_ccfs2, v_grid=v_grid_orders_ccfs, order_labels=map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2)),
                      output_path=output_dir * target_subdir, save_fig=save_plot(pipeline_plan,:ccf_orders)   )
   display(plt)
end

order_ccfs.-order_ccfs2

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

if true || need_to(pipeline_plan,:scalpels)
   rvs_scalpels2 = map(n->Scalpels.clean_rvs_scalpels(rvs_ccf, ccfs2, num_basis=n), 1:5)
   #rvs_scalpels2 = map(n->Scalpels.clean_rvs_scalpels(rvs_ccf, reshape(order_ccfs[:,order_to_use,:],155*25,44), num_basis=n), 1:5)
   println("RMS RVs cleaned by Scalpels: ",std.(rvs_scalpels2) )
   dont_need_to!(pipeline_plan,:scalpels)
end

make_plot!(pipeline_plan,:scalpels)
if true || make_plot(pipeline_plan, :scalpels)
   @assert !need_to(pipeline_plan, :rvs_ccf_total)
   @assert !need_to(pipeline_plan, :ccf_total)
   plt = make_plots_scalpels(rvs_ccf2, ccfs2, max_num_basis=2, v_grid=v_grid2, times=order_list_timeseries.times, output_path=output_dir, save_fig = save_plot(pipeline_plan, :scalpels) )
   display(plt)
end
