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

line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 30e3)

(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_scale_factor=10.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true)

line_width = RvSpectML.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)

line_list_df = prepare_line_list_pass1(linelist_for_ccf_filename, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width)

if make_plot(pipeline_plan, :ccf_total)
   using Plots
   t_idx = 20
   plt = plot(v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
   if isdefined(Main,:ccfs_expr)   scatter!(plt,v_grid,ccfs_expr[:,t_idx]./maximum(ccfs_expr[:,t_idx],dims=1),markersize=1.2,label=:none)  end
   #plt = plot(v_grid,ccfs[:,t_idx],label=:none)
   #if isdefined(Main,:ccfs_expr)   scatter!(plt,v_grid_expr,ccfs_expr[:,t_idx],markersize=1.2,label=:none)   end
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


if need_to(pipeline_plan,:scalpels)
   rvs_scalpels = map(n->Scalpels.clean_rvs_scalpels(rvs_ccf, ccfs, num_basis=n), 1:5)
   println("RMS RVs cleaned by Scalpels: ",std.(rvs_scalpels) )
   dont_need_to!(pipeline_plan,:scalpels)
end

make_plot!(pipeline_plan,:scalpels)
if make_plot(pipeline_plan, :scalpels)
   @assert !need_to(pipeline_plan, :rvs_ccf_total)
   @assert !need_to(pipeline_plan, :ccf_total)
   include("../scripts/plots/scalpels.jl")
   plt = Scalpels.make_plots_scalpels(rvs_ccf, ccfs, max_num_basis=2, v_grid=v_grid, times=order_list_timeseries.times, output_path="examples/output/figures")
   display(plt)
end


if make_plot(pipeline_plan, :rvs_ccf_total)
   using Plots
   rvs_ccf .-= mean(rvs_ccf)
   plt = scatter(rvs_ccf,markersize=3,label="RVs CCF Std", legend=:bottomright)
   if isdefined(Main, :rvs_ccf_expr)
      rvs_ccf_expr .-= mean(rvs_ccf_expr)
      scatter!(plt,rvs_ccf_expr,markersize=3,label="RVs CCF Expr")
   end
   if isdefined(Main, :rvs_scalpels)
      for i in 1:length(rvs_scalpels)
         rvs_scalpels[i] .-= mean(rvs_scalpels[i])
         global plt
         scatter!(plt,rvs_scalpels[i],markersize=3,color=i+2,label="RVs post-Scalpels " * string(i+1))
         plot!(plt,rvs_scalpels[i],markersize=3,color=i+2,label=:none)
      end
   end
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
   #diff = rvs_ccf.-rvs_yale_cbc
   #diff = rvs_ccf2.-rvs_yale_cbc
   #diff = rvs_ccf.-rvs_yale_ccf
   #diff = rvs_ccf2.-rvs_yale_ccf
end
if make_plot(pipeline_plan, :rvs_ccf_total)
   if isdefined(Main, :rvs_ccf_expr)
      diff = rvs_ccf.-rvs_ccf_expr
      println(std(diff))
      plt = scatter(order_list_timeseries.times,diff,markersize=4,label="Delta RV Expr")
   end
   if isdefined(Main, :rvs_scalpels)
      diff = rvs_ccf.-rvs_scalpels[max(2,length(rvs_scalpels))]
      println(std(diff))
      plt = scatter(order_list_timeseries.times,diff,markersize=4,label="Delta RV Scalpels")
   end
   ylabel!("Δv (m/s) (Two mask shapes)")
   xlabel!("Time (d)")
   if save_plot(pipeline_plan,:rvs_ccf_total)   savefig(plt,joinpath(output_dir,target_subdir * "_rvs_ccf_sum.png"))   end
   display(plt)
end


#need_to!(pipeline_plan, :ccf_orders)
(order_ccfs, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_df, pipeline_plan)

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
