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
 write_ccf_to_csv = true
 write_rvs_to_csv = true
 write_template_to_csv = true
 write_spectral_grid_to_jld2 = true
 write_dcpca_to_csv = true

df_files = make_manifest(target_subdir, EXPRES )
  eval(code_to_include_param_jl())

  if verbose println("# Reading in FITS files.")  end
  @time expres_data = map(EXPRES.read_data,eachrow(df_files_use))

  if verbose println("# Extracting order list timeseries from spectra.")  end
   order_list_timeseries = RvSpectML.make_order_list_timeseries(expres_data)
   order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
   RvSpectML.normalize_spectra!(order_list_timeseries,expres_data);

 if verbose println("# Reading line list for CCF: ", linelist_for_ccf_filename, ".")  end
 lambda_range_with_good_data = get_λ_range(expres_data)
 espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks",linelist_for_ccf_filename)
 espresso_df = RvSpectML.read_linelist_espresso(espresso_filename)
 line_list_df = EXPRES.filter_line_list(espresso_df,first(expres_data).inst)

    # Compute CCF's & measure RVs
if verbose println("# Computing CCF.")  end
mask_shape = RvSpectML.CCF.TopHatCCFMask(order_list_timeseries.inst, scale_factor=4*tophap_ccf_mask_scale_factor)
line_list = RvSpectML.CCF.BasicLineList(line_list_df.lambda, line_list_df.weight)
 ccf_plan = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity)
 v_grid = RvSpectML.CCF.calc_ccf_v_grid(ccf_plan)
 @time ccfs = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
 # Write CCFs to file
  if write_ccf_to_csv
    using CSV
    CSV.write(target_subdir * "_ccfs_old.csv",Tables.table(ccfs',header=Symbol.(v_grid)))
    #CSV.write(target_subdir * "_ccfs_orig.csv",Tables.table(ccfs',header=Symbol.(v_grid)))
  end

mask_shape_new = RvSpectML.CCF.GaussianCCFMask(100000.0,0.01433)
ccf_plan_new = RvSpectML.CCF.BasicCCFPlan(mask_shape = mask_shape_new, line_list=line_list, midpoint=ccf_mid_velocity)
@time ccfs_expr = RvSpectML.CCF.calc_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan_new)
  # Write CCFs to file
   if write_ccf_to_csv
     using CSV
     CSV.write(target_subdir * "_ccfs_new.csv",Tables.table(ccfs_expr',header=Symbol.(v_grid)))
   end

using Plots
ccfs_old_disk = CSV.read(target_subdir * "_ccfs_orig.csv",DataFrame)
ccfs_old_disk = convert(Array{Float64,2},ccfs_old_disk)'

plot(ccfs_old_disk,label=:none)
scatter!(ccfs,markersize=1.5,label=:none)
plot(ccfs_old_disk,markersize=1.5,label=:none)
scatter!(ccfs_expr,markersize=1.5,label=:none)

plot(ccfs,markersize=1.5,label=:none)
scatter!(ccfs_expr,markersize=1.5,label=:none)

1/mean(ccfs_expr./ccfs)

1/mean(ccfs_expr./ccfs_old_disk)
1/mean(ccfs./ccfs_old_disk)

plot(ccfs,label=:none)

plot(ccfs_old_disk.-maximum(ccfs_old_disk,dims=1),label=:none)
scatter!(ccfs_expr.-maximum(ccfs_expr,dims=1),markersize=1.5,label=:none)
scatter!(ccfs.-maximum(ccfs,dims=1),markersize=1.5,label=:none)


mean(maximum(ccfs_old_disk,dims=1)-minimum(ccfs_old_disk,dims=1))
mean(maximum(ccfs_expr,dims=1)-minimum(ccfs_expr,dims=1))

maximum(ccfs_expr,dims=1)

ccfs_expr


scatter(ccfs_expr.-ccfs,markersize=1.5,label=:none)
maximum(abs.(ccfs_expr.-ccfs))
