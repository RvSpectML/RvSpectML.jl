include("neid_solar_read.jl")

using DataFrames, CSV, Query
using Statistics

lambda_range_with_data = (min = maximum(d->minimum(d.λ),solar_data), max = minimum(d->maximum(d.λ),solar_data) )

#vald_filename = "VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesoriginal_depthsoriginal_nbin1depth0.mas"
vald_filename = joinpath(ancilary_data_path,"VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas")
vald_df = RvSpectML.read_mask_vald(vald_filename)
#vald_df = CSV.read(vald_filename,header=["lambda_lo","lambda_hi","depth_vald"])
line_list_df = vald_df |>
  @filter(lambda_range_with_data.min <= _.lambda_lo ) |>
  @filter( _.lambda_hi < lambda_range_with_data.max) |>
#  @filter( _.lambda_lo >6157 || _.lambda_hi < 6155  ) |>   # "Line" w/ large variability
  DataFrame

# TODO: Move adjust width of chunks into read_mask_vald
chunk_size_factor = 3       # TODO: Figure out what value to use
 max_vald_line_offset = 0.0       # km/s
 line_width = RvSpectML.predict_line_width(5780,v_rot=1.8) # # km/s
 Δλoλ_fit_line = (max_vald_line_offset+chunk_size_factor*line_width)*1000/RvSpectML.speed_of_light_mps
 println("# Δλ/λ = ",Δλoλ_fit_line)
 line_list_df.lambda_hi .= line_list_df.lambda*(1 + Δλoλ_fit_line)
 line_list_df.lambda_lo .= line_list_df.lambda/(1 + Δλoλ_fit_line)

chunk_list_df = line_list_df

# Check if chunks overlap
if any(line_list_df.lambda_hi[1:end-1] .>= line_list_df.lambda_lo[2:end])
    idx_overlap = line_list_df.lambda_hi[1:end-1] .>= line_list_df.lambda_lo[2:end]
    println("# Overlapping chunks: ",length(findall(idx_overlap)))
end
minimum(line_list_df.lambda_lo), maximum(line_list_df.lambda_hi)

chunk_list_df = RvSpectML.merge_lines(line_list_df)
size(chunk_list_df)
#chunk_list_df

#solar_data[1].λ[:,1:90]

#@time time_series_of_chunk_lists = map(spec->make_chunk_list(spec,line_list_df),solar_data)
@time time_series_of_chunk_lists = map(spec->RvSpectML.make_chunk_list(spec,chunk_list_df),solar_data)
chunk_list_timeseries = ChunkListTimeseries(df_files_use.bjd,time_series_of_chunk_lists)

chunk_list_timeseries.chunk_list[5].data[10].flux

#chunk_list_timeseries.chunk_list[2].data[10].flux
#solar_data[5].flux ./= 100

# If want to keep orders as one big chunk
#=
neid_orders_to_use = 1:90
pixels_to_use_neid = fill(min_col_neid_default:max_col_neid_default,length(neid_orders_to_use))
if  73<maximum(neid_orders_to_use)   pixels_to_use_neid[73]=1940:max_col_neid_default    end

if true # If want to break orders into 1024 pixel chunks
neid_orders_to_use = 1:60
pixels_to_use_neid = repeat(map(i->560+i*1024:560+(i+1)*1024,0:7),length(neid_orders_to_use))
neid_orders_to_use = mapreduce(i->repeat([i],8),vcat,neid_orders_to_use)
end

=#
solar_data[1]
@time time_series_of_order_lists = map( spec->RvSpectML.make_orders_into_chunks(spec,NEID.NEID2D()), solar_data)
                #orders_to_use=neid_orders_to_use, pixels_to_use=pixels_to_use_neid) ,solar_data)
order_list_timeseries = ChunkListTimeseries(df_files_use.bjd,time_series_of_order_lists)

# Check that no NaN's included
#map(spec->any(map(o->any(isnan.(spec.λ[pixels_to_use_neid[o],o])),1:60)),solar_data)

bad_pixels = map(x->(x[1],x[2]),findall(isnan.(solar_data[5].flux)))
not_in_bad_collum_idx = findall(x->(x[1]<439 || x[1]>450),bad_pixels)
bad_pixels[not_in_bad_collum_idx]

println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )
(chunk_list_timeseries, chunk_list_df) = RvSpectML.filter_bad_chunks(chunk_list_timeseries,chunk_list_df)
println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )

chunk_list_df[1,:]




#normalize_spectra!(chunk_list_timeseries,solar_data)
RvSpectML.normalize_spectra!(order_list_timeseries,solar_data);

solar_data[5].flux[4000:6000]
chunk_list_timeseries.chunk_list[5].data[100].flux

using Plots

#chunk_list_timeseries.chunk_list[1].data[1].λ

order_idx = 20:22
plt = RvSpectML.plot_spectrum_chunks(order_list_timeseries, order_idx)

chunk_idx = 13:22
RvSpectML.plot_spectrum_chunks(chunk_list_timeseries, chunk_idx, plt=plt)

#= Goofing around to see how well chunks line up
chunk_idx = 13:20
xmin = minimum(chunk_list_df.lambda_lo[chunk_idx])
xmax = maximum(chunk_list_df.lambda_hi[chunk_idx])
plt = plot(legend=:none)
#xlims!(xmin,xmax)
for c in chunk_idx
    t = 1
    #if(sum(chunk_list_df.line_depths[c])<0.25) continue end
    λ_mid = sqrt(chunk_list_df.lambda_hi[c]*chunk_list_df.lambda_lo[c])
    println("c= ",c , " λs= ",chunk_list_df.line_λs[c]," depths= ",chunk_list_df.line_depths[c])
    #println("  λlo= ",chunk_list_df.lambda_lo[c]," λhi= ",chunk_list_df.lambda_hi[c], " Δλ= ",chunk_list_df.lambda_hi[c]-chunk_list_df.lambda_lo[c])
    plot!(plt,chunk_list_timeseries.chunk_list[t].data[c].λ.-λ_mid,chunk_list_timeseries.chunk_list[t].data[c].flux)
end
#plot!(plt,solar_data[1].λ,solar_data[1].flux)
#xlims!(4560,4565)
display(plt)
=#


chunk_grids = map(c->RvSpectML.make_grid_for_chunk(chunk_list_timeseries,c), 1:num_chunks(chunk_list_timeseries) )

# Doesn't work because looking up lambda_min and lambda_max, but haven't been updated since splitting order into segments
#chunk_grids = map(c->make_grid_for_chunk(order_list_timeseries,c), 1:num_chunks(order_list_timeseries) )


mean(chunk_grids[1])

# TODO:  RESUME here
@time (fm, vm, λv, cm) = RvSpectML.pack_chunks_into_matrix(chunk_list_timeseries,chunk_grids)

size(fm)
fm_mean = sum(fm./vm,dims=2)./sum(1.0./vm,dims=2)
#fm[cm[1]],λv[cm[1]]

#fm
#fm_mean = mean(fm,dims=2)
deriv_old = copy(deriv)
deriv = calc_mean_deriv(fm,vm,λv,cm)
idx_plt = cm[13]
rvs = vec(sum((fm[idx_plt,:].-fm_mean[idx_plt]).*deriv[idx_plt],dims=1)./sum(abs2.(deriv[idx_plt]))).*speed_of_light_mps
#rvs = vec(sum((fm[idx_plt,:].-fm_mean[idx_plt]).*deriv[idx_plt],dims=1)./sum(abs2.(deriv[idx_plt]))).*speed_of_light_mps
rvs .-= mean(rvs)
println(rvs)
plot(λv[idx_plt] ,fm_mean[idx_plt].-1.0, label="mean")
plot!(λv[idx_plt], deriv[idx_plt]./std(deriv[idx_plt]), label="deriv")
plot!(λv[idx_plt], deriv_old[idx_plt]./std(deriv_old[idx_plt]), label="deriv old")

# Compute RVs
zs = vec(sum((fm.-fm_mean).*deriv,dims=1)./sum(abs2.(deriv)))
rvs = (zs.-mean(zs)).*speed_of_light_mps

σ_rvs = sqrt.(vec(sum(vm.*abs.(deriv),dims=1)./sum(abs2.(deriv))))./length(deriv).*speed_of_light_mps

std(rvs[idx_good])

idx_good = ( (chunk_list_timeseries.times.-minimum(chunk_list_timeseries.times))*hours_per_day .>= 1.0 )
                  # .&  (-3 .<= pca_out[1,:] .<= 3 )

hours_per_day = 24
plot((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,rvs[idx_good],
        yerr=σ_rvs[idx_good],xlabel="Time (hr)",ylabel="RV (m/s)", legend=:none)

obs_per_bin = 5
num_obs_binned = floor(Int,length(findall(idx_good))//obs_per_bin)
idx_binned = map(i->1+(i-1)*obs_per_bin:obs_per_bin*i,1:num_obs_binned)
rvs_binned = map(i->mean(rvs[findall(idx_good)][i]),idx_binned)
times_binned = (map(i->mean(chunk_list_timeseries.times[findall(idx_good)][i]),idx_binned).-minimum(chunk_list_timeseries.times))*hours_per_day
scatter!(times_binned,rvs_binned)

rms_rvs_binned = std(rvs_binned)

fm_perp = fm .- zs'.* deriv
M = fit(PCA, fm_perp[:,idx_good]; maxoutdim=12)
pca_out = MultivariateStats.transform(M,fm_perp[:,idx_good])

scatter(1.0.-cumsum(principalvars(M))./tvar(M), xlabel="Number of PCs", ylabel="Frac Variance Unexplained")

idx_plt = cm[13]
plt0 = plot((fm_mean[idx_plt].-1.0)./std(fm_mean[idx_plt]),legend=:none)
plt0 = plot!(deriv[idx_plt]./std(deriv[idx_plt]),legend=:none)
plt0 = plot!(M.proj[idx_plt,1]./std(M.proj[idx_plt,1]),legend=:none)
plt0 = plot!(M.proj[idx_plt,2]./std(M.proj[idx_plt,2]),legend=:none)
#plt1 = plot(M.proj[idx_plt,1],legend=:none)
plt1 = plot(M.proj[idx_plt,1]./std(M.proj[idx_plt,1]),legend=:none)
plt1 = plot!(M.proj[idx_plt,2]./std(M.proj[idx_plt,2]),legend=:none)
plt2 = plot(M.proj[idx_plt,2],legend=:none)
plt3 = plot(M.proj[idx_plt,3],legend=:none)
plt4 = plot(M.proj[idx_plt,4],legend=:none)
plot(plt0,plt1,plt2,plt3,plt4, layout = (5,1) )

plot((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,rvs[idx_good],
        yerr=σ_rvs[idx_good],xlabel="Time (hr)",ylabel="RV (m/s)", legend=:none)
plot!((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,vec(pca_out[1,:]),
        xlabel="Time (hr)",ylabel="PC1", legend=:none)
plot!((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,vec(pca_out[2,:]),
        xlabel="Time (hr)",ylabel="PC2", legend=:none)

plt1 = scatter(rvs[idx_good],vec(pca_out[1,:]),xlabel="RV",ylabel="PC1",legend=:none)
plt2 = scatter(rvs[idx_good],vec(pca_out[2,:]),xlabel="RV",ylabel="PC2",legend=:none)
plt3 = scatter(rvs[idx_good],vec(pca_out[3,:]),xlabel="RV",ylabel="PC3",legend=:none)
plt4 = scatter(vec(pca_out[1,:]),vec(pca_out[2,:]),xlabel="PC1",ylabel="PC2",legend=:none)
plot(plt1,plt2,plt3,plt4,layout=(2,2))

M2 = fit(PCA, pca_out; maxoutdim=1)

pca_out2 = vec(MultivariateStats.transform(M2, pca_out ))


scatter(rvs[idx_good],pca_out2,xlabel="RV",ylabel="PC1",legend=:none)
plot!(rvs[idx_good],pca_out2,xlabel="RV",ylabel="PC1",legend=:none)



if false
   ts = chunk_list_timeseries
   t = 1
   c = 1
   interp = LinearInterpolation(ts.chunk_list[t].data[c].λ, ts.chunk_list[t].data[c].flux)

    xobs = (ts.chunk_list[t].data[c].λ)
    yobs = copy(ts.chunk_list[t].data[c].flux)
    obs_var = (ts.chunk_list[t].data[c].var)
    med_y = median(yobs)

    # Choose the length-scale and variance of the process.
    l = 0.1
    σ² = 1.0 * med_y^2
    # Construct a kernel with this variance and length scale.
    k = σ² * stretch(Matern52(), 1 / l)

    # Specify a zero-mean GP with this kernel. Don't worry about the GPC object.

    f = GP(med_y, k, GPC())
    fx = f(xobs, obs_var)
    f_posterior = f | Obs(fx, yobs)
    #mean(f_posterior(chunk_grids[c]))

   grid = chunk_grids[c]
   plt1 = plot(grid,mean(f_posterior(grid)))
   plot!(plt1,grid,interp.(grid))
   #plot!(grid,mean(f_posterior(grid)))
   scatter!(plt1,ts.chunk_list[t].data[c].λ, ts.chunk_list[t].data[c].flux)
   plt2 = scatter(ts.chunk_list[t].data[c].λ, mean(f_posterior(ts.chunk_list[t].data[c].λ)).-ts.chunk_list[t].data[c].flux)
   #scatter!(plt2,ts.chunk_list[t].data[c].λ, interp.(ts.chunk_list[t].data[c].λ).-ts.chunk_list[t].data[c].flux)
   plot(plt1,plt2,layout=(2,1)) #display(plot)

end


if false

    map(c->allequal(map(t->length(chunk_list_timeseries.chunk_list[t].data[c].λ),1:length(chunk_list_timeseries))),1:10)

    c=2
    println(map(t->length(chunk_list_timeseries.chunk_list[t].data[c].λ),1:length(chunk_list_timeseries)))
end


if false
    θ_init = [-0.02, 0.1]
xobs = sol_wave[idx_cols,order] .- line_center
yobs = convert(Array{Float64,1},sol_flux[idx_cols,order] )
sigmaobs = convert(Array{Float64,1},sol_var[idx_cols,order] )
(ypred, coeff) = fit_and_predict(θ_init, x=xobs, yobs=yobs, sigmaobs=sigmaobs, degree_poly=1)
println("RMS = ",calc_rms_error(ypred,yobs), " χ^2 = ", calc_chi_sq(ypred,yobs,sigmaobs), " dof = ", length(xobs))
println("coeff = ", coeff)
scatter(xobs,yobs,legend=:none)
plot!(xobs,ypred,legend=:none)
end
