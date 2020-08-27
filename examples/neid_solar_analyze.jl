make_plots_orig = isdefined(Main,:make_plots) ? make_plots : true
make_plots = false
include("neid_solar_2_extract_chunks.jl")
make_plots = make_plots_orig
if make_plots
   using Plots
end

chunk_grids = map(c->RvSpectML.make_grid_for_chunk(chunk_list_timeseries,c), 1:num_chunks(chunk_list_timeseries) )

# Doesn't work because looking up lambda_min and lambda_max, but haven't been updated since splitting order into segments
#chunk_grids = map(c->make_grid_for_chunk(order_list_timeseries,c), 1:num_chunks(order_list_timeseries) )


mean(chunk_grids[1])

# TODO:  RESUME here
#@time (fm, vm, λv, cm) = RvSpectML.pack_chunks_into_matrix(chunk_list_timeseries,chunk_grids)
spectra_matrix = RvSpectML.pack_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_grids)


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
