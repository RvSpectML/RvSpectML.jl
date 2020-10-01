using RvSpectML
 using Statistics
 using Dates

make_plots = false
include("expres_1_read.jl")
order_list_timeseries = RvSpectML.make_order_list_timeseries(expres_data)
order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries,verbose=true)
lambda_range_with_data = (min = maximum(d->minimum(d.λ[expres_data[1].metadata[:excalibur_mask]]),expres_data), max = minimum(d->maximum(d.λ[expres_data[1].metadata[:excalibur_mask]]),expres_data) )
RvSpectML.normalize_spectra!(order_list_timeseries,expres_data);

using StaticArrays
using Stheno, TemporalGPs
using Optim
function optimize_gp_fit_to_data(xobs::AA1, yobs::AA2, xpred::AA3;  # ; kernel::Function = matern52_sparse_kernel,
			sigmasq_obs::AA4 = 1e-16*ones(length(xobs))  #= sigmasq_cor::Real=1.0, rho::Real=1.0)  =# ) where {
				T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T3<:Real, AA3<:AbstractArray{T3,1}, T4<:Real, AA4<:AbstractArray{T4,1} }
	# global ncalls += 1
  	println("# size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
	#=
	kobs = make_kernel_data(xobs, kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  	kobs_pred = make_kernel_obs_pred(xobs, xpred, kernel, sigmasq_cor=sigmasq_cor, rho=rho)
  	alpha = kobs \ yobs
	pred_mean = kobs_pred' * alpha
	=#
	function pack(;σ²=σ²,l=l)
		return [log(σ²-1e-6), log(l-1e-6)]
	end
	function unpack(θ)
	    σ² = exp(θ[1]) + 1e-6
	    l = exp(θ[2]) + 1e-6
	    return σ², l
	end

	function neg_log_marg_likel(θ)
		σ², l = unpack(θ)
	    k = σ² * stretch(Matern52(), 1 / l)
	    f_naive = GP(k, GPC())
		#f = to_sde(f_naive)
		f = to_sde(f_naive, SArrayStorage(Float64))
	    return -logpdf(f(xobs, sigmasq_obs), yobs)
	end

	#logyobs = log.(yobs)
	σ_guess = 0.5*std(yobs)
	length_guess = 0.16
	θ_guess = pack(σ²=σ_guess^2,l=length_guess)

#	return neg_log_marg_likel(θ_guess)

	results = Optim.optimize(neg_log_marg_likel, θ_guess, NelderMead()) #, show_trace=true)
	σ²_ml, l_ml = unpack(results.minimizer)
	return σ²_ml, l_ml

	#=k = σ²_ml * stretch(Matern52(), 1 / l_ml);
	f_naive = GP(k, GPC());
	#f = to_sde(f_naive)
	f = to_sde(f_naive, SArrayStorage(Float64))
	fx = f(xobs, sigmasq_obs./yobs.^2)
	f_posterior_ml = posterior(fx, log.(order_list_timeseries.chunk_list[1][ch].flux) )
	marginals(f_posterior_ml(xobs))
	#plot!(plt, f_posterior_ml(xobs); samples=10, color=:green, label="");
	#display(plt);
	=#
	#=
	@time f = GP( 0.05^2 * stretch(Matern52(), 1 / 0.01) , GPC())
	@time f_tgp = to_sde(f)
	@time fx_tgp = f_tgp(order_list_timeseries.chunk_list[1][ch].λ, order_list_timeseries.chunk_list[1][ch].var./order_list_timeseries.chunk_list[1][ch].flux.^2)
	@time f_post_tgp = posterior(fx_tgp, log.(order_list_timeseries.chunk_list[1][ch].flux) )
	@time f_pred_mean_tgp = marginals(f_post(order_list_timeseries.chunk_list[1][ch].λ))
	@time plot(order_list_timeseries.chunk_list[1][ch].λ,mean.(f_pred_mean_tgp); label="", color=:blue)
	scatter!(order_list_timeseries.chunk_list[1][ch].λ,log.(order_list_timeseries.chunk_list[1][ch].flux),markersize=1.5)
	xlims!(5329,5332)
	=#

end

function calc_marginal_gp(gp_param, xobs::AA1, yobs::AA2, xpred::AA3;  # ; kernel::Function = matern52_sparse_kernel,
				sigmasq_obs::AA4 = 1e-16*ones(length(xobs))  #= sigmasq_cor::Real=1.0, rho::Real=1.0)  =# ) where {
					T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T3<:Real, AA3<:AbstractArray{T3,1}, T4<:Real, AA4<:AbstractArray{T4,1} }
	σ², l = gp_param
	k = σ² * stretch(Matern52(), 1 / l)
	f_naive = GP(k, GPC())
	#f = to_sde(f_naive)
	f = to_sde(f_naive, SArrayStorage(Float64))
	fx = f(xobs, sigmasq_obs)
	f_posterior_ml = posterior(fx, (order_list_timeseries.chunk_list[1][ch].flux) )
	f_marg = marginals(f_posterior_ml(xpred))
	return f_marg
	#plot!(plt, f_posterior_ml(xobs); samples=10, color=:green, label="");
	#display(plt)
end

function neg_log_marg_likel(θ)
	σ², l = θ
	k = σ² * stretch(Matern52(), 1 / l)
	f_naive = GP(k, GPC())
	#f = to_sde(f_naive)
	f = to_sde(f_naive, SArrayStorage(Float64))
	xobs = order_list_timeseries.chunk_list[1][ch].λ
	yobs = (order_list_timeseries.chunk_list[1][ch].flux)
	sigmasq_obs = order_list_timeseries.chunk_list[1][ch].var
	return -logpdf(f(xobs, sigmasq_obs), yobs)
end

ch = 10
res = optimize_gp_fit_to_data(log.(order_list_timeseries.chunk_list[1][ch].λ),
  order_list_timeseries.chunk_list[1][ch].flux,log.(order_list_timeseries.chunk_list[1][ch].λ),sigmasq_obs=
  order_list_timeseries.chunk_list[1][ch].var)
  gp_param = res

@time neg_log_marg_likel(gp_param)
gp_param = [std(order_list_timeseries.chunk_list[1][ch].flux), 0.05]

using Plots
ch = 17
 #gp_param = [.1639394167390819*4, 0.14584100679829712/5615*4]
 gp_param = [ 0.4890909216856761, 5.800274590507981e-5 ]  .* 2
 println("gp_param = ", gp_param, " -log L = ", neg_log_marg_likel(gp_param))
 #xplt = range(log.(minimum(order_list_timeseries.chunk_list[1][ch].λ)),stop=log.(maximum(order_list_timeseries.chunk_list[1][ch].λ)),length=2*length(order_list_timeseries.chunk_list[1][ch].λ))
 xplt = log.(order_list_timeseries.chunk_list[1][ch].λ)
 marg = calc_marginal_gp(gp_param,log.(order_list_timeseries.chunk_list[1][ch].λ),
  order_list_timeseries.chunk_list[1][ch].flux,xplt,sigmasq_obs=
  order_list_timeseries.chunk_list[1][ch].var)
 plot(exp.(xplt), mean.(marg), #= yerr=std.(marg).*mean.(marg) =# )
 scatter!(order_list_timeseries.chunk_list[1][ch].λ, order_list_timeseries.chunk_list[1][ch].flux .- 0.0*mean.(marg), markersize=1.5)
 #xlims!(5295,5300)
 #xlims!(5625,5627)
 #ylims!(0.9,1.01)

nxplt = length(xplt)
 σ², l = gp_param
 k = σ² * stretch(Matern52(), 1 / l)
 #f_naive = GP(k, GPC())
 f_naive = GP(k, GPC())
 #f = to_sde(f_naive)
 f = to_sde(f_naive, SArrayStorage(Float64))
 fx = f(log.(order_list_timeseries.chunk_list[1][ch].λ), order_list_timeseries.chunk_list[1][ch].var./order_list_timeseries.chunk_list[1][ch].flux.^2)
 f_posterior_ml = posterior(fx, log.(order_list_timeseries.chunk_list[1][ch].flux) )

m = (mean.(marginals(f_posterior_ml((xplt)))))
 dmdlnl = zeros(size(m))
 dmdlnl[2:end-1] .= (m[3:end].-m[1:end-2])./(xplt[3]-xplt[1])
 d2mdlnl2[1] = dmdlnl[2]
 d2mdlnl2[end] = dmdlnl[end-1]
 d2mdlnl2 = zeros(size(m))
 d2mdlnl2[2:end-1] .= (m[3:end].+m[1:end-2].-2.0.*m[2:end-1])./(xplt[2]-xplt[1]).^2
 d2mdlnl2[1] = d2mdlnl2[2]
 d2mdlnl2[end] = d2mdlnl2[end-1]
 plot(exp.(xplt),exp.(m),color=:black)
 plot!(exp.(xplt),dmdlnl./std(dmdlnl),color=:blue)
 plot!(exp.(xplt),d2mdlnl2./std(d2mdlnl2),color=:red)
 #idx_line_cand = findall(x->x>0.5,d2mdlnl2./std(d2mdlnl2))
 idx_line_cand = findall(x->x>2e8,d2mdlnl2)
 mean_idx = convert_boundaries_to_means(idx_line_cand)
 line_ctrs = exp.(xplt[mean_idx])
 idx_good_line_ctrs = findall(x->abs(x)<1e3,dmdlnl[convert_boundaries_to_means(idx_line_cand)])
 good_line_ctrs = exp.(xplt[mean_idx[idx_good_line_ctrs]])
 scatter!(line_ctrs,zeros(length(line_ctrs)))
 scatter!(good_line_ctrs,0.2.+zeros(length(good_line_ctrs)))
 scatter!(order_list_timeseries.chunk_list[1][ch].λ,order_list_timeseries.chunk_list[1][ch].flux,markersize=1.5)
ylims!(-0.05,1.25)
xlims!(5635.15,5635.65)

log(5635.65/5635.15)
log(5635.65/5635.15)/log(order_list_timeseries.chunk_list[1][ch].λ[2941]/order_list_timeseries.chunk_list[1][ch].λ[2940])

std(d2mdlnl2)
xlims!(5601,5603.5)
xlims!(5588,5594)
xlims!(5594,5600)
xlims!(5600,5606)
ylims!(0.8,1.25)
xlims!(5640,5644)
 #ylims!(-1,2)
xlims!(5625,5627)
xlims!(5585,5595)
xlims!(5591,5593)
xlims!(5615,5620)
ylims!(0.4,0.65)
xlims!(5635,5645)
xlims!(5635,5636)

histogram(log10.(abs.(dmdlnl[convert_boundaries_to_means(idx_line_cand)])),nbins=60)

findall(x->abs(x)<3e3,dmdlnl[convert_boundaries_to_means(idx_line_cand)])
findall(x->abs(x)<1e3,dmdlnl[convert_boundaries_to_means(idx_line_cand)])

function how_many_lines_found(threshold)
	idx_line_cand = findall(x->x>threshold,d2mdlnl2)
	line_ctrs = exp.(xplt[convert_boundaries_to_means(idx_line_cand)])
	length(line_ctrs)
end

thresholds = [1e10, 3e9, 2e9, 1e9, 5e8, 3e8, 2e8, 1.5e8, 1e8, 8e7, 5e7, 3e7, 2e7, 1e7, 3e6, 1e6]
scatter(log10.(thresholds), how_many_lines_found.(thresholds))

function convert_boundaries_to_means(x::AbstractArray{T,1}) where { T<:Real }
	idx_start = 1
	idx_stop = 1
	centers = Int64[]
	for i in 2:length(x)
		if x[i]-x[i-1] != 1
			idx_stop = i-1
			if idx_start!=1
				push!(centers,floor(Int,(x[idx_start]+x[idx_stop])/2 ))
			end
			idx_start = i
		end
	end
	return centers
end


idx_boundaries = findall(x->x!=0.75,idx_line_cand[2:end].-idx_line_cand[1:end-1])

line_ctrs = map(i->(exp(xplt[floor(Int,idx_boundaries[i-1])])+exp(xplt[floor(Int,idx_boundaries[i])]))/2, 2:length(idx_boundaries))

line_ctrs

exp(xplt[idx_boundaries][end])
exp(maximum(xplt))
size(m)

#=

@time f = GP(0.5^2 * stretch(Matern52(), 1 / 0.1) + 0.5^2 * stretch(Matern52(), 1 / 50) , GPC())
ch = 4
 @time fx = f(order_list_timeseries.chunk_list[1][ch].λ, order_list_timeseries.chunk_list[1][ch].var./order_list_timeseries.chunk_list[1][ch].flux.^2)
 @time f_posterior = f | Obs(fx, log.(order_list_timeseries.chunk_list[1][ch].flux) )
 @time f_poster_on_grid = f_posterior(order_list_timeseries.chunk_list[1][ch].λ)
 @time fpred_mean = marginals(f_poster_on_grid);
 @time plot!(f_poster_on_grid; label="", color=:red)
 scatter!(order_list_timeseries.chunk_list[1][ch].λ,log.(order_list_timeseries.chunk_list[1][ch].flux))

using TemporalGPs
#@time f = GP(0.1^2 * stretch(Matern52(), 1 / 0.1) + 0.1^2 * stretch(Matern52(), 1 / 50) , GPC())
std(log.(order_list_timeseries.chunk_list[1][ch].flux))
@time f = GP( 0.05^2 * stretch(Matern52(), 1 / 0.01) , GPC())
@time f_tgp = to_sde(f)
@time fx_tgp = f_tgp(order_list_timeseries.chunk_list[1][ch].λ, order_list_timeseries.chunk_list[1][ch].var./order_list_timeseries.chunk_list[1][ch].flux.^2)
@time f_post_tgp = posterior(fx_tgp, log.(order_list_timeseries.chunk_list[1][ch].flux) )
@time f_pred_mean_tgp = marginals(f_post(order_list_timeseries.chunk_list[1][ch].λ))
@time plot(order_list_timeseries.chunk_list[1][ch].λ,mean.(f_pred_mean_tgp); label="", color=:blue)
scatter!(order_list_timeseries.chunk_list[1][ch].λ,log.(order_list_timeseries.chunk_list[1][ch].flux),markersize=1.5)
xlims!(5329,5332)

@time logpdf(f_post_tgp(order_list_timeseries.chunk_list[1][ch].λ), log.(order_list_timeseries.chunk_list[1][ch].flux))


=#
