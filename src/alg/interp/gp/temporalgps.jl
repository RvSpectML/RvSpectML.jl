"""
Author: Eric Ford
Adapted from: https://github.com/eford/RvSpectraKitLearn.jl/blob/master/src/deriv_spectra_gp.jl

GPs using Stheno.jl, TemporalGPs., KernelFunctions.jl, etc.
"""

"""
Module for interpolating via Gaussian Process Regression based on Stheno and TemporalGPs packages.
"""
module TemporalGPInterpolation
using LinearAlgebra
using PDMats
using StaticArrays
using Stheno, TemporalGPs
using Dates
using ..RvSpectML

#export make_kernel_data, make_kernel_obs_pred,
export gp_marginal
export predict_mean, predict_deriv, predict_deriv2, predict_mean_and_deriv, predict_mean_and_derivs

import Stheno: AbstractGP
import Distributions: AbstractMvNormal

using Distributions
const Fx_PosteriorType = Distribution{Multivariate,Continuous}

function construct_gp(; smooth_factor::Real = 1 ) #xobs::AA1, yobs::AA2, xpred::AA3; #= kernel::Function = matern52_sparse_kernel, =# sigmasq_obs::AA4 = 1e-16*ones(length(xobs)) #=, sigmasq_cor::Real=1.0, rho::Real=1.0)  =#
			#, use_logx::Bool = true, use_logy::Bool = true )   where {
				#T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T3<:Real, AA3<:AbstractArray{T3,1}, T4<:Real, AA4<:AbstractArray{T4,1} }
	#gp_param_default = [.1639394167390819, 0.01/5615] .* smooth_factor   # TODO Generalize.  Values fr fit to one order of one EXPRES spectra
	#gp_param_default = [.1639394167390819, 0.01/5615 .* smooth_factor  ] # TODO Generalize.  Values fr fit to one order of one EXPRES spectra
	#gp_param_default = [.1639394167390819, 0.14584100679829712/5615] .* smooth_factor   # TODO Generalize.  Values fr fit to one order of one EXPRES spectra
	#println("smooth_factor = ",smooth_factor)
	gp_param_default = [ 0.4890909216856761, 5.800274590507981e-5] .* smooth_factor   # TODO Generalize.  Values fr fit to one order of one EXPRES spectra
	#gp_param_default = [ 0.4890909216856761, 5.800274590507981e-5 * smooth_factor ]
	σ², l = gp_param_default
	k = σ² * stretch(Matern52(), 1 / l)
	f_naive = GP(k, GPC())
	#f = to_sde(f_naive)   # if develop issues with StaticArrays could revert to this
	f = to_sde(f_naive, SArrayStorage(Float64))
end

function construct_gp_posterior(xobs::AA1, yobs::AA2, xpred::AA3; #= kernel::Function = matern52_sparse_kernel, =# sigmasq_obs::AA4 = 1e-16*ones(length(xobs)) #=, sigmasq_cor::Real=1.0, rho::Real=1.0)  =#
			, use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1, boost_factor::Real = 1 )   where {
				T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T3<:Real, AA3<:AbstractArray{T3,1}, T4<:Real, AA4<:AbstractArray{T4,1} }
	f = construct_gp( smooth_factor=smooth_factor ) # xobs,yobs,xpred,sigmasq_obs=sigmasq_obs)
	xobs_trans = use_logx ? log.(xobs./boost_factor) : xobs./boost_factor
    yobs_trans = use_logy ? log.(yobs) : yobs
    sigmasq_obs_trans = use_logy ? sigmasq_obs./yobs.^2 : sigmasq_obs
	#=
	if xobs == xpred
		xpred_trans = xobs_trans
	else
    	xpred_trans = use_logx ? log.(xpred) : xpred
	end
	=#
	fx = f(xobs_trans, sigmasq_obs_trans)
	f_posterior = posterior(fx, yobs_trans )
end

function predict_gp(gp::AGP,  xpred::AA3, use_logx::Bool = true, use_logy::Bool = true )   where {
		AGP<:AbstractGP, T3<:Real, AA3<:AbstractArray{T3,1} }
		m = mean.(gp(xpred))
		s = sqrt.(var.(gp(xpred)))
		return m, s
end

function predict_mean(gp::AGP ; use_logy::Bool = true )   where { AGP<:Fx_PosteriorType }
	#xpred_trans = use_logx ? log.(xpred) : xpred
	gp_marginals = marginals(gp) # (xpred_trans))
	output = use_logy ? exp.(mean.(gp_marginals)) : mean.(gp_marginals)
end

function predict_deriv(gp::AGP, xpred::AA3; use_logx::Bool = true, use_logy::Bool = true )   where { AGP<:Fx_PosteriorType,
				T3<:Real, AA3<:AbstractArray{T3,1}  }
  #kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  #kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  #alpha = kobs \ yobs
  #pred_deriv = kobs_pred_deriv' * alpha
  xpred_trans = use_logx ? log.(xpred) : xpred
  m = predict_mean(gp, use_logy=use_logy)
  #if use_logy   m .= exp(m)   end
  dfluxdlnλ = zeros(size(m))
  dfluxdlnλ[1] = (m[2]-m[1])/(xpred_trans[2]-xpred_trans[1])
  dfluxdlnλ[2:end-1] .= (m[3:end].-m[1:end-2])./(xpred_trans[3:end].-xpred_trans[1:end-2]) # exp.(m[2:end-1]).*
  dfluxdlnλ[end] = (m[end]-m[end-1])/(xpred_trans[end]-xpred_trans[end-1])
  return dfluxdlnλ

end

function predict_deriv2(gp::AGP, xpred::AA3; use_logx::Bool = true, use_logy::Bool = true )   where { AGP<:Fx_PosteriorType,
				T3<:Real, AA3<:AbstractArray{T3,1}  }
	#kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
	#kobs_pred_deriv2 = make_kernel_obs_pred(xobs,xpred, kernel=d2kerneldx2, sigmasq_cor=sigmasq_cor, rho=rho)
	#alpha = kobs \ yobs
	#pred_deriv = kobs_pred_deriv2' * alpha

	xpred_trans = use_logx ? log.(xpred) : xpred
	m = predict_mean(gp,use_logy=use_logy)
	#if use_logy   m .= exp(m)   end
	d2fluxdlnλ2 = zeros(size(m))
	d2fluxdlnλ2[2:end-1] .= (m[3:end].+m[1:end-2].-2.0.*m[2:end-1])./(xpred_trans[3:end].-xpred_trans[1:end-2]).^2
	d2fluxdlnλ2[1] = d2fluxdlnλ2[2]
	d2fluxdlnλ2[end] = d2fluxdlnλ2[end-1]

	# TODO add some sort of check based on second derivative and dx that the finite difference is good enough?
	return d2fluxdlnλ2
end

function predict_mean(xobs::AA1, yobs::AA2, xpred::AA3;	sigmasq_obs::AA4 = 1e-16*ones(length(xobs)) #=, sigmasq_cor::Real=1.0, rho::Real=1.0)  =#
			, use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1, boost_factor::Real = 1  )   where {
				T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T3<:Real, AA3<:AbstractArray{T3,1}, T4<:Real, AA4<:AbstractArray{T4,1} }
	# global ncalls += 1
	@assert size(xobs) == size(yobs) == size(sigmasq_obs)
  	#println("# predict_mean (TemporalGPs): size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
	tstart = now()
	f_posterior = construct_gp_posterior(xobs,yobs,xpred,sigmasq_obs=sigmasq_obs, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor, boost_factor=boost_factor )
	#println("typeof(f_posterior) = ",typeof(f_posterior))
	#println("f_posterior <: AbstractGP = ",typeof(f_posterior) <: AbstractGP )
	#println("f_posterior <: AbstractMvNormal = ",typeof(f_posterior) <: AbstractMvNormal )
	xpred_trans = use_logx ? log.(xpred) : xpred
	fx_posterior = f_posterior(xpred_trans)
	#println("typeof(fx_posterior) = ",typeof(fx_posterior))
	#println("f_posterior <: AbstractGP = ",typeof(fx_posterior) <: AbstractGP )
	#println("f_posterior <: AbstractMvNormal = ",typeof(fx_posterior) <: AbstractMvNormal )
	#println("f_posterior <: Fx_PosteriorType = ",typeof(fx_posterior) <: Fx_PosteriorType )
	#output = predict_mean(f_posterior(xpred_trans), xpred_trans ) #, use_logx=use_logx,use_logy=use_logy)
	output = predict_mean(f_posterior(xpred_trans), use_logy=use_logy)
	#println("# predict_mean (TemporalGPs) runtime: ", now()-tstart)
	return output
end

# NEED TO TEST
function predict_deriv(xobs::AA, yobs::AA, xpred::AA; sigmasq_obs::AA = 1e-16*ones(length(xobs))
						, use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1  )   where { T<:Real, AA<:AbstractArray{T,1} }
#			kernel::Function = matern52_sparse_kernel, dkerneldx::Function = dkerneldx_matern52_sparse,
#			sigmasq_cor::Real=1.0, rho::Real=1.0)
  #kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  #kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  #alpha = kobs \ yobs
  #pred_deriv = kobs_pred_deriv' * alpha
  println("# predict_deriv (TemporalGPs): size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
  xpred_trans = use_logx ? log.(xpred) : xpred
  tstart = now()
  f_posterior = construct_gp_posterior(xobs,yobs,xpred,sigmasq_obs=sigmasq_obs, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor)
  output = predict_deriv(f_posterior(xpred_trans), vec(xpred_trans), use_logx=use_logx,use_logy=use_logy)
  println("# predict_deriv (TemporalGPs) runtime: ", now()-tstart)
	return output
end

# NEED TO TEST
function predict_deriv2(xobs::AA, yobs::AA, xpred::AA;sigmasq_obs::AA = 1e-16*ones(length(xobs))
						, use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1  )  where { T<:Real, AA<:AbstractArray{T,1} }
#			kernel::Function = matern52_sparse_kernel, d2kerneldx2::Function = d2kerneldx2_matern52_sparse,
#			,	sigmasq_cor::Real=1.0, rho::Real=1.0
  #kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  #kobs_pred_deriv2 = make_kernel_obs_pred(xobs,xpred, kernel=d2kerneldx2, sigmasq_cor=sigmasq_cor, rho=rho)
  #alpha = kobs \ yobs
  #pred_deriv = kobs_pred_deriv2' * alpha

  println("# predict_deriv2 (TemporalGPs): size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
  tstart = now()
  xpred_trans = use_logx ? log.(xpred) : xpred
  f_posterior = construct_gp_posterior(xobs,yobs,xpred,sigmasq_obs=sigmasq_obs, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
  output = predict_deriv2(f_posterior(xpred_trans), xpred_trans, use_logx=use_logx,use_logy=use_logy)
  println("# predict_deriv2 (TemporalGPs) runtime: ", now()-tstart)
    return output
end

function predict_mean_and_deriv(xobs::AA, yobs::AA, xpred::AA;sigmasq_obs::AA = 1e-16*ones(length(xobs))
	 							, use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1 )   where { T<:Real, AA<:AbstractArray{T,1} }
	#		kernel::Function = matern52_sparse_kernel, dkerneldx::Function = dkerneldx_matern52_sparse,
	#		sigmasq_cor::Real=1.0, rho::Real=1.0)
  #kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  #alpha = kobs \ yobs
  #kobs_pred = make_kernel_obs_pred(xobs,xpred, kernel=kernel, sigmasq_cor=sigmasq_cor, rho=rho)
  #pred_mean = kobs_pred' * alpha
  #kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  #pred_deriv = kobs_pred_deriv' * alpha
  println("# predict_mean_and_deriv (TemporalGPs): size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
  tstart = now()
  f_posterior = construct_gp_posterior(xobs,yobs,xpred,sigmasq_obs=sigmasq_obs, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
	# TODO Opt:  Avoid repeated calculating of mean
  pred_mean = predict_mean(f_posterior(xpred_trans), use_logx=use_logx,use_logy=use_logy)
  pred_deriv = predict_deriv(f_posterior(xpred_trans), use_logx=use_logx,use_logy=use_logy)
  println("# predict_mean_and_deriv (TemporalGPs) runtime: ", now()-tstart)

  return (mean=pred_mean, deriv=pred_deriv)
end

function predict_mean_and_derivs(xobs::AA, yobs::AA, xpred::AA; sigmasq_obs::AA = 1e-16*ones(length(xobs))
								, use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1  ) where { T<:Real, AA<:AbstractArray{T,1} }
	#		kernel::Function = matern52_sparse_kernel, dkerneldx::Function = dkerneldx_matern52_sparse, d2kerneldx2::Function = d2kerneldx2_matern52_sparse,
	#		,	sigmasq_cor::Real=1.0, rho::Real=1.0)
  #=
  kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  kobs_pred = make_kernel_obs_pred(xobs,xpred, kernel=kernel, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_mean = kobs_pred' * alpha
  kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv = kobs_pred_deriv' * alpha
  kobs_pred_deriv2 = make_kernel_obs_pred(xobs,xpred, kernel=d2kerneldx2, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv2 = kobs_pred_deriv2' * alpha
  =#
  #println("# predict_mean_and_derivs (TemporalGPs): size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
  xobs_trans = use_logx ? log.(xobs) : xobs
  yobs_trans = use_logy ? log.(yobs) : yobs
  xpred_trans = use_logx ? log.(xpred) : xpred
  sigmasq_obs_trans = use_logy ? sigmasq_obs./yobs.^2 : sigmasq_obs

  tstart = now()
  f_posterior = construct_gp_posterior(xobs_trans,yobs_trans,xpred_trans,sigmasq_obs=sigmasq_obs_trans, use_logx=false, use_logy=false, smooth_factor=smooth_factor )
  pred_mean = predict_mean(f_posterior(xpred_trans), use_logy=false)  # doesn't need use_logx
  if use_logy
	pred_mean = exp.(pred_mean)
  end
  pred_deriv = calc_dfluxdlnlambda(pred_mean,xpred_trans)
  pred_deriv2 = calc_d2fluxdlnlambda2(pred_mean,xpred_trans)
  #pred_deriv = predict_deriv(f_posterior, xpred_trans, use_logx=false,use_logy=false)
  #pred_deriv2 = predict_deriv2(f_posterior, xpred_trans, use_logx=false,use_logy=false)
  #println("# predict_mean_and_derivs (TemporalGPs) runtime: ", now()-tstart)
  return (mean=pred_mean, deriv=pred_deriv, deriv2=pred_deriv2)
end

function gp_marginal(xobs::AA, yobs::AA #, kernel::Function;
			; sigmasq_obs::AA = 1e-16*ones(length(xobs)),
			use_logx::Bool = true, use_logy::Bool = true, smooth_factor::Real = 1  )  where { T<:Real, AA<:AbstractArray{T,1} }
  	#kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  	# -0.5*( invquad(kobs, yobs) + logdet(kobs) + length(xobs)*log(2pi) )
	xobs_trans = use_logx ? log.(xobs) : xobs
    yobs_trans = use_logy ? log.(yobs) : yobs
	# TODO Opt check if logs are being recalculated twice or otpimized away by compiler
    sigmasq_obs_trans = use_logy ? sigmasq_obs./yobs.^2 : sigmasq_obs
	f_posterior = construct_gp_posterior(xobs,yobs,xobs,sigmasq_obs=sigmasq_obs, use_logx=use_logx, use_logy=use_logy, smooth_factor=smooth_factor )
	return -logpdf(f_posterior(xobs_trans, sigmasq_obs_trans), yobs_trans)

end

end # module
