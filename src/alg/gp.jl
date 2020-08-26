"""
 gp.jl
Basic code for brute-force Gaussian Process Regression
Adapted from: https://github.com/eford/RvSpectraKitLearn.jl/blob/master/src/deriv_spectra_gp.jl
"""
module GPs
using PDMats

function kernel_matern32(d::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  x = abs(sqrt(3)*d/rho)
  sigmasq * (1+x) * exp(-x)
end

function dkerneldx_matern32(d::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  x = abs(sqrt(3)*d/rho)
  -sigmasq * (-3*d/rho^2) * exp(-x)
end

function d2kerneldx2_matern32(d::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  x = abs(sqrt(3)*d/rho)
  -sigmasq * (3/rho^2) * (1-x)*exp(-x)
end

function kernel_matern52(d::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  x = abs(sqrt(5)*d/rho)
  sigmasq * (1+x*(1+x/3)) * exp(-x)
end

# TODO: Check arithmetic, esp. signs
function dkerneldx_matern52(d::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  x = abs(sqrt(5)*d/rho)
  -sigmasq/3 * (-5*d/rho^2) *(1+x)* exp(-x)
end

# TODO: Check arithmetic, esp. signs
function d2kerneldx2_matern52(d::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  x = abs(sqrt(5)*d/rho)
  -sigmasq/3 * (5/rho^2) * (1-x*(1+x))*exp(-x)
end


function nuttall_kernel(x::T; rho::T = 1.0) where { T<:Real }
   const   a0 = 0.355768
   const   a1 = 0.487396
   const   a2 = 0.144232
   const   a3 = 0.012604
   if abs(x)>rho return zero(T) end
   cpx = cospi(x/rho)
   spx = sinpi(x/rho)
   c2px = cpx*cpx-spx*spx
   return a0+a1*cpx+(a2+a3*cpx)*cpx*cpx-(a2+3*a3*cpx)*spx*spx
end

function nuttall_dkerneldx(x::T; rho::T = 1.0) where { T<:Real }
   const   a0 = 0.355768
   const   a1 = 0.487396
   const   a2 = 0.144232
   const   a3 = 0.012604
   if abs(x)>rho return zero(T) end
   cpx = cospi(x/rho)
   spx = sinpi(x/rho)
   c2px = cpx*cpx-spx*spx
   s2px = 2*cpx*spx
   s3px = spx*c2px+cpx*s2px
   return -pi/rho*(a1*spx+2*a2*s2px+3*a3*s3px)
end

function nuttall_d2kerneldx2(x::T; rho::T = 1.0) where { T<:Real }
   const   a0 = 0.355768
   const   a1 = 0.487396
   const   a2 = 0.144232
   const   a3 = 0.012604
   if abs(x)>rho return zero(T) end
   cpx = cospi(x/rho)
   spx = sinpi(x/rho)
   c2px = cpx*cpx-spx*spx
   c3px = cpx*(cpx*cpx-3*spx*spx)
   return -(pi/rho)^2*(a1*cpx+4*a2*c2px+9*a3*c3px)
end


function matern32_sparse_kernel(x::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  km = kernel_matern32(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*4)
  km*kn
end

function dkerneldx_matern32_sparse(x::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  km = kernel_matern32(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*4)
  dkmdx = dkerneldx_matern32(x,rho=rho,sigmasq=sigmasq)
  dkndx = nuttall_dkerneldx(x,rho=rho*4)
  km*dkndx+kn*dkmdx
end

function d2kerneldx2_matern32_sparse(x::T; rho::T = 1.0, sigmasq::T = 1.0) where { T<:Real }
  km = kernel_matern32(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*4)
  dkmdx = dkerneldx_matern32(x,rho=rho,sigmasq=sigmasq)
  dkndx = nuttall_dkerneldx(x,rho=rho*4)
  d2kmdx2 = d2kerneldx2_matern32(x,rho=rho,sigmasq=sigmasq)
  d2kndx2 = nuttall_d2kerneldx2(x,rho=rho*4)
  km*d2kndx2+2*dkndx*dkmdx+kn*d2kmdx2
end

function make_kernel_data(x::AbstractArray{T,1}, kernel::Function;
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  @assert length(x) == length(sigmasq_obs)
  K = diagm(sigmasq_obs)
  for i in 1:size(K,1)
      for j in 1:size(K,2)
	      K[i,j] += kernel(x[i]-x[j], sigmasq=sigmasq_cor, rho=rho)
	  end
  end
  return PDMat(K)
end

function make_kernel_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1}, kernel::Function;
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  K = zeros(length(xobs),length(xpred))
  for i in 1:length(xobs)
      for j in 1:length(xpred)
	      K[i,j] += kernel(xobs[i]-xpred[j], sigmasq=sigmasq_cor, rho=rho)
	  end
  end
  return K
end

function make_kernel_matern32_data(x::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_data(x, kernel_matern32, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, kernel_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_dmatern32dx_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, dkerneldx_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_d2matern32dx2_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, d2kerneldx2_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_sparse_data(x::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_data(x, matern32_sparse_kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_sparse_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, matern32_sparse_kernel, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_dmatern32dx_sparse_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, dkerneldx_matern32_sparse, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_d2matern32dx2_sparse_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, d2kerneldx2_matern32_sparse, sigmasq_cor=sigmasq_cor, rho=rho)
end

function predict_mean(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  kobs = make_kernel_matern32_sparse_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred = make_kernel_matern32_sparse_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_mean = kobs_pred' * alpha
end

function predict_deriv(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  kobs = make_kernel_matern32_sparse_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred_deriv = make_kernel_dmatern32dx_sparse_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_deriv = kobs_pred_deriv' * alpha
end

function predict_deriv2(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred_deriv2 = make_kernel_d2matern32dx2_sparse_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_deriv = kobs_pred_deriv2' * alpha
end

function predict_mean_and_derivs(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  kobs = make_kernel_matern32_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  kobs_pred = make_kernel_matern32_sparse_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_mean = kobs_pred' * alpha
  kobs_pred_deriv = make_kernel_dmatern32dx_sparse_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv = kobs_pred_deriv' * alpha
  kobs_pred_deriv2 = make_kernel_d2matern32dx2_sparse_obs_pred(xobs,xpred, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv2 = kobs_pred_deriv2' * alpha
  return (pred_mean, pred_deriv, pred_deriv2)
end

function gp_marginal(xobs::AbstractArray{T,1}, yobs::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(xobs)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  kobs = make_kernel_matern32_sparse_data(xobs, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  -0.5*( invquad(kobs, yobs) + logdet(kobs) + length(xobs)*log(2pi) )
end

#=
function calc_gp_marginal_on_segments(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)  where { T<:Real }
  @assert length(lambda) == length(flux) == length(sigmasq_obs)
  println("# sigmasq_obs[1,1] = ", sigmasq_obs[1,1], " sigmasq_cor= ", sigmasq_cor, " rho= ", rho)
  output = 0.0
  num_seg = convert(Int64,ceil(length(lambda)/half_chunck_size)-1)
  for i in 1:num_seg
    idx_begin = 1+half_chunck_size*(i-1)
    idx_end = min(half_chunck_size*(i+1), length(lambda))
    write_idx_begin = idx_begin + div(half_chunck_size,2)
    write_idx_end = idx_end - div(half_chunck_size,2)
    if i==1 write_idx_begin=1 end
    if i==num_seg
      idx_begin = max(1,idx_end-2*half_chunck_size)
      write_idx_end=length(lambda)
    end
    #println("# i= ",i,": ", idx_begin, " - ", idx_end, " -> ", write_idx_begin, " - ", write_idx_end)
    #output[write_idx_begin:write_idx_end] = predict_gp(view(lambda,idx_begin:idx_end), view(flux,idx_begin:idx_end), view(lambda,write_idx_begin:write_idx_end), sigmasq_obs=view(sigmasq_obs,idx_begin:idx_end), sigmasq_cor=sigmasq_cor, rho=rho)
    output += gp_marginal(lambda[idx_begin:idx_end], flux[idx_begin:idx_end], sigmasq_obs=sigmasq_obs[idx_begin:idx_end], sigmasq_cor=sigmasq_cor, rho=rho)
  end
  output
end

function gp_marginal_wrapper(param::Vector)
   @assert length(param) == 2
   idx_min = 40000
   idx_max = 45000
   println("# wrapper: ", param)
   calc_gp_marginal_on_segments(lambda[idx_min:idx_max],vec(mean(obs[idx_min:idx_max],2)), sigmasq_obs = ones(length(lambda[idx_min:idx_max]))/150000,
                    sigmasq_cor = exp(param[1]), rho = exp(param[2]) )
end





function calc_doppler_component_gp(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)  where { T<:Real }
   lambda.*calc_gp_on_segments(predict_deriv,lambda, flux, sigmasq_obs=sigmasq_obs,	sigmasq_cor=sigmasq_cor, rho=rho, half_chunck_size=half_chunck_size)
end

function calc_gp_on_segments(predict_gp::Function, lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)  where { T<:Real }
  @assert length(lambda) == length(flux) == length(sigmasq_obs)
  output = Array{eltype(flux)}(undef,length(lambda))
  num_seg = convert(Int64,ceil(length(lambda)/half_chunck_size)-1)
  for i in 1:num_seg
    idx_begin = 1+half_chunck_size*(i-1)
    idx_end = min(half_chunck_size*(i+1), length(lambda))
    write_idx_begin = idx_begin + div(half_chunck_size,2)
    write_idx_end = idx_end - div(half_chunck_size,2)
    if i==1 write_idx_begin=1 end
    if i==num_seg
      idx_begin = max(1,idx_end-2*half_chunck_size)
      write_idx_end=length(lambda)
    end
    #println("# i= ",i,": ", idx_begin, " - ", idx_end, " -> ", write_idx_begin, " - ", write_idx_end)
    #output[write_idx_begin:write_idx_end] = predict_gp(view(lambda,idx_begin:idx_end), view(flux,idx_begin:idx_end), view(lambda,write_idx_begin:write_idx_end), sigmasq_obs=view(sigmasq_obs,idx_begin:idx_end), sigmasq_cor=sigmasq_cor, rho=rho)
    output[write_idx_begin:write_idx_end] = predict_gp(lambda[idx_begin:idx_end], flux[idx_begin:idx_end], lambda[write_idx_begin:write_idx_end], sigmasq_obs=sigmasq_obs[idx_begin:idx_end], sigmasq_cor=sigmasq_cor, rho=rho)
  end
  output
end

=#

#=
function calc_doppler_component_gp(lambda::AbstractArray{T,1}, flux::AbstractArray{T,2};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  doppler_basis = calc_doppler_component_gp(lambda,vec(mean(flux,2)),sigmasq_obs=sigmasq_obs,sigmasq_cor=sigmasq_cor,rho=rho)
end

function calc_doppler_quadratic_term_gp(lambda::AbstractArray{T,1}, flux::AbstractArray{T,1};
                                sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0,
                                half_chunck_size::Integer = 100)  where { T<:Real }
   0.5*lambda.^2.*calc_gp_on_segments(predict_deriv2,lambda, flux, sigmasq_obs=sigmasq_obs,	sigmasq_cor=sigmasq_cor, rho=rho, half_chunck_size=half_chunck_size)
end

function calc_doppler_quadratic_term_gp(lambda::AbstractArray{T,1}, flux::AbstractArray{T,2};
                                     sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(lambda)),	sigmasq_cor::T = 1.0, rho::T = 1.0)  where { T<:Real }
  doppler_quad_term = calc_doppler_quadratic_term_gp(lambda,vec(mean(flux,2)),sigmasq_obs=sigmasq_obs,sigmasq_cor=sigmasq_cor,rho=rho)
end
=#

#=
function test_gp_deriv()
  #lambda_small = collect(linspace(1000.0,1100.0,101))
  #yobs_small = sin(2pi*lambda/40)+ 0.3*randn(length(lambda_small))
  #lambda_pred = collect(linspace(1000.0,1100.0,201))
  lambda_small = lambda[10000:10100]
  yobs_small = obs[10000:10100,1]
  lambda_pred = lambda_small
  (f,df,df2) = predict_mean_and_derivs(lambda_small,yobs_small,lambda_pred, sigmasq_obs=1.0/150000.0*ones(length(lambda_small)), sigmasq_cor=0.1,rho = 0.01)
  plot(lambda_small,(yobs_small),"r.")
  plot(lambda_pred,f,"b-")
  plot(lambda_pred,df,"g-")
  #plot(lambda_pred,df2,"m-")
end
=#

end
