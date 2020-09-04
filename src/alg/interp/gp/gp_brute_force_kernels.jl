"""
   Kernels for GPs

Author: Eric Ford
Adapted from: https://github.com/eford/RvSpectraKitLearn.jl/blob/master/src/deriv_spectra_gp.jl

Could update to use KernelFunctions.jl, etc.
If upgrade to make use of other packages, please keep this so we ahve somet that gets the job done with minimal dependancies.
Author: Eric Ford
"""

function kernel_matern32(d::T; rho::T=one(T), sigmasq::T =one(T)) where { T<:Real }
  x = abs(sqrt(3)*d/rho)
  sigmasq * (1+x) * exp(-x)
end

function dkerneldx_matern32(d::T; rho::T=one(T), sigmasq::T =one(T)) where { T<:Real }
  x = abs(sqrt(3)*d/rho)
  -sigmasq * (3*d/rho^2) * exp(-x)
end

function d2kerneldx2_matern32(d::T; rho::T=one(T), sigmasq::T =one(T)) where { T<:Real }
  x = abs(sqrt(3)*d/rho)
  -sigmasq * (3/rho^2) * (1-x)*exp(-x)
end

function kernel_matern52(d::T; rho::T=one(T), sigmasq::T =one(T)) where { T<:Real }
  x = abs(sqrt(5)*d/rho)
  sigmasq * (1+x*(1+x/3)) * exp(-x)
end

function dkerneldx_matern52(d::T; rho::T=one(T), sigmasq::T =one(T)) where { T<:Real }
  x = abs(sqrt(5)*d/rho)
  -sigmasq/3 * (5*d/rho^2) *(1+x)*exp(-x)
end

function d2kerneldx2_matern52(d::T; rho::T=one(T), sigmasq::T =one(T)) where { T<:Real }
  x = abs(sqrt(5)*d/rho)
  -sigmasq/3 * (5/rho^2) * (1+x*(1-x))*exp(-x)
end

function nuttall_kernel(x::T; rho::T=one(T)) where { T<:Real }
      a0 = 0.355768
      a1 = 0.487396
      a2 = 0.144232
      a3 = 0.012604
   if abs(x)>rho return zero(T) end
   cpx = cospi(x/rho)
   spx = sinpi(x/rho)
   c2px = cpx*cpx-spx*spx
   return a0+a1*cpx+(a2+a3*cpx)*cpx*cpx-(a2+3*a3*cpx)*spx*spx
end

function nuttall_dkerneldx(x::T; rho::T=one(T)) where { T<:Real }
      a0 = 0.355768
      a1 = 0.487396
      a2 = 0.144232
      a3 = 0.012604
   if abs(x)>rho return zero(T) end
   cpx = cospi(x/rho)
   spx = sinpi(x/rho)
   c2px = cpx*cpx-spx*spx
   s2px = 2*cpx*spx
   s3px = spx*c2px+cpx*s2px
   return -pi/rho*(a1*spx+2*a2*s2px+3*a3*s3px)
end

function nuttall_d2kerneldx2(x::T; rho::T=one(T)) where { T<:Real }
      a0 = 0.355768
      a1 = 0.487396
      a2 = 0.144232
      a3 = 0.012604
   if abs(x)>rho return zero(T) end
   cpx = cospi(x/rho)
   spx = sinpi(x/rho)
   c2px = cpx*cpx-spx*spx
   c3px = cpx*(cpx*cpx-3*spx*spx)
   return -(pi/rho)^2*(a1*cpx+4*a2*c2px+9*a3*c3px)
end


function matern32_sparse_kernel(x::T; rho::T=one(T), bandwidth::T=4*one(T), sigmasq::T =one(T)) where { T<:Real }
  km = kernel_matern32(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*bandwidth)
  km*kn
end

function dkerneldx_matern32_sparse(x::T; rho::T=one(T), bandwidth::T=4*one(T), sigmasq::T =one(T)) where { T<:Real }
  km = kernel_matern32(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*bandwidth)
  dkmdx = dkerneldx_matern32(x,rho=rho,sigmasq=sigmasq)
  dkndx = nuttall_dkerneldx(x,rho=rho*bandwidth)
  km*dkndx+kn*dkmdx
end

function d2kerneldx2_matern32_sparse(x::T; rho::T=one(T), bandwidth::T=4*one(T), sigmasq::T =one(T)) where { T<:Real }
  km = kernel_matern32(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*bandwidth)
  dkmdx = dkerneldx_matern32(x,rho=rho,sigmasq=sigmasq)
  dkndx = nuttall_dkerneldx(x,rho=rho*bandwidth)
  d2kmdx2 = d2kerneldx2_matern32(x,rho=rho,sigmasq=sigmasq)
  d2kndx2 = nuttall_d2kerneldx2(x,rho=rho*bandwidth)
  km*d2kndx2+2*dkndx*dkmdx+kn*d2kmdx2
end

function matern52_sparse_kernel(x::T; rho::T=one(T), bandwidth::T=4*one(T), sigmasq::T =one(T)) where { T<:Real }
  km = kernel_matern52(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*bandwidth)
  km*kn
end

function dkerneldx_matern52_sparse(x::T; rho::T=one(T), bandwidth::T=4*one(T), sigmasq::T =one(T)) where { T<:Real }
  km = kernel_matern52(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*bandwidth)
  dkmdx = dkerneldx_matern52(x,rho=rho,sigmasq=sigmasq)
  dkndx = nuttall_dkerneldx(x,rho=rho*bandwidth)
  km*dkndx+kn*dkmdx
end

function d2kerneldx2_matern52_sparse(x::T; rho::T=one(T), bandwidth::T=4*one(T), sigmasq::T =one(T)) where { T<:Real }
  km = kernel_matern52(x,rho=rho,sigmasq=sigmasq)
  kn = nuttall_kernel(x,rho=rho*bandwidth)
  dkmdx = dkerneldx_matern52(x,rho=rho,sigmasq=sigmasq)
  dkndx = nuttall_dkerneldx(x,rho=rho*bandwidth)
  d2kmdx2 = d2kerneldx2_matern52(x,rho=rho,sigmasq=sigmasq)
  d2kndx2 = nuttall_d2kerneldx2(x,rho=rho*bandwidth)
  km*d2kndx2+2*dkndx*dkmdx+kn*d2kmdx2
end

#=
function make_kernel_matern32_data(x::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_data(x, kernel_matern32, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, kernel_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_dmatern32dx_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, dkerneldx_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_d2matern32dx2_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, d2kerneldx2_matern32, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_sparse_data(x::AbstractArray{T,1};
			sigmasq_obs::AbstractArray{T,1} = zeros(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_data(x, matern32_sparse_kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_matern32_sparse_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, matern32_sparse_kernel, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_dmatern32dx_sparse_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, dkerneldx_matern32_sparse, sigmasq_cor=sigmasq_cor, rho=rho)
end

function make_kernel_d2matern32dx2_sparse_obs_pred(xobs::AbstractArray{T,1}, xpred::AbstractArray{T,1};
			#sigmasq_obs::AbstractArray{T,1} = 1e-16*ones(length(x)),
			sigmasq_cor::T = 1.0, rho::T=one(T))  where { T<:Real }
  make_kernel_obs_pred(xobs, xpred, d2kerneldx2_matern32_sparse, sigmasq_cor=sigmasq_cor, rho=rho)
end

=#
