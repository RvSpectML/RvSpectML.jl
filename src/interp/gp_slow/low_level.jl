include("gp_brute_force_kernels.jl")

function make_kernel_data(x::AA1, kernel::Function; sigmasq_obs::AA2 = zeros(length(x)), sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T1<:Real, AA1<:AbstractArray{T1,1},  T2<:Real, AA2<:AbstractArray{T2,1} }
  @assert length(x) == length(sigmasq_obs)
  K = diagm(sigmasq_obs)
  for i in 1:size(K,1)
      for j in 1:size(K,2)
	      K[i,j] += kernel(x[i]-x[j], sigmasq=sigmasq_cor, rho=rho)
	  end
  end
  return PDMat(K)
end

function make_kernel_obs_pred(xobs::AA1, xpred::AA2,
			kernel::Function; sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T<:Real }
  K = zeros(length(xobs),length(xpred))
  for i in 1:length(xobs)
      for j in 1:length(xpred)
	      K[i,j] += kernel(xobs[i]-xpred[j], sigmasq=sigmasq_cor, rho=rho)
	  end
  end
  return K
end

ncalls = 0
function reset_ncalls()
	global ncalls = 0
end

function predict_mean(xobs::AA1, yobs::AA2, xpred::AA3; kernel::Function = matern52_sparse_kernel,
			sigmasq_obs::AA4 = 1e-16*ones(length(xobs)), sigmasq_cor::Real=1.0, rho::Real=1.0)  where {
				T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}, T3<:Real, AA3<:AbstractArray{T3,1}, T4<:Real, AA4<:AbstractArray{T4,1} }
	# global ncalls += 1
  	println("# size(xobs) = ",size(xobs), "  size(xpred) = ", size(xpred))
	kobs = make_kernel_data(xobs, kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  	kobs_pred = make_kernel_obs_pred(xobs, xpred, kernel, sigmasq_cor=sigmasq_cor, rho=rho)
  	alpha = kobs \ yobs
	pred_mean = kobs_pred' * alpha
end

function predict_deriv(xobs::AA, yobs::AA, xpred::AA;
			kernel::Function = matern52_sparse_kernel, dkerneldx::Function = dkerneldx_matern52_sparse,
			sigmasq_obs::AA = 1e-16*ones(length(xobs)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T<:Real, AA<:AbstractArray{T,1} }
  kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_deriv = kobs_pred_deriv' * alpha
end

function predict_deriv2(xobs::AA, yobs::AA, xpred::AA;
			kernel::Function = matern52_sparse_kernel, d2kerneldx2::Function = d2kerneldx2_matern52_sparse,
			sigmasq_obs::AA = 1e-16*ones(length(xobs)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T<:Real, AA<:AbstractArray{T,1} }
  kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  kobs_pred_deriv2 = make_kernel_obs_pred(xobs,xpred, kernel=d2kerneldx2, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  pred_deriv = kobs_pred_deriv2' * alpha
end

function predict_mean_and_deriv(xobs::AA, yobs::AA, xpred::AA;
			kernel::Function = matern52_sparse_kernel, dkerneldx::Function = dkerneldx_matern52_sparse,
			sigmasq_obs::AA = 1e-16*ones(length(xobs)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T<:Real, AA<:AbstractArray{T,1} }
  kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  kobs_pred = make_kernel_obs_pred(xobs,xpred, kernel=kernel, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_mean = kobs_pred' * alpha
  kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv = kobs_pred_deriv' * alpha
  return (mean=pred_mean, deriv=pred_deriv)
end

function predict_mean_and_derivs(xobs::AA, yobs::AA, xpred::AA;
			kernel::Function = matern52_sparse_kernel, dkerneldx::Function = dkerneldx_matern52_sparse, d2kerneldx2::Function = d2kerneldx2_matern52_sparse,
			sigmasq_obs::AA = 1e-16*ones(length(xobs)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T<:Real, AA<:AbstractArray{T,1} }
  kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  alpha = kobs \ yobs
  kobs_pred = make_kernel_obs_pred(xobs,xpred, kernel=kernel, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_mean = kobs_pred' * alpha
  kobs_pred_deriv = make_kernel_obs_pred(xobs,xpred, kernel=dkerneldx, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv = kobs_pred_deriv' * alpha
  kobs_pred_deriv2 = make_kernel_obs_pred(xobs,xpred, kernel=d2kerneldx2, sigmasq_cor=sigmasq_cor, rho=rho)
  pred_deriv2 = kobs_pred_deriv2' * alpha
  return (mean=pred_mean, deriv=pred_deriv, deriv2=pred_deriv2)
end

function gp_marginal(xobs::AA, yobs::AA, kernel::Function;
			sigmasq_obs::AA = 1e-16*ones(length(xobs)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T<:Real, AA<:AbstractArray{T,1} }
  kobs = make_kernel_data(xobs, kernel=kernel, sigmasq_obs=sigmasq_obs, sigmasq_cor=sigmasq_cor, rho=rho)
  -0.5*( invquad(kobs, yobs) + logdet(kobs) + length(xobs)*log(2pi) )
end




#=
function calc_gp_marginal_on_segments(lambda::AA, flux::AA;
                                sigmasq_obs::AA = 1e-16*ones(length(lambda)),	sigmasq_cor::Real=1.0, rho::Real=1.0,
                                half_chunk_size::Integer = 100)  where { T<:Real, AA<:AbstractArray{T,1} }
  @assert length(lambda) == length(flux) == length(sigmasq_obs)
  println("# sigmasq_obs[1,1] = ", sigmasq_obs[1,1], " sigmasq_cor= ", sigmasq_cor, " rho= ", rho)
  output = 0.0
  num_seg = convert(Int64,ceil(length(lambda)/half_chunk_size)-1)
  for i in 1:num_seg
    idx_begin = 1+half_chunk_size*(i-1)
    idx_end = min(half_chunk_size*(i+1), length(lambda))
    write_idx_begin = idx_begin + div(half_chunk_size,2)
    write_idx_end = idx_end - div(half_chunk_size,2)
    if i==1 write_idx_begin=1 end
    if i==num_seg
      idx_begin = max(1,idx_end-2*half_chunk_size)
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





function calc_doppler_component_gp(lambda::AA, flux::AA;
                                sigmasq_obs::AA = 1e-16*ones(length(lambda)),	sigmasq_cor::Real=1.0, rho::Real=1.0,
                                half_chunk_size::Integer = 100)  where { T<:Real, AA<:AbstractArray{T,1} }
   lambda.*calc_gp_on_segments(predict_deriv,lambda, flux, sigmasq_obs=sigmasq_obs,	sigmasq_cor=sigmasq_cor, rho=rho, half_chunk_size=half_chunk_size)
end

function calc_gp_on_segments(predict_gp::Function, lambda::AA, flux::AA;
                                sigmasq_obs::AA = 1e-16*ones(length(lambda)),	sigmasq_cor::Real=1.0, rho::Real=1.0,
                                half_chunk_size::Integer = 100)  where { T<:Real, AA<:AbstractArray{T,1} }
  @assert length(lambda) == length(flux) == length(sigmasq_obs)
  output = Array{eltype(flux)}(undef,length(lambda))
  num_seg = convert(Int64,ceil(length(lambda)/half_chunk_size)-1)
  for i in 1:num_seg
    idx_begin = 1+half_chunk_size*(i-1)
    idx_end = min(half_chunk_size*(i+1), length(lambda))
    write_idx_begin = idx_begin + div(half_chunk_size,2)
    write_idx_end = idx_end - div(half_chunk_size,2)
    if i==1 write_idx_begin=1 end
    if i==num_seg
      idx_begin = max(1,idx_end-2*half_chunk_size)
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
function calc_doppler_component_gp(lambda::AA, flux::AbstractArray{T,2};
                                sigmasq_obs::AA = 1e-16*ones(length(lambda)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where { T<:Real, AA<:AbstractArray{T,1} }
  doppler_basis = calc_doppler_component_gp(lambda,vec(mean(flux,2)),sigmasq_obs=sigmasq_obs,sigmasq_cor=sigmasq_cor,rho=rho)
end

function calc_doppler_quadratic_term_gp(lambda::AA, flux::AA;
                                sigmasq_obs::AA = 1e-16*ones(length(lambda)),	sigmasq_cor::Real=1.0, rho::Real=1.0,
                                half_chunk_size::Integer = 100)  where { T<:Real, AA<:AbstractArray{T,1} }
   0.5*lambda.^2.*calc_gp_on_segments(predict_deriv2,lambda, flux, sigmasq_obs=sigmasq_obs,	sigmasq_cor=sigmasq_cor, rho=rho, half_chunk_size=half_chunk_size)
end

function calc_doppler_quadratic_term_gp(lambda::AA, flux::AbstractArray{T,2};
                                     sigmasq_obs::AA = 1e-16*ones(length(lambda)),	sigmasq_cor::Real=1.0, rho::Real=1.0)  where rho*bandwidth{ T<:Real, AA<:AbstractArray{T,1} }
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
