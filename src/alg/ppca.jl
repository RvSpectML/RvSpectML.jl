"""
Module for performing a Probilistic PCA analysis

Author: Eric Ford
"""
module PPCA

using LinearAlgebra: dot, Symmetric
using Statistics: mean, std
using Flux                           # For optimization
using PDMats, BandedMatrices         # For regularizers
using Stheno #, TemporalGPs

abstract type AbstractPPCA end

struct PPCAGridded <: AbstractPPCA
  W::Array{Float64,2}
  b::Array{Float64,1}
  x::Array{Float64,2}
  lambda::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}

  function PPCAGridded(W::AA1, b::AA2, x::AA3, lambda_range::R) where { T<:Real, AA1<:AbstractArray{T,2}, AA2<:AbstractArray{T,1}, AA3<:AbstractArray{T,2}, R<:AbstractRange{T} }
     @assert size(W,1) == size(b,1)
     @assert size(W,2) == size(x,1)
     @assert length(lambda_range) == size(b,1)
     new(W, b, x, lambda_range )
  end

  function PPCAGridded(in::PPCAGridded)
     @assert size(in.W,1) == size(in.b,1)
     @assert size(in.W,2) == size(in.x,1)
     @assert length(in.lambda) == size(in.b,1)
    new(in.W,in.b,in.x,in.lambda)
  end

end

import Base.copy
copy(in::PPCAGridded) = PPCAGridded(in)

function PPCAGridded(in::Integer, out::Integer, nobs::Integer, lambda_range::R) where {R<:AbstractRange}
    @assert length(lambda_range) == out
    PPCAGridded(randn(out, in), randn(out), randn(in,nobs), lambda_range)
end

get_num_inputs(m::PPCAGridded) = size(m.W,2)
get_num_outputs(m::PPCAGridded) = size(m.W,1)
get_num_param(m::PPCAGridded) = get_num_outputs(m)*(get_num_inputs(m)+1)

Flux.@functor PPCAGridded (W,b,x)

# Overload call, so the object can be used as a function
function (m::PPCAGridded)(x::AbstractArray{T})::Array{T,2} where { T<:Real }
    m.W * m.x .+ m.b
end

struct SpectraSet1D
    lambda_obs::Array{Float64,2}
    lambda_ssb::Array{Float64,2}
    flux::Array{Float64,2}
    flux_var::Array{Float64,2}
end

get_num_outputs(s::SpectraSet1D) = size(s.flux,1)
get_num_obs(s::SpectraSet1D) = size(s.flux,2)

function chisq(x::AbstractArray{T,2}, m::PPCAGridded, ss::SpectraSet1D)::T where { T<:Real }
    yhat = m(x) # our predictions for fluxes's
    sum(((ss.flux.-yhat).^2 ./ss.flux_var))
end

chisq(m::PPCAGridded, ss::SpectraSet1D) = chisq(m.x,m,ss)

mutable struct Regularizer
    covar::AbstractPDMat{Float64}
    scale_input::Float64
    scale_output::Float64
    num_bands::Int64
    kernel::Function
end

const sqrt5 = sqrt(5)
matern52(x) = (one(x)+x*(sqrt5+x*sqrt5/3))*exp(-sqrt5*x)
expsq(x) = exp(-0.5*x^2)

function make_covar_banded(n::Integer, b::Integer, kernel::Function; lambda::Real = 1.0)
    @assert b == 5  # TODO: Generalize for arbitrary number of bands
    lags = 1:b
    covar_at_lag = kernel.(lags./(lambda*b))
    PDMat(Symmetric(BandedMatrix((0=>ones(n),1=>covar_at_lag[1]*ones(n-1),2=>covar_at_lag[2]*ones(n-2),3=>covar_at_lag[3]*ones(n-3),4=>covar_at_lag[4]*ones(n-4),5=>covar_at_lag[5]*ones(n-5)), (n,n), (b,b) )))
end

function update_lambda!(r::Regularizer, lambda::Real)
    r.scale_input = lambda
    r.covar = make_covar_banded(size(r.covar,1), r.num_bands, r.kernel, lambda=lambda)
    return r
end

function update_sigma!(r::Regularizer, sigma::Real)
    r.scale_output = sigma
    return r
end

function (r::Regularizer)(x::AbstractArray{T,1})::T where { T<:Real }
    r.scale_output*invquad(r.covar,x./r.scale_input)
end

function (r::Regularizer)(x::AbstractArray{T,2})::T where { T<:Real }
    sum(r(x[:,i]) for i in 1:size(x,2) )
end

function Matern52Regularizer(nobs::Integer, nband::Integer; lambda::Float64 = 1.0, sigma::Float64 = 1.0)
    b = make_covar_banded(nobs,nband, matern52, lambda=lambda)
    return Regularizer(b, lambda, sigma, nband, matern52)
end

function ExpSqRegularizer(nobs::Integer, nband::Integer; lambda::Float64 = 1.0, sigma::Float64 = 1.0)
    b = make_covar_banded(nobs,nband, expsq, lambda=lambda)
    return Regularizer(b, lambda, sigma, nband, expsq)
end

function regularization_L2norm(x::AbstractArray{T}) where {T<:Real}
    sum(abs2.(x))
end

function regularization_L1norm(x::AbstractArray{T}) where {T<:Real}
    sum(abs.(x))
end

function regularization_deriv(x::AbstractArray{T,2})  where { T<:Real }
    sum(abs2.(x[2:end].-x[1:end-1]))
end

function regularization_deriv2(x::AbstractArray{T,2})  where { T<:Real }
    sum(abs2.(-x[1:end-2].+2*x[2:end-1].-x[3:end]))
end

function regularization_scores(x::AbstractArray{T,2})  where { T<:Real }
    regularization_L2norm(x)
end

function regularization_mean(x::AbstractArray{T,1}; mu=zero(T) )  where { T<:Real }
    (gp_reg.scale_output == zero(T)) ? zero(T) : gp_reg(x.-mu)
end

function regularization_basis(x::AbstractArray{T,2})  where { T<:Real }
    (gp_reg.scale_output == zero(T)) ? zero(T) : gp_reg(x)
end

function regularization(x::AbstractArray{T,2}, m::PPCAGridded) #=, ss::SpectraSet1D)=#  where { T<:Real }
    regularization_scores(x) + regularization_basis(m.W) + regularization_mean(m.b, mu=1.0)
end

function regularization(m::PPCAGridded) #=, ss::SpectraSet1D) =# where { T<:Real }
    regularization_scores(m.x) + regularization_basis(m.W) + regularization_mean(m.b, mu=1.0)
end

function loss(x::AbstractArray{T,2}, m::M, ss::SpectraSet1D)  where { T<:Real, M<:AbstractPPCA }
    chisq(x,m,ss) + regularization(x,m)
end

function loss(m::M, ss::SpectraSet1D)  where { M<:AbstractPPCA }
    chisq(m.x,m,ss) + regularization(m.x,m)
end

# random initial parameters
num_outputs = 50
num_inputs = 3
num_obs = 50
snr_per_pixel = 50
snr_per_pixel = 200

lambda_min = 5000
lambda_max = 5025
#lambda_pad = (lambda_max-lambda_min)/num_outputs
#lambda_range_extra = range(lambda_min-lambda_pad,stop=lambda_max+lambda_pad,length=num_outputs+2)
lambda_range = range(lambda_min,stop=lambda_max,length=num_outputs)
# ppca = PPCAGridded(num_inputs,num_outputs, lambda_range)

# Create a PPCA model to create "true" data
ppca_mean = 1.01 .- 0.01*4.0.*((lambda_range.-mean(lambda_range))./(maximum(lambda_range)-minimum(lambda_range))).^2
ppca_basis = zeros(num_outputs,num_inputs+1)
ppca_basis[:,1] .= sin.(2*pi.*lambda_range/20.0)
ppca_basis[:,2] .= sin.(2*pi.*lambda_range/10.0)
ppca_basis[:,3] .= sin.(2*pi.*lambda_range/4.0)
ppca_basis[:,4] .= 0.0; # sin.(2*pi.*lambda_range/6.0)

x_true = randn(num_inputs+1,num_obs)        # Generate random weights for each basis vector at each epoch
ppca_true = PPCAGridded(ppca_basis, ppca_mean, x_true, lambda_range)
ytrue = ppca_true(x_true)
#ytrue = ppca_basis*x_true.+ppca_mean
yobs = ytrue .+ 1.0/snr_per_pixel * randn(num_outputs,num_obs)  # Add measurement noise
x = copy(x_true[1:num_inputs,1:num_obs])

lambda_ranges = reshape(repeat(collect(lambda_range),num_obs),num_outputs,num_obs)
doppler_factor = 1.001
spec = SpectraSet1D(lambda_ranges.*doppler_factor, lambda_ranges, yobs, ones(size(yobs))./snr_per_pixel )

@time chisq(ppca_true,spec)

using Plots

plot(ppca_true.lambda,ppca_true.W[:,1], color=:blue, label="W1", title="Basis vectors (true)")
#scatter!(ppca_true.lambda,ppca.W[:,1]./maximum(ppca.W[:,1]), ms=2, color=:blue)
plot!(ppca_true.lambda,ppca_true.W[:,2], color=:red, label="W2")
#scatter!(ppca_true.lambda,ppca.W[:,2]./maximum(ppca.W[:,2]), ms=2,color=:red)
plot!(ppca_true.lambda,ppca_true.W[:,3], color=:green, label="W3")
plot!(ppca_true.lambda,ppca_true.W[:,4], color=:cyan, label="W4")
#scatter!(ppca_true.lambda,ppca.W[:,3]./maximum(ppca.W[:,3]), ms=2, color=:green)
plot!(ppca_true.lambda,ppca_true.b, color=:black, label="b")

plot(lambda_range,ppca_true(x_true)[:,1:5],lw=2, label=:none, title="Sample Obs")
#scatter!(lambda_range,ppca(x_true)[:,1], label="Data (init)", color=:lightblue)
xlabel!("λ") # * L"\A")

# Choose Regularization method
gp_reg = Matern52Regularizer(num_outputs, 5, lambda=0.4, sigma=0.0)
#gp_reg = ExpSqRegularizer(num_outputs,5,lambda=0.25, sigma=0.0)
#gp_reg.covar[1:7,1:7]

ppca = PPCAGridded(num_inputs,num_outputs, num_obs, lambda_range) # Create a PPCA model with random initialization
update_sigma!(gp_reg, 0.0)                               # No regularization
num_training_itterations = 11  # First test that it's working
dataset = Base.Iterators.repeated((ppca,spec),num_training_itterations)
#dataset = Base.Iterators.repeated((ppca,spec),num_training_itterations)
opt = NADAM()          # Choose Optimization algorithm
loss_tol = 0.01
itterations = 0
last_loss = Inf
consecutive_passes = 0
evalcb = function ()   # Callback to show loss and stop optimization
    global last_loss, loss_tol, consecutive_passes
    cs = chisq(ppca,spec)
    reg = regularization(ppca)
    l = cs + reg
    global itterations += 1
    if itterations%10 == 1
        println( string(itterations) * ": χ^2 = ", cs, " reg = ", reg, " loss = ", l)
    end
    num_dof = size(spec.flux,1)*size(spec.flux,2) - get_num_param(ppca) - get_num_inputs(ppca)*get_num_obs(spec)
    #if cs < min( 0.1*num_dof, num_dof - 2.0*sqrt(num_dof) )
    if (last_loss < l+loss_tol) && ( l <= last_loss)
        consecutive_passes += 1
        println( string(itterations) * ": χ^2 = ", cs, " reg = ", reg, " loss = ", l)
        if consecutive_passes >=10
            println("Hurray!")
            Flux.stop()
        end
    else
        consecutive_passes = 0
    end
    last_loss = l
end
#ppca_par = Flux.params(ppca,x)

ppca_par = Flux.params(ppca)
delete!(ppca_par, ppca.x)      # Fix x at true values
ppca.x .= x_true[1:3,:]
Flux.train!(loss, ppca_par, dataset, opt, cb = evalcb)

mean(abs.(ppca.x.-x_true[1:3,:]))

update_sigma!(gp_reg, 0.0)
num_training_itterations = 10000
dataset = Base.Iterators.repeated((ppca,spec),num_training_itterations)
Flux.train!(loss, ppca_par, dataset, opt, cb = evalcb)

plt1 = plot(ppca_true.lambda,ppca_true.b, label="true", title="Mean Vector")
plot!(plt1,ppca.lambda,ppca.b, label="est")
plt2 = plot(ppca.lambda,ppca.b.-ppca_true.b, label="resid")
plot(plt1, plt2, layout =  @layout [a ; b ] )

idx_plt = 1:(floor(Int64,size(ppca_true.W,1)/2))
plt3 = plot(ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,1]./(sqrt(2)*std(ppca_true.W[idx_plt,1])), color=:blue, title="Basis Vectors", legend=:none)
scatter!(plt3, ppca.lambda[idx_plt],ppca.W[idx_plt,1]./(sqrt(2)*std(ppca.W[idx_plt,1])), ms=2, color=:blue)
plot!(plt3,ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,2]./(sqrt(2)*std(ppca_true.W[idx_plt,2])), color=:red)
scatter!(plt3,ppca.lambda[idx_plt],ppca.W[idx_plt,2]./(sqrt(2)*std(ppca.W[idx_plt,2])), ms=2,color=:red)
plot!(plt3,ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])), color=:green)
scatter!(plt3,ppca.lambda[idx_plt],ppca.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])), ms=2, color=:green)
plt4 = plot( ppca.lambda[idx_plt],ppca.W[idx_plt,1]./(sqrt(2)*std(ppca.W[idx_plt,1])).-ppca_true.W[idx_plt,1]./(sqrt(2)*std(ppca_true.W[idx_plt,1])), label="W1", color=:blue)
plot!(plt4,ppca.lambda[idx_plt],ppca.W[idx_plt,2]./(sqrt(2)*std(ppca.W[idx_plt,2])).-ppca_true.W[idx_plt,2]./(sqrt(2)*std(ppca_true.W[idx_plt,2])), label="W2", color=:red)
plot!(plt4,ppca.lambda[idx_plt],ppca.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])).-ppca_true.W[idx_plt,3]./(sqrt(2)*std(ppca_true.W[idx_plt,3])), label="W3", color=:green)
plot(plt3, plt4, layout =  @layout [a ; b ] )

std(ppca.b.-ppca_true.b), [std(ppca.W[:,i]./(sqrt(2)*std(ppca.W[:,i])).-ppca_true.W[:,i]./(sqrt(2)*std(ppca_true.W[:,i]))) for i in 1:3]

ppca_par = Flux.params(ppca)    # Make x active again
update_sigma!(gp_reg, 0.0)
num_training_itterations = 10000
dataset = Base.Iterators.repeated((ppca,spec),num_training_itterations)
Flux.train!(loss, ppca_par, dataset, opt, cb = evalcb)

std(ppca.b.-ppca_true.b), [std(ppca.W[:,i]./(sqrt(2)*std(ppca.W[:,i])).-ppca_true.W[:,i]./(sqrt(2)*std(ppca_true.W[:,i]))) for i in 1:3]

mean(abs.(ppca.x.-x_true[1:3,:]))

plt1 = plot(ppca_true.lambda,ppca_true.b, label="true", title="Mean Vector")
plot!(plt1,ppca.lambda,ppca.b, label="est")
plt2 = plot(ppca.lambda,ppca.b.-ppca_true.b, label="resid")
plot(plt1, plt2, layout =  @layout [a ; b ] )

idx_plt = 1:(floor(Int64,size(ppca_true.W,1)/2))
plt3 = plot(ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,1]./(sqrt(2)*std(ppca_true.W[idx_plt,1])), color=:blue, title="Basis Vectors", legend=:none)
scatter!(plt3, ppca.lambda[idx_plt],ppca.W[idx_plt,1]./(sqrt(2)*std(ppca.W[idx_plt,1])), ms=2, color=:blue)
plot!(plt3,ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,2]./(sqrt(2)*std(ppca_true.W[idx_plt,2])), color=:red)
scatter!(plt3,ppca.lambda[idx_plt],ppca.W[idx_plt,2]./(sqrt(2)*std(ppca.W[idx_plt,2])), ms=2,color=:red)
plot!(plt3,ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])), color=:green)
scatter!(plt3,ppca.lambda[idx_plt],ppca.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])), ms=2, color=:green)
plt4 = plot( ppca.lambda[idx_plt],ppca.W[idx_plt,1]./(sqrt(2)*std(ppca.W[idx_plt,1])).-ppca_true.W[idx_plt,1]./(sqrt(2)*std(ppca_true.W[idx_plt,1])), label="W1", color=:blue)
plot!(plt4,ppca.lambda[idx_plt],ppca.W[idx_plt,2]./(sqrt(2)*std(ppca.W[idx_plt,2])).-ppca_true.W[idx_plt,2]./(sqrt(2)*std(ppca_true.W[idx_plt,2])), label="W2", color=:red)
plot!(plt4,ppca.lambda[idx_plt],ppca.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])).-ppca_true.W[idx_plt,3]./(sqrt(2)*std(ppca_true.W[idx_plt,3])), label="W3", color=:green)
plot(plt3, plt4, layout =  @layout [a ; b ] )

std(ppca.b.-ppca_true.b), [std(ppca.W[:,i]./(sqrt(2)*std(ppca.W[:,i])).-ppca_true.W[:,i]./(sqrt(2)*std(ppca_true.W[:,i]))) for i in 1:3]

#update_lambda!(gp_reg, 0.2)
update_sigma!(gp_reg, 3)
num_training_itterations = 10000
#ppca = PPCAGridded(num_inputs,num_outputs, num_obs, lambda_range)  # Uncomment to start fresh
dataset = Base.Iterators.repeated((ppca,spec),num_training_itterations)
ppca_par = params(ppca)
itterations = 0
#loss_tol = 0.01
last_loss = Inf
Flux.train!(loss, ppca_par, dataset, opt, cb = evalcb)

std(ppca.b.-ppca_true.b), [std(ppca.W[:,i]./(sqrt(2)*std(ppca.W[:,i])).-ppca_true.W[:,i]./(sqrt(2)*std(ppca_true.W[:,i]))) for i in 1:3]

idx_plt = 1:(floor(Int64,size(ppca_true.W,1)/1))
plt7 = plot(ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,1], color=:blue, label=:none)
scatter!(plt7, ppca_true.lambda[idx_plt],ppca.W[idx_plt,1]./(sqrt(2)*std(ppca.W[idx_plt,1])), ms=2, color=:blue)
plot!(plt7,ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,2], color=:red, label=:none)
scatter!(plt7,ppca_true.lambda[idx_plt],ppca.W[idx_plt,2]./(sqrt(2)*std(ppca.W[idx_plt,2])), ms=2,color=:red)
plot!(plt7,ppca_true.lambda[idx_plt],ppca_true.W[idx_plt,3], color=:green, label=:none)
scatter!(plt7,ppca_true.lambda[idx_plt],ppca.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])), ms=2, color=:green)
plt8 = plot( ppca_true.lambda[idx_plt],ppca.W[idx_plt,1]./(sqrt(2)*std(ppca.W[idx_plt,1])).-ppca_true.W[idx_plt,1]./(sqrt(2)*std(ppca_true.W[idx_plt,1])),  color=:blue, label=:none, title="Basis Vector Residuals (w/ reg)")
plot!(plt8,ppca_true.lambda[idx_plt],ppca.W[idx_plt,2]./(sqrt(2)*std(ppca.W[idx_plt,2])).-ppca_true.W[idx_plt,2]./(sqrt(2)*std(ppca_true.W[idx_plt,2])), color=:red, label=:none)
plot!(plt8,ppca_true.lambda[idx_plt],ppca.W[idx_plt,3]./(sqrt(2)*std(ppca.W[idx_plt,3])).-ppca_true.W[idx_plt,3]./(sqrt(2)*std(ppca_true.W[idx_plt,3])),  color=:green, label=:none)
title!(plt4,"Basis Vector Residuals (w/o reg)")
plot(plt4, plt8, layout =  @layout [a ; b ] )

plt5 = plot(ppca_true.lambda,ppca_true.b, label=:none) # "true")
plot!(plt5,ppca.lambda,ppca.b, label=:none) #"est")
plt6 = plot(ppca.lambda,ppca.b.-ppca_true.b, label=:none) # "resid")
plot(plt1, plt2, plt5, plt6, layout =  @layout [a b ; c d ] )

std(ppca.b.-ppca_true.b), [std(ppca.W[:,i]./(sqrt(2)*std(ppca.W[:,i])).-ppca_true.W[:,i]./(sqrt(2)*std(ppca_true.W[:,i]))) for i in 1:3]

mean(abs.(ppca.x.-x_true[1:3,:]))

function chisq_versus_reg_const(reg_const)
    ppca_tmp = copy(ppca)
    update_sigma!(gp_reg, reg_const)
    num_training_itterations = 10000
    #ppca = PPCAGridded(num_inputs,num_outputs, lambda_range)
    #dataset = Base.Iterators.repeated((x,ppca,spec),num_training_itterations)
    dataset = Base.Iterators.repeated((ppca_tmp,spec),num_training_itterations)
    #ppca_par = params(ppca,x)
    ppca_par = params(ppca_tmp)
    itterations = 0
    loss_tol = 0.05
    last_loss = Inf
    Flux.train!(loss, ppca_par, dataset, opt, cb = evalcb)
    #cs = chisq(x,ppca,spec)
    cs = chisq(ppca.x,ppca_tmp,spec)
end

reg_const_list = [0.3, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20. ]
#update_lambda!(gp_reg, 0.2)
chisq_list = chisq_versus_reg_const.(reg_const_list)
reg_const_list[1] = 0.1 # kludge for log
chisq_list

plt_idx = 1:length(reg_const_list) - 0
scatter(log10.(reg_const_list[plt_idx]),log10.(chisq_list[plt_idx]),legend=:none)
xlabel!("log_10 Regularization Constant")
ylabel!("log_10 chi^2")



struct PPCAInterpolated <: AbstractPPCA
  W::Array{Float64,2}
  b::Array{Float64,1}
  x::Array{Float64,2}

  lambda::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
  w_lo_array::Array{Float64,2}
  idx_lo_array::Array{Int64,2}

  function PPCAInterpolated(W::AA2, b::AA1, x::AA2, lambda_range::R; num_pad::Integer = 0) where { T<:Real, AA1<:AbstractArray{T,1}, AA2<:AbstractArray{T,2}, R<:AbstractRange }
     @assert size(W,1) == size(b,1)
     @assert size(W,2) == size(x,1)
     @assert length(lambda_range) == size(b,1)
     num_out = size(b,1) + 2 * num_pad
     num_obs = size(x,2)
     #@assert size(w_lo_array,1) ==  size(b,1)
     #@assert size(w_lo_array,2) ==  size(b,1)
     #@assert size(idx_lo_array) == size(w_lo_array)
     if num_pad == 0
        return new(W, b, x, lambda_range, zeros(Int64, num_out,num_obs), zeros(num_out,num_obs) )
     else
        @error "Need to write code to automatically pad PPCAInterpoalted's grid"
     end

  end


  function PPCAInterpolated(in::PPCAInterpolated)
     @assert size(in.W,1) == size(in.b,1)
     @assert size(in.W,2) == size(in.x,1)
     @assert length(in.lambda) == size(in.b,1)
    new(in.W,in.b,in.x,in.lambda,in.w_lo_array,in.idx_lo_array)
  end
end


function PPCAInterpolated(num_in::Integer, num_out::Integer, lambda_range::R, num_obs::Integer) where {R<:AbstractRange}
    @assert length(lambda_range) == num_out
    PPCAInterpolated(randn(num_out, num_in), randn(num_out), lambda_range, num_obs )
end

get_num_inputs(m::PPCAInterpolated) = size(m.W,2)
get_num_outputs(m::PPCAInterpolated) = size(m.W,1)
get_num_obs(m::PPCAInterpolated) = size(m.x,2)
get_num_param(m::PPCAInterpolated) = get_num_outputs(m)*(get_num_inputs(m)+1)+get_num_obs(m)*get_num_inputs(m)

Flux.@functor PPCAInterpolated (W,b)

function myfindlast(x::T; r::R) where {T<:Real, R<:StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}} }
    idx_lo = min(max(1,ceil(Int64,(x-r.ref.hi)/r.step.hi)),length(r)-1)
end

function myfindlast(x::T; r::R) where { T<:Real, R<:AbstractArray{T} }
    idx_lo = min(max(1,ceil(Int64,convert(T,(x-r.ref)/r.step))),length(r)-1)
end


function calc_interp_w_lo(x, idx_lo; r::AR) where {T<:Real, AR<:AbstractRange }
    idx_hi = idx_lo + 1
    denom = r.step
    w = (r[idx_hi]-x)/denom
end

function calc_interp_w_lo(x, idx_lo; r::AA) where {T<:Real, AA<:AbstractArray{T} }
    idx_hi = idx_lo + 1
    denom = (r[idx_hi]-r[idx_lo])
    w = (r[idx_hi]-x)/denom
end

function update_weights!(m::PPCAInterpolated, lambda::Array{T,2}) where { T<:Real }
    m.idx_lo_array .= myfindlast.(lambda, r=m.lambda)
    m.w_lo_array   .= calc_interp_w_lo.(lambda, m.idx_lo_array, r=m.lambda)
    return m
end

function predict(m::PPCAInterpolated, x::AbstractArray{T,2}, i::Integer) where {T<:Real}
    @assert length(m.w_lo_array) >= 1
    idx = view(m.idx_lo_array,:,i)
    w = view(m.w_lo_array,:,i)
    #=
    Wtmp = w.*m.W[idx,:] + (1.0.-w).*m.W[idx.+1,:]
    btmp = w.*m.b[idx]   + (1.0.-w).*m.b[idx.+1]
    Wtmp*x[:,i].+btmp
    =#
          w .*(m.W[idx,:]   *x[:,i].+m.b[idx]   ) +
    (1.0.-w).*(m.W[idx.+1,:]*x[:,i].+m.b[idx.+1])
end

function predict(m::PPCAInterpolated, i::Integer) where {T<:Real}
    predict(m,m.x,i)
end

function (m::PPCAInterpolated)(x::AbstractArray{T,2})  where {T<:Real}
    mapreduce(i->predict(m,x,i), hcat,1:get_num_obs(m))  # WARNING: mapreduce can't be Autodiffed by Flux
end

function chisq(x::AbstractArray{T,2}, m::PPCAInterpolated, ss::SpectraSet1D)  where { T<:Real }
    @assert get_num_obs(m) == get_num_obs(ss)
    χ2 = zero(T)
    for i in 1:get_num_obs(m)
        yhat = predict(m,i)
        χ2 += sum(((ss.flux[:,i].-yhat).^2 ./ss.flux_var[:,i]))
    end
    return χ2
end

chisq(m::PPCAInterpolated, ss::SpectraSet1D) = chisq(m.x,m,ss)

function regularization(x::AbstractArray{T,2}, m::PPCAInterpolated, ss::SpectraSet1D)  where { T<:Real }
    regularization_scores(x) + regularization_basis(m.W) + regularization_mean(m.b, mu=1.0)
end

regularization(m::PPCAInterpolated, ss::SpectraSet1D) = regularization(m.x,m,ss)


ppca_interp = PPCAInterpolated(ppca.W,ppca.b,ppca.x,ppca.lambda, num_pad=0 )
update_weights!(ppca_interp, spec.lambda_ssb)

maximum(abs.(ppca_interp(ppca_interp.x).-ppca(ppca.x)))

chisq(ppca.x,ppca,spec), chisq(ppca_interp.x,ppca_interp,spec)

lmin, lmax = ppca.lambda[1], ppca.lambda[end]
lambda_pert = ppca.lambda .+ convert(Float64,5*ppca.lambda.step) .* sin.(2π*(ppca.lambda.-lmin)./(lmax-lmin))
spec.lambda_ssb .= lambda_pert
update_weights!(ppca_interp, spec.lambda_ssb)

chisq(ppca.x,ppca,spec), chisq(ppca_interp.x,ppca_interp,spec)

update_sigma!(gp_reg, 4.0)                               # No regularization
num_training_itterations = 11  # First test that it's working
dataset = Base.Iterators.repeated((ppca_interp,spec),num_training_itterations)
opt = NADAM()          # Choose Optimization algorithm
loss_tol = 0.01
itterations = 0
last_loss = Inf
consecutive_passes = 0
evalcb = function ()   # Callback to show loss and stop optimization
    global last_loss, loss_tol, consecutive_passes
    cs = chisq(ppca_interp.x,ppca_interp,spec)
    reg = regularization(ppca_interp.x,ppca_interp,spec)
    l = cs + reg
    global itterations += 1
    if itterations%10 == 1
        println( string(itterations) * ": χ^2 = ", cs, " reg = ", reg, " loss = ", l)
    end
    num_dof = size(spec.flux,1)*size(spec.flux,2) - get_num_param(ppca_interp) - get_num_inputs(ppca_interp)*get_num_obs(spec)
    #if cs < min( 0.1*num_dof, num_dof - 2.0*sqrt(num_dof) )
    if (last_loss < l+loss_tol) && ( l <= last_loss)
        consecutive_passes += 1
        println( string(itterations) * ": χ^2 = ", cs, " reg = ", reg, " loss = ", l)
        if consecutive_passes >=10
            println("Hurray!")
            Flux.stop()
        end
    else
        consecutive_passes = 0
    end
    last_loss = l
end
ppcai_par = params(ppca_interp)
Flux.train!(loss, ppcai_par, dataset, opt, cb = evalcb)

num_training_itterations = 3000
dataset = Base.Iterators.repeated((x,ppca_interp,spec),num_training_itterations)
Flux.train!(loss, ppcai_par, dataset, opt, cb = evalcb)

plt1b = plot(ppca_true.lambda,ppca_true.b, label="true", title="Mean")
plot!(plt1b,ppca_interp.lambda,ppca_interp.b, label="est")
plt2b = plot(ppca_true.lambda,ppca_interp.b.-ppca_true.b, label="resid")
plot(plt1b, plt2b, layout =  @layout [a ; b ] )

plt3b = plot(ppca_true.lambda,ppca_true.W[:,1], color=:blue, title="Basis Vectors (after λ transform)")
scatter!(plt3b,ppca.lambda,ppca.W[:,1]./(sqrt(2)*std(ppca.W[:,1])), ms=2, color=:blue)
plot!(plt3b,ppca_true.lambda,ppca_true.W[:,2], color=:red)
scatter!(plt3b,ppca.lambda,ppca.W[:,2]./(sqrt(2)*std(ppca.W[:,2])), ms=2,color=:red)
plot!(plt3b,ppca_true.lambda,ppca_true.W[:,3], color=:green)
scatter!(plt3b,ppca.lambda,ppca.W[:,3]./(sqrt(2)*std(ppca.W[:,3])), ms=2, color=:green)
plt4b = plot( ppca_true.lambda,ppca.W[:,1]./(sqrt(2)*std(ppca.W[:,1])).-ppca_true.W[:,1]./(sqrt(2)*std(ppca_true.W[:,1])), label="W1", color=:blue)
plot!(plt4b,ppca_true.lambda,ppca.W[:,2]./(sqrt(2)*std(ppca.W[:,2])).-ppca_true.W[:,2]./(sqrt(2)*std(ppca_true.W[:,2])), label="W2", color=:red)
plot!(plt4b,ppca_true.lambda,ppca.W[:,3]./(sqrt(2)*std(ppca.W[:,3])).-ppca_true.W[:,3]./(sqrt(2)*std(ppca_true.W[:,3])), label="W3", color=:green)
plot(plt3b, plt4b, layout =  @layout [a ; b ] )

plt5b = plot(ppca_true.lambda,ppca_true.b, label=:none)
plot!(plt5b,ppca.lambda,ppca_interp.b, label=:none)
plt6b = plot(ppca.lambda,ppca_interp.b.-ppca_true.b, label=:none)
title!(plt1b, "Mean before transform")
title!(plt5b, "Mean after transform")
plot(plt1b, plt2b, plt5b, plt6b, layout =  @layout [a b ; c d ] )

plt7b = plot(ppca_true.lambda,ppca_true.W[:,1], color=:blue, label=:none, title="Basis Vectors")
scatter!(plt7b, ppca_true.lambda,ppca_interp.W[:,1]./(sqrt(2)*std(ppca_interp.W[:,1])), ms=2, color=:blue)
plot!(plt7b,ppca_true.lambda,ppca_true.W[:,2], color=:red, label=:none)
scatter!(plt7b,ppca_true.lambda,ppca_interp.W[:,2]./(sqrt(2)*std(ppca_interp.W[:,2])), ms=2,color=:red)
plot!(plt7b,ppca_true.lambda,ppca_true.W[:,3], color=:green, label=:none)
scatter!(plt7b,ppca_true.lambda,ppca_interp.W[:,3]./(sqrt(2)*std(ppca_interp.W[:,3])), ms=2, color=:green)
plt8b = plot( ppca_true.lambda,ppca_interp.W[:,1]./(sqrt(2)*std(ppca_interp.W[:,1])).-ppca_true.W[:,1]./(sqrt(2)*std(ppca_true.W[:,1])),  color=:blue, label=:none)
plot!(plt8b,ppca_true.lambda,ppca_interp.W[:,2]./(sqrt(2)*std(ppca_interp.W[:,2])).-ppca_true.W[:,2]./(sqrt(2)*std(ppca_true.W[:,2])), color=:red, label=:none)
plot!(plt8b,ppca_true.lambda,ppca_interp.W[:,3]./(sqrt(2)*std(ppca_interp.W[:,3])).-ppca_true.W[:,3]./(sqrt(2)*std(ppca_true.W[:,3])),  color=:green, label=:none)
title!(plt4b,"Basis Vector Residuals before transform")
title!(plt8b,"Basis Vector Residuals after transform")
plot(plt4b, plt8b, layout =  @layout [a ; b ] )



struct PPCAGPBasis <: AbstractPPCA
  W::Array{<:Stheno.AbstractGP,1}
  b::Stheno.AbstractGP
  x::Array{Float64,2}
  W_prior::Array{<:Stheno.AbstractGP,1}
  b_prior::Stheno.AbstractGP
  lambda
#=
  function PPCAGPBasis(W::AbstractArray{WelT,1}, b::BT, x::AA2, lambda_range::R; num_pad::Integer = 0) where {
        WelT<:AbsGp,1}, BT<:AbsGp, R<:AbstractRange }
     @assert length(first(W)) == length(b)
     @assert size(W,1) == size(x,1)
     #@assert length(lambda_range) == length(b)
     num_out = length(b)
     num_obs = size(x,2)
     return new(W, b, x, lambda_range )
  end


  function PPCAGPBasis(in::PPCAGPBasis)
    new(in.W,in.b,in.x,in.lambda)
  end
    =#
end


function PPCAGPBasis(num_basis::Integer, lambda_range::R, num_obs::Integer) where {R<:AbstractRange}
    @assert 3 <= length(lambda_range) <= 512
    AbsGP = Stheno.AbstractG
    num_out = length(lambda_range)
    W = Array{AbsGP,1}(undef,num_basis)
    b = GP(k, GPC())
    W_prior = Array{AbsGP,1}(GP(k, GPC()),num_basis)
    b_prior = GP(k, GPC())
    x = randn(num_basis,num_obs)
    PPCAGPBasis(W,b,x,lambda_range, W_prior, b_prior)
end

get_num_inputs(m::PPCAGPBasis) = length(m.lambda) # or length(first(m.W))
get_num_outputs(m::PPCAGPBasis) = length(m.W)
get_num_obs(m::PPCAGPBasis) = size(m.x,2)
get_num_param(m::PPCAGPBasis) = get_num_outputs(m)*(get_num_inputs(m)+1)+get_num_obs(m)*get_num_outputs(m)

Flux.@functor PPCAGPBasis (W,b)

function predict(m::PPCAGPBasis, x::AbstractArray{T,2}, i::Integer, lambda::AbstractArray{T,1} ) where {T<:Real}
    xtmp = view(x,:,i)
    Wtmp = mapreduce(j->mean(m.W[j](lambda)), hcat,1:length(m.W) )
    btmp = mean(m.b(lambda))
    ypred = Wtmp*xtmp+btmp
end

function predict(m::PPCAGPBasis, i::Integer, lambda::AbstractArray{T,1} ) where {T<:Real}
    predict(m,m.x,i, lambda)
end

function (m::PPCAGPBasis)(x::AbstractArray{T,2})  where {T<:Real}
    mapreduce(i->predict(m,x,i), hcat,1:get_num_obs(m))  # WARNING: mapreduce can't be Autodiffed by Flux
end

function chisq(x::AbstractArray{T,2}, m::PPCAGPBasis, ss::SpectraSet1D)  where { T<:Real }
    @assert get_num_obs(m) == get_num_obs(ss)
    χ2 = zero(T)
    for i in 1:get_num_obs(m)
        yhat = predict(m,i,ss.lambda_ssb[:,i])
        χ2 += sum(((ss.flux[:,i].-yhat).^2 ./ss.flux_var[:,i]))
    end
    return χ2
end

chisq(m::PPCAGPBasis, ss::SpectraSet1D) = chisq(m.x,m,ss)


function regularization_mean(m::PPCAGPBasis; mu=zero(T), lambda = m.lambda )  where { T<:Real, GPT<:Stheno.AbstractGP }
    -logpdf(m.b_prior(lambda), mean(m.b(lambda)).-mu)
end

function regularization_basis(m::PPCAGPBasis, i::Integer, lambda = m.lambda)  where { GPT<:Stheno.AbstractGP }
    -logpdf(m.W_prior[i](lambda), mean(m.W[i](lambda)) )
end

function regularization_basis(m::PPCAGPBasis)  where { GPT<:Stheno.AbstractGP }
    sum( [regularization_basis(m,i) for i in 1:get_num_outputs(m) ] )
end

function regularization(x::AbstractArray{T,2}, m::PPCAGPBasis #=, ss::SpectraSet1D =#)  where { T<:Real }
    regularization_scores(x) + regularization_basis(m) + regularization_mean(m, mu=1.0)
end

function regularization(m::PPCAGPBasis #=, ss::SpectraSet1D=# )  where { T<:Real }
    regularization_scores(m.x) + regularization_basis(m) + regularization_mean(m, mu=1.0)
end

#= Old demo code

ktmp = 1.0 * stretch(Matern52(), 1 /1.0)
bGP = GP(1.0, ktmp, GPC())
b_prior = GP(1.0, ktmp, GPC())
W_prior = [GP(1.0, ktmp, GPC()) for j in 1:size(ppca.W,2)]
bGPx = bGP(ppca.lambda, 0.01)
bGPpost = bGP | Obs( bGPx , ppca.b)
WGPpost = [bGP | Obs( bGPx , ppca.W[:,j]) for j in 1:size(ppca.W,2)]
#ppca_gp = PPCAGPBasis(ppca.W,ppca.b,ppca.x,ppca.lambda)
#update_weights!(ppca_interp, spec.lambda_ssb)
#scatter(ppca.lambda, ppca.W[:,1])
plt = plot(ppca.lambda, mean(bGPpost(ppca.lambda)), color=1, label="b")
scatter!(plt,ppca.lambda, ppca.b, label=:none)

for j in 1:size(ppca.W,2)
    plot!(plt,ppca.lambda, mean(WGPpost[j](ppca.lambda)), color = j+1,  label = "W"*string(j))
    scatter!(plt,ppca.lambda, ppca.W[:,j], color=j+1, label=:none)
end
display(plt)

ppca_gp = PPCAGPBasis(WGPpost,bGPpost,ppca.x, W_prior,b_prior,lambda_range)

plt = plot(ppca_gp.lambda, mean(ppca_gp.b(ppca_gp.lambda)), color=1, label="b")
scatter!(plt,ppca.lambda, ppca.b, label=:none)
#=
for j in 1:size(ppca.W,2)
    plot!(plt,ppca.lambda, mean(WGPpost[j](ppca.lambda)), color = j+1,  label = "W"*string(j))
    scatter!(plt,ppca.lambda, ppca.W[:,j], color=j+1, label=:none)
end
=#
display(plt)

mean(ppca_gp.b(ppca_gp.lambda))

@time cs = chisq(ppca_gp.x,ppca_gp,spec)

#update_sigma!(gp_reg, 0.0)                               # No regularization
num_training_itterations = 1  # First test that it's working
dataset = Base.Iterators.repeated((ppca_gp,spec),num_training_itterations)
opt = NADAM()          # Choose Optimization algorithm
loss_tol = 0.01
itterations = 0
last_loss = Inf
consecutive_passes = 0
evalcb = function ()   # Callback to show loss and stop optimization
    global last_loss, loss_tol, consecutive_passes
    cs = chisq(ppca_gp.x,ppca_gp,spec)
    reg = regularization(ppca_gp.x,ppca_gp)
    l = cs + reg
    global itterations += 1
    if itterations%10 == 1
        println( string(itterations) * ": χ^2 = ", cs, " reg = ", reg, " loss = ", l)
    end
    num_dof = size(spec.flux,1)*size(spec.flux,2) - get_num_param(ppca_gp) - get_num_inputs(ppca_gp)*get_num_obs(spec)
    #if cs < min( 0.1*num_dof, num_dof - 2.0*sqrt(num_dof) )
    if (last_loss < l+loss_tol) && ( l <= last_loss)
        consecutive_passes += 1
        println( string(itterations) * ": χ^2 = ", cs, " reg = ", reg, " loss = ", l)
        if consecutive_passes >=10
            println("Hurray!")
            Flux.stop()
        end
    else
        consecutive_passes = 0
    end
    last_loss = l
end
ppca_gp_par = params(ppca_gp)
Flux.train!(loss, ppca_gp_par, dataset, opt, cb = evalcb)

num_training_itterations = 10
dataset = Base.Iterators.repeated((x,ppca_gp,spec),num_training_itterations)
@time Flux.train!(loss, ppca_gp_par, dataset, opt, cb = evalcb)

plt1b = plot(ppca_true.lambda,ppca_true.b, label="true", title="Mean")
plot!(plt1b,ppca_gp.lambda,ppca_gp.b, label="est")
plt2b = plot(ppca_true.lambda,ppca_gp.b.-ppca_true.b, label="resid")
plot(plt1b, plt2b, layout =  @layout [a ; b ] )

plt3b = plot(ppca_true.lambda,ppca_true.W[:,1], color=:blue, title="Basis Vectors (after λ transform)")
scatter!(plt3b,ppca.lambda,ppca.W[:,1]./(sqrt(2)*std(ppca.W[:,1])), ms=2, color=:blue)
plot!(plt3b,ppca_true.lambda,ppca_true.W[:,2], color=:red)
scatter!(plt3b,ppca.lambda,ppca.W[:,2]./(sqrt(2)*std(ppca.W[:,2])), ms=2,color=:red)
plot!(plt3b,ppca_true.lambda,ppca_true.W[:,3], color=:green)
scatter!(plt3b,ppca.lambda,ppca.W[:,3]./(sqrt(2)*std(ppca.W[:,3])), ms=2, color=:green)
plt4b = plot( ppca_true.lambda,ppca.W[:,1]./(sqrt(2)*std(ppca.W[:,1])).-ppca_true.W[:,1]./(sqrt(2)*std(ppca_true.W[:,1])), label="W1", color=:blue)
plot!(plt4b,ppca_true.lambda,ppca.W[:,2]./(sqrt(2)*std(ppca.W[:,2])).-ppca_true.W[:,2]./(sqrt(2)*std(ppca_true.W[:,2])), label="W2", color=:red)
plot!(plt4b,ppca_true.lambda,ppca.W[:,3]./(sqrt(2)*std(ppca.W[:,3])).-ppca_true.W[:,3]./(sqrt(2)*std(ppca_true.W[:,3])), label="W3", color=:green)
plot(plt3b, plt4b, layout =  @layout [a ; b ] )

plt5b = plot(ppca_true.lambda,ppca_true.b, label=:none)
plot!(plt5b,ppca.lambda,ppca_interp.b, label=:none)
plt6b = plot(ppca.lambda,ppca_interp.b.-ppca_true.b, label=:none)
title!(plt1b, "Mean before transform")
title!(plt5b, "Mean after transform")
plot(plt1b, plt2b, plt5b, plt6b, layout =  @layout [a b ; c d ] )

plt7b = plot(ppca_true.lambda,ppca_true.W[:,1], color=:blue, label=:none, title="Basis Vectors")
scatter!(plt7b, ppca_true.lambda,ppca_interp.W[:,1]./(sqrt(2)*std(ppca_interp.W[:,1])), ms=2, color=:blue)
plot!(plt7b,ppca_true.lambda,ppca_true.W[:,2], color=:red, label=:none)
scatter!(plt7b,ppca_true.lambda,ppca_interp.W[:,2]./(sqrt(2)*std(ppca_interp.W[:,2])), ms=2,color=:red)
plot!(plt7b,ppca_true.lambda,ppca_true.W[:,3], color=:green, label=:none)
scatter!(plt7b,ppca_true.lambda,ppca_interp.W[:,3]./(sqrt(2)*std(ppca_interp.W[:,3])), ms=2, color=:green)
plt8b = plot( ppca_true.lambda,ppca_interp.W[:,1]./(sqrt(2)*std(ppca_interp.W[:,1])).-ppca_true.W[:,1]./(sqrt(2)*std(ppca_true.W[:,1])),  color=:blue, label=:none)
plot!(plt8b,ppca_true.lambda,ppca_interp.W[:,2]./(sqrt(2)*std(ppca_interp.W[:,2])).-ppca_true.W[:,2]./(sqrt(2)*std(ppca_true.W[:,2])), color=:red, label=:none)
plot!(plt8b,ppca_true.lambda,ppca_interp.W[:,3]./(sqrt(2)*std(ppca_interp.W[:,3])).-ppca_true.W[:,3]./(sqrt(2)*std(ppca_true.W[:,3])),  color=:green, label=:none)
title!(plt4b,"Basis Vector Residuals before transform")
title!(plt8b,"Basis Vector Residuals after transform")
plot(plt4b, plt8b, layout =  @layout [a ; b ] )

=#


end # module
