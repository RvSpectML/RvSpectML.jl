""" `numerical_deriv( x, y)`
Calculate simple estimate of numerical derivative
"""
function numerical_deriv( x::AA1, y::AA2 )   where { T1<:Real, AA1<:AbstractArray{T1,1}, T2<:Real, AA2<:AbstractArray{T2,1}  }
	@assert length(x) == length(y)
	dfluxdlnλ = zeros(size(y))
	dfluxdlnλ[1] = (y[2]-y[1])/(x[2]-x[1])
	dfluxdlnλ[2:end-1] .= (y[3:end].-y[1:end-2])./(x[3:end].-x[1:end-2])
	dfluxdlnλ[end] = (y[end]-y[end-1])/(x[end]-x[end-1])
	return dfluxdlnλ
end

""" Estimate numerical derivative of fluxes given wavelengths. """
function calc_dfluxdlnlambda(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
    @assert size(flux) == size(λ)
    @assert length(flux) >= 3
    dfdlogλ = Array{T1,1}(undef,length(flux))
    calc_dfluxdlnlambda!(dfdlogλ,flux,λ)
    #dfdlogλ[1] = 0.5 * (flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    #dfdlogλ[2:end-1] .= 0.5 .* (flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[1:end-2])
    #dfdlogλ[end] = 0.5 * (flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

""" Estimate numerical derivative of fluxes given wavelengths. """
function calc_dfluxdlnlambda!(dfdlogλ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, λ::AbstractArray{T3,1}) where { T1<:Real, T2<:Real, T3<:Real }
    @assert size(flux) == size(λ)
    @assert size(dfdlogλ) == size(flux)
    @assert length(flux) >= 3
    #dfdlogλ = Array{T1,1}(undef,length(flux))
    dfdlogλ[1] = 0.5 * (flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    dfdlogλ[2:end-1] .= 0.5 .* (flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[1:end-2])
    dfdlogλ[end] = 0.5 * (flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

""" Estimate numerical second derivative of fluxes given wavelengths. """
function calc_d2fluxdlnlambda2(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
    @assert size(flux) == size(λ)
    @assert length(flux) >= 3
    logλ = log.(λ)
    d2fdlogλ2 = Array{T1,1}(undef,length(flux))
    #d2fdlogλ2[2:end-1] .= 0.25*(flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1])./(λ[3:end].+λ[1:end-2].-2.0*λ[2:end-1]).*(λ[3:end].+λ[end-2]).^2
    #d2fdlogλ2[2:end-1] .= 0.5 * (flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1]).* ((λ[3:end].+λ[1:end-2])./(logλ[3:end].-logλ[1:end-2])).^2
    #d2fdlogλ2[end] = d2fdlogλ2[end-1]
    #d2fdlogλ2[1] = d2fdlogλ2[2]
    calc_d2fluxdlnlambda2!(d2fdlogλ2,flux,λ)
    return d2fdlogλ2
end

""" Estimate numerical second derivative of fluxes given wavelengths. """
function calc_d2fluxdlnlambda2!(d2fdlogλ2::AbstractArray{T1,1}, flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real, T3<:Real }
    @assert size(flux) == size(λ)
    @assert size(d2fdlogλ2) == size(flux)
    @assert length(flux) >= 3
    #logλ = log.(λ)
    #d2fdlogλ2 = Array{T1,1}(undef,length(flux))
    #d2fdlogλ2[2:end-1] .= 0.25*(flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1])./(λ[3:end].+λ[1:end-2].-2.0*λ[2:end-1]).*(λ[3:end].+λ[end-2]).^2
    d2fdlogλ2[2:end-1] .= 0.5 * (flux[3:end].+flux[1:end-2].-2.0.*flux[2:end-1]).* ((λ[3:end].+λ[1:end-2]).^2 ./ (2.0.*log.(λ[3:end]./λ[1:end-2])))
    d2fdlogλ2[end] = d2fdlogλ2[end-1]
    d2fdlogλ2[1] = d2fdlogλ2[2]
    return d2fdlogλ2
end

function calc_weights_for_running_mean(n::Integer; decay_factor::Float64 = 2.0)
	@assert 1<=n<=100
	@assert 1 <= decay_factor < Inf
	w = zeros(n)
	if isodd(n)
		center = convert(Int,1 + (n-1)//2)
		w[center] = 1
		for i in 1:(center-1)
			w[center+i] = w[center-i] = decay_factor^(-i)
		end
	else
		center = convert(Int,n//2)
		for i in 1:center
			w[center+i] = w[center+1-i] = decay_factor^(-(i-0.5))
		end
	end
	w ./= sum(w)
	return w
end

""" Smooth array via weighted running mean of n consecutive points """
function calc_running_mean_weighted(x::AbstractVector{T1}; n::Integer = 5, decay_factor::Float64 = 2.0, 
				weights::AbstractVector{T2} = calc_weights_for_running_mean(n, decay_factor=decay_factor) ) where { T1<:Real, T2<:Real }
	@assert 2 <= n < length(x)/2
	@assert isodd(n)   # Could generalize, but haven't
	smoothed = Array{T1,1}(undef,length(x))
	n = length(weights)
	m = floor(Int,n//2)
	for i in 1:m
		start_idx = 2+m-i
		smoothed[i] = sum(view(x,1:i+m).*view(weights,start_idx:n)) /sum(view(weights,start_idx:n))
	end
	for i in (m+1):(length(x)-m-1)
		smoothed[i] = sum(view(x,i-m:i+m) .* weights)
	end
	for i in (length(x)-m):length(x)
		stop_idx = length(x)-(i-m)+1
		smoothed[i] = sum(view(x,i-m:length(x)).*view(weights,1:stop_idx)) /sum(view(weights,1:stop_idx))
	end
	return smoothed
end
