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
