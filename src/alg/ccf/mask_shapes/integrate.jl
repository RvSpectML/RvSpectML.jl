using QuadGK
using StaticArrays

abstract type AbstractIntegrateFixedOrder end

struct IntegrateFixedOrder{FuncT,Order,NumGW,RetT} <: AbstractIntegrateFixedOrder  where { FuncT, Order<:Integer, RetT<:Real }
    x::SVector{Order,RetT}  # Technically, order = Order-1
    w::SVector{Order,RetT}
    gw::SVector{NumGW,RetT}
    f::FuncT
end

function IntegrateFixedOrder(f::FuncT,order::Integer)  where { FuncT } #, RetT }
    @assert 1<=order <=10
    #ret_type = RetT
    ret_type = typeof(f(0.0))
    x,w,gw = QuadGK.cachedrule(ret_type,order)
    s_x = SVector{order+1,ret_type}(x)
    s_w = SVector{order+1,ret_type}(w)
    num_gw = floor(Int,(order+1)//2)
    s_gw = SVector{num_gw,ret_type}(gw)
    #s_gw = @SVector zeros(floor(Int,(order+1)//2))
    #s_gw[1:length(gw)] .= gw   # since s_gw is generally too big

    IntegrateFixedOrder{FuncT,order+1,num_gw,ret_type}(s_x,s_w,s_gw,f)
end

function integrate_scalar_fixed_order_gauss(ifo::T, a::Real, b::Real)  where { T<:AbstractIntegrateFixedOrder }
    s = convert(eltype(ifo.x), 0.5) * (b-a)
    n1 = 1 - (length(ifo.x) & 1) # 0 if even order, 1 if odd order
    # unroll first iterationof loop to get correct type of Ik and Ig
    fg = ifo.f(a + (1+ifo.x[2])*s) + ifo.f(a + (1-ifo.x[2])*s)
    Ig = fg * ifo.gw[1]
    for i = 2:length(ifo.gw)-n1
        fg = ifo.f(a + (1+ifo.x[2i])*s) + ifo.f(a + (1-ifo.x[2i])*s)
        Ig += fg * ifo.gw[i]
    end
    if n1 != 0 # even: Gauss rule does not include x == 0
        f0 = ifo.f(a + s)
        Ig += f0 * ifo.gw[end]
    end
    Ig_s = Ig * s # new variable since this may change the type
    return Ig_s
end

function integrate_scalar_fixed_order_kronrod(ifo::T, a::Real, b::Real)  where { T<:AbstractIntegrateFixedOrder }
    s = convert(eltype(ifo.x), 0.5) * (b-a)
    n1 = 1 - (length(ifo.x) & 1) # 0 if even order, 1 if odd order
    # unroll first iterationof loop to get correct type of Ik and Ig
    fg = ifo.f(a + (1+ifo.x[2])*s) + ifo.f(a + (1-ifo.x[2])*s)
    fk = ifo.f(a + (1+ifo.x[1])*s) + ifo.f(a + (1-ifo.x[1])*s)
    Ik = fg * ifo.w[2] + fk * ifo.w[1]
    for i = 2:length(ifo.gw)-n1
        fg = ifo.f(a + (1+ifo.x[2i])*s) + ifo.f(a + (1-ifo.x[2i])*s)
        fk = ifo.f(a + (1+ifo.x[2i-1])*s) + ifo.f(a + (1-ifo.x[2i-1])*s)
        Ik += fg * ifo.w[2i] + fk * ifo.w[2i-1]
    end
    if n1 == 0 # even: Gauss rule does not include x == 0
        Ik += f(a + s) * ifo.w[end]
    else # odd: don't count x==0 twice in Gauss rule
        f0 = f(a + s)
        Ik += f0 * ifo.w[end] +
            (f(a + (1+ifo.x[end-1])*s) + f(a + (1-ifo.x[end-1])*s)) * ifo.w[end-1]
    end
    Ik_s = Ik * s # new variable since this may change the type
    return Ik_s
end

function (ifo::IntegrateFixedOrder)(a::Real, b::Real)
    integrate_scalar_fixed_order_gauss(ifo,a,b)
end

function error_estimate(ifo::T, a::Real, b::Real)  where { T<:AbstractIntegrateFixedOrder }
    Ig = integrate_scalar_fixed_order_gauss(ifo,a,b)
    Ik = integrate_scalar_fixed_order_kronrod(ifo,a,b)
    return abs(Ik - Ig)
end
