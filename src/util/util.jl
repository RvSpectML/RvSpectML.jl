const speed_of_light_mps = 299792458.0 # TODO: Update value
"""Compute Doppler boost factor (non-relativistic) for rv in km/s"""
calc_doppler_factor(rv::Real) = one(rv) + rv/speed_of_light_mps

"""Compute Doppler boost factor (relativistic) for rv and v_perp in km/s"""
calc_doppler_factor(rv::Real, v_perp::Real) = (one(rv) + rv/speed_of_light_mps)/(one(rv) - (rv^2+v_perp^2)/speed_of_light_mps^2)

"""Multiply spectra's λs by doppler_factor and update spectra metadata, so doppler_factor knows how to undo the transform."""
function apply_doppler_factor!(spectra::S,doppler_factor::Real) where {S<:AbstractSpectra}
    spectra.λ .*= doppler_factor
    if haskey(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor
    end
    return spectra
end


"""Estimate line width based on stellar Teff (K) and optionally v_rot (km/s).  Output in km/s."""
function predict_line_width(Teff::Real; v_rot::Real=zero(Teff))
    @assert 3000 < Teff < 10000 # K
    @assert 0 <= v_rot <=100 # km/s
    line_width_thermal = 13*sqrt(Teff/1e4) # km/s
    line_width = sqrt(v_rot^2+line_width_thermal^2) # km/s
end


# Find indicies for pixels around lines

const Δλoλ_fit_line_default = 5*(1.8*1000/speed_of_light_mps)
const Δλoλ_edge_pad_default = 0*(1.8*1000/speed_of_light_mps)

function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_center::Real; Δ::Real = Δλoλ_fit_line_default) where T<:Real
    findfirst(x->x>=line_center*(1-Δ),wavelengths):findlast(x->x<=line_center*(1+Δ),wavelengths)
end

function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_lo::Real, line_hi::Real; Δ::Real = Δλoλ_edge_pad_default) where T<:Real
    @assert line_lo < line_hi
    findfirst(x->x>=line_lo*(1-Δ),wavelengths):findlast(x->x<=line_hi*(1+Δ),wavelengths)
end

function find_orders_with_line(goal::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal<=order_max(i), 1:size(lambda,2) )
end

function find_orders_with_line(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal_lo && goal_hi<=order_max(i), 1:size(lambda,2) )
end

function findall_line(goal::Real,lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal <= lambda[end,end]
    orders = find_orders_with_line(goal,lambda)
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal,Δ=Δ),order=o), orders)
    locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal_lo < goal_hi <= lambda[end,end]
    orders = find_orders_with_line(goal_lo,goal_hi,lambda)
    #if ! (length(orders) >= 1) return end
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal_lo, goal_hi,Δ=Δ),order=o), orders)
    locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra}
    findall_line(goal,spectra.λ,spectra.var, Δ=Δ)
end

function findall_line(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra}
    findall_line(goal_lo,goal_hi,spectra.λ,spectra.var, Δ=Δ)
end

function calc_snr(flux::AbstractArray{T1},var::AbstractArray{T2}) where {T1<:Real, T2<:Real}
    @assert size(flux) == size(var)
    sum(flux./var)/sum(1.0 ./ var)
end

function find_line_best(goal::Real,lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal,lambda,var,Δ=Δ)
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal_lo,goal_hi,lambda,var,Δ=Δ)
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra}
    find_line_best(goal,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end

function find_line_best(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra}
    find_line_best(goal_lo,goal_hi,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end

# Manipulation data in chunks of spectrua

#=
function make_chunck_list(spectra::AS, line_list::AA; Δ::Real=Δλoλ_fit_line_default) where { AS<:AbstractSpectra, R<:Real, AA<:AbstractArray{R,1} }
   ChunckList(map(l->ChunkOfSpectra(spectra,find_line_best(l,spectra,Δ=Δ)), line_list) )
end
=#

function make_chunck_list(spectra::AS, line_list::DataFrame; Δ::Real=Δλoλ_edge_pad_default) where { AS<:AbstractSpectra }
    @assert haskey(line_list,:lambda_lo)
    @assert haskey(line_list,:lambda_hi)
    ChunckList(map(row->ChunkOfSpectra(spectra,find_line_best(row.lambda_lo,row.lambda_hi,spectra,Δ=Δ)), eachrow(line_list) ))
end

function make_orders_into_chunks(spectra::AS; orders_to_use=1:size(spectra.flux,2),
        min_col::Integer=min_col_neid_default, max_col::Integer=max_col_neid_default ,
        pixels_to_use=fill(min_col:max_col,length(orders_to_use)) ) where {AS<:AbstractSpectra, }
    ChunckList(map(order->ChunkOfSpectra(spectra,(pixels=pixels_to_use[order],order=order)), orders_to_use ))
end

function filter_bad_chunks(chunk_list_timeseries::ACLT, line_list::DataFrame; verbose::Union{Int,Bool} = false) where { ACLT<:AbstractChunckListTimeseries }
    @assert(haskey(line_list,:lambda_lo))
    @assert(haskey(line_list,:lambda_hi))
    @assert(length(chunk_list_timeseries)>=1)
    idx_keep = trues(num_chunks(chunk_list_timeseries))
    for t in 1:length(chunk_list_timeseries)
        idx_bad_λ = findall(c->any(isnan.(chunk_list_timeseries.chuck_list[t].data[c].λ)),1:num_chunks(chunk_list_timeseries))
        idx_bad_flux = findall(c->any(isnan.(chunk_list_timeseries.chuck_list[t].data[c].flux)),1:num_chunks(chunk_list_timeseries))
        idx_bad_var = findall(c->any(isnan.(chunk_list_timeseries.chuck_list[t].data[c].var)),1:num_chunks(chunk_list_timeseries))
        idx_keep[idx_bad_λ] .= false
        idx_keep[idx_bad_flux] .= false
        idx_keep[idx_bad_var] .= false
        if verbose && (length(idx_bad_λ)+length(idx_bad_flux)+length(idx_bad_var) > 0)
               println("# Removing chunks", vcat(idx_bad_λ,idx_bad_flux,idx_bad_var), " at time ", t, " due to NaNs (",length(idx_bad_λ),",",length(idx_bad_flux),",",length(idx_bad_var),").")
        end

    end
    chunks_to_remove = findall(.!idx_keep)
    if length(chunks_to_remove) == 0
        println("# No lines to remove.")
        return (chunk_timeseries=chunk_list_timeseries, line_list=line_list)
    else
        println("# Removing ", length(chunks_to_remove), " chunks due to NaNs.")
        map(c->println("# ",c,": ",line_list.lambda_lo[c]," - ",line_list.lambda_hi[c]),chunks_to_remove)
        new_line_list = line_list[findall(idx_keep),:]
        new_chunk_list_timeseries = [ChunckList(chunk_list_timeseries.chuck_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        return (chunk_timeseries=ChunckListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries), line_list=new_line_list)
    end
end

""" Create a range with equal spacing between points with end points set based on union of all chunks in timeseries.
    Takes ChunkListTimeseries, chunk index, and optional an oversample_factor (1). """
function make_grid_for_chunck(timeseries::ACLT, c::Integer; oversample_factor::Real = 1.0 ) where { ACLT<:AbstractChunckListTimeseries }
    num_obs = length(timeseries.chuck_list)
    λ_min = maximum(minimum(timeseries.chuck_list[t].data[c].λ) for t in 1:num_obs)
    λ_max = minimum(maximum(timeseries.chuck_list[t].data[c].λ) for t in 1:num_obs)
    Δλ_grid_obs = median((timeseries.chuck_list[t].data[c].λ[end]-
            timeseries.chuck_list[t].data[c].λ[1])/
              (length(timeseries.chuck_list[t].data[c].λ)-1) for t in 1:num_obs)
    num_pixels_obs = mean(length(timeseries.chuck_list[t].data[c].λ) for t in 1:num_obs)
    num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
    Δλ_grid_gen = (λ_max-λ_min)/ (num_pixels_gen-1)
    range(λ_min,stop=λ_max,step=Δλ_grid_gen)
end



# Maniuplation Line lists

function merge_lines(line_list::DataFrame)
    chunk_list_df = DataFrame(:lambda_lo=>Float64[],:lambda_hi=>Float64[],
                                :line_λs=>Array{Float64,1}[],:line_depths=>Array{Float64,1}[])
    num_lines = size(line_list,1)
    @assert num_lines >= 2
    lambda_lo_last = line_list[1,:lambda_lo]
    lambda_hi_last = line_list[1,:lambda_hi]
    line_λs = [line_list[1,:lambda_mid]]
    line_depths = [line_list[1,:depth_vald]]
    for i in 2:num_lines
        lambda_lo = line_list[i,:lambda_lo]
        lambda_hi = line_list[i,:lambda_hi]
        if lambda_lo>lambda_hi_last
            push!(chunk_list_df, (lambda_lo_last, lambda_hi_last, line_λs, line_depths))
            (lambda_lo_last, lambda_hi_last) = (lambda_lo, lambda_hi)
            line_λs = [line_list[i,:lambda_mid]]
            line_depths = [line_list[i,:depth_vald]]
        else
            lambda_hi_last = lambda_hi
            push!(line_λs,line_list[i,:lambda_mid])
            push!(line_depths,line_list[i,:depth_vald])
        end
    end
    if chunk_list_df[end,:lambda_hi] != lambda_hi_last
        #push!(chunk_list_df, (lambda_lo_last, lambda_hi_last))
        push!(chunk_list_df, (lambda_lo_last, lambda_hi_last, line_λs, line_depths))
    end
    return chunk_list_df
end

# Manipulate spectra
""" Calc normalization of spectra based on average flux in a ChunkList. """
function calc_normalization(chunk_list::ACL) where { ACL<:AbstractChunckList}
    total_flux = sum(sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end

function normalize_spectrum!(spectrum::ST, scale_fac::Real) where { ST<:AbstractSpectra }
    @assert 0 < scale_fac < Inf
    @assert !isnan(scale_fac^2)
    spectrum.flux .*= scale_fac
    spectrum.var .*= scale_fac^2
    return spectrum
end

function normalize_spectra!(timeseries::ACLT, spectra::AS) where { ACLT<:AbstractChunckListTimeseries, ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    @assert length(timeseries) == length(spectra)
    for t in 1:length(timeseries)
        scale_fac = calc_normalization(timeseries.chuck_list[t])
        println("# t= ",t, " scale_fac= ", scale_fac)
        normalize_spectrum!(spectra[t], scale_fac)
    end
    return timeseries
end


# Unorganized code
function make_interpolator_linear_flux(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra}
    LinearInterpolation(spectra.λ, spectra.flux)
end

function make_interpolator_linear_var(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra}
    LinearInterpolation(spectra.λ, spectra.var)
end

#=
using Stheno

function make_interpolator_gp(spectra::Union{AS,AC}; length_scale::Real = 0.1, σ_scale::Real = 1.0) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra}
    xobs = spectra.λ
    yobs = copy(spectra.flux)
    obs_var = spectra.var
    med_y = mean(yobs)

    # Choose the length-scale and variance of the process.
    σ² = (σ_scale * med_y)^2
    # Construct a kernel with this variance and length scale.
    k = σ² * stretch(Matern52(), 1 / length_scale)

    # Specify a zero-mean GP with this kernel. Don't worry about the GPC object.
    f = GP(med_y, k, GPC())
    fx = f(xobs, obs_var)
    f_posterior = f | Obs(fx, yobs)
end

function interp_to_grid(spectra::Union{AS,AC}, grid::AR) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra, AR<:AbstractRange}
   #grid = chunk_grids[c]
   #make_interpolator_linear(spectra).(grid)
   f_posterior = make_interpolator_gp(spectra,length_scale=6e-5*mean(spectra.λ))
   (mean=mean(f_posterior(grid)),  std=std.(marginals(gp_interp(grid))))
end


function pack_chunks_into_matrix(timeseries::ACLT, chunk_grids::AR) where { ACLT<:AbstractChunckListTimeseries, RT<:AbstractRange, AR<:AbstractArray{RT,1} }
   num_obs = length(timeseries)
   num_λ = sum(length.(chunk_grids))
   flux_matrix = Array{Float64,2}(undef,num_λ,num_obs)
   var_matrix = Array{Float64,2}(undef,num_λ,num_obs)
   λ_vec = Array{Float64,1}(undef,num_λ)
   chunk_map = Array{UnitRange{Int64}}(undef,length(chunk_grids))

   for t in 1:num_obs
        idx_start = 0
        for c in 1:length(chunk_grids)
            if length(chunk_grids[c]) <= 512
                gp_interp = make_interpolator_gp(timeseries.chuck_list[t].data[c],length_scale=1e-4*mean(chunk_grids[c]))
                idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
                flux_matrix[idx,t] .= mean(gp_interp(chunk_grids[c]))
                var_matrix[idx,t] .= var.(marginals(gp_interp(chunk_grids[c])))
                if t == 1
                    λ_vec[idx] .= chunk_grids[c]
                    chunk_map[c] = idx
                    #= if 8500<idx_start <8700
                        println(c,": ",idx)
                    end
                    =#
                end
                idx_start += length(chunk_grids[c])
            else
                lin_interp_flux = make_interpolator_linear_flux(timeseries.chuck_list[t].data[c])
                lin_interp_var = make_interpolator_linear_var(timeseries.chuck_list[t].data[c])
                idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
                flux_matrix[idx,t] .= lin_interp_flux.(chunk_grids[c])
                var_matrix[idx,t] .= lin_interp_var.(chunk_grids[c])
                if t == 1
                    λ_vec[idx] .= chunk_grids[c]
                    chunk_map[c] = idx
                    #= if 8500<idx_start <8700
                        println(c,": ",idx)
                    end
                    =#
                end
                idx_start += length(chunk_grids[c])
            end
        end
    end

    #return flux_matrix
    return (flux=flux_matrix, var=var_matrix, λ=λ_vec, chunk_map=chunk_map)
end

=#

function calc_deriv(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
    @assert size(flux) == size(λ)
    @assert length(flux) >= 3
    dfdlogλ = Array{T1,1}(undef,length(flux))
    dfdlogλ[1] = 0.5*(flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    dfdlogλ[2:end-1] .= 0.5*(flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[end-2])
    dfdlogλ[end] = 0.5*(flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

function calc_mean_deriv(flux::AbstractArray{T1,2}, var::AbstractArray{T1,2}, λ::AbstractArray{T3,1},
        chunk_map::AbstractArray{URT,1}) where
    { T1<:Real, T2<:Real, T3<:Real, URT<:AbstractUnitRange} #, RT<:AbstractRange }
    flux_mean = sum(flux./vm,dims=2)./sum(1.0./vm,dims=2)
    deriv = Array{T1,1}(undef,length(flux_mean))
    map(c->deriv[c] .= calc_deriv(flux_mean[c],λv[c]),chunk_map )
    return deriv
end

@inline function allequal(x::AbstractArray{T,1}) where {T<:Real}
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end
