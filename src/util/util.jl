"""
Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

# Constants
const speed_of_light_mps = 299792458.0 # TODO: Update value

"""Compute Doppler boost factor (non-relativistic) for rv in km/s."""
calc_doppler_factor(rv::Real) = one(rv) + rv/speed_of_light_mps

"""Compute Doppler boost factor (relativistic) for rv and v_perp in km/s"""
calc_doppler_factor(rv::Real, v_perp::Real) = (one(rv) + rv/speed_of_light_mps)/(one(rv) - (rv^2+v_perp^2)/speed_of_light_mps^2)

"""Multiply spectra's λs by doppler_factor and update spectra metadata, so doppler_factor knows how to undo the transform."""
function apply_doppler_factor!(spectra::S,doppler_factor::Real) where {S<:AbstractSpectra}
    spectra.λ .*= doppler_factor
    if haskey(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] ./= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor
    end
    return spectra
end

""" Calculate total SNR in (region of) spectra. """
function calc_snr(flux::AbstractArray{T1},var::AbstractArray{T2}) where {T1<:Real, T2<:Real}
    @assert size(flux) == size(var)
    sqrt(sum(flux./var)) # /sum(1.0 ./ var)
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

""" Return a range of columns indices with wavelengths within Δ of line_center """
function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_center::Real; Δ::Real = Δλoλ_fit_line_default) where T<:Real
    @assert Δ >= zero(Δ)
    findfirst(x->x>=line_center*(1-Δ),wavelengths):findlast(x->x<=line_center*(1+Δ),wavelengths)
end

""" Return a range of columns indices with wavelengths between line_lo and line_hi """
function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_lo::Real, line_hi::Real; Δ::Real = Δλoλ_edge_pad_default) where T<:Real
    @assert line_lo < line_hi
    findfirst(x->x>=line_lo*(1-Δ),wavelengths):findlast(x->x<=line_hi*(1+Δ),wavelengths)
end

""" Return list of all orders that contain a pixel with wavelength lambda """
function find_orders_with_line(goal::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal<=order_max(i), 1:size(lambda,2) )
end

""" Return list of all orders that contain all pixels with wavelengths between goal_lo and goal_hi """
function find_orders_with_line(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal_lo && goal_hi<=order_max(i), 1:size(lambda,2) )
end

""" Return list of (pixels, order) pairs that contain pixels with desireed wavelengths.
    Excludes locations that contain any pixels with var == NaN.
"""
function findall_line end

function findall_line(goal::Real,lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal <= lambda[end,end]
    @assert size(lambda) == size(var)
    @assert Δ >= zero(Δ)
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

""" Return (pixels, order) pair that contain "best" region of spectra, based on highest SNR. """
function find_line_best end

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
function make_chunk_list(spectra::AS, line_list::AA; Δ::Real=Δλoλ_fit_line_default) where { AS<:AbstractSpectra, R<:Real, AA<:AbstractArray{R,1} }
   ChunkList(map(l->ChuckOfSpectrum(spectra,find_line_best(l,spectra,Δ=Δ)), line_list) )
end
=#

""" Return a ChunkList of best regions of spectrum with lines in line_line.
    line_list is a DataFrame containing :lambda_lo and :lambda_hi.
    Pads edges by Δ.
"""
function make_chunk_list(spectra::AS, line_list::DataFrame; Δ::Real=Δλoλ_edge_pad_default) where { AS<:AbstractSpectra }
    @assert haskey(line_list,:lambda_lo)
    @assert haskey(line_list,:lambda_hi)
    ChunkList(map(row->ChuckOfSpectrum(spectra,find_line_best(row.lambda_lo,row.lambda_hi,spectra,Δ=Δ)), eachrow(line_list) ))
end

""" Return a ChunkList with a region of spectrum from each order in orders_to_use.
    inst:  Instrument trait that provides orders_to_use and pixels_to_use.
    or
    orders_to_use: Range or Array (orders_to_use(inst))
    pixels_to_use: Array of Ranges (each from min_col to max_col)
    min_col: (min_col_default(inst))
    max_col: (max_col_default(inst))
"""
function make_orders_into_chunks
end

function make_orders_into_chunks(spectra::AS, inst::AbstractSpectra,
        min_col::Integer = min_col_default(inst), # min_col_neid_default,
        max_col::Integer = max_col_default(inst); #max_col_neid_default ,
        orders_to_use = orders_to_use(inst) # 1:size(spectra.flux,2),
        ) where {AS<:AbstractSpectra }
    @assert eltype(orders_to_use) <: Integer
    @assert all( min_order(inst) .<= orders_to_use .<= max_order(inst) )
    @assert min_col >= min_pixel_in_order(inst)
    @assert max_col <= max_pixel_in_order(inst)
    pixels_to_use=fill(min_col:max_col,length(orders_to_use))
    make_orders_into_chunks(spectra,inst,orders_to_use=orders_to_use, pixels_to_use=pixels_to_use)
end

function make_orders_into_chunks(spectra::AS;
        orders_to_use::Union{AR,AA1}, pixels_to_use::AAR ) where {
         AS<:AbstractSpectra, AR<:AbstractRange{Int64}, AA1<:AbstractArray{Int64,1}, AAR<:AbstractArray{AR,1} }
    @assert eltype(orders_to_use) <: Integer
    #@assert all( min_order(inst) .<= orders_to_use .<= max_order(inst) )
    #@assert minimum(pixels_to_use) >= min_pixel_in_order(inst)
    #@assert maximum(pixels_to_use) <= max_pixel_in_order(inst)
    ChunkList( map(order->
                    ChuckOfSpectrum(spectra,(pixels=pixels_to_use[order],order=order)),
                    orders_to_use ))
end

""" Return (chunk_timeseries, line_list) that have been trimmed of any chunks that are bad based on any spectra in the chunk_timeseries.
    chunk_timeseries: ChunkListTimeseries
    line_linst:  DataFrame w/ lambda_lo, lambda_hi
    verbose: print debugging info (false)
"""
function filter_bad_chunks(chunk_list_timeseries::ACLT, line_list::DataFrame; verbose::Union{Int,Bool} = false) where { ACLT<:AbstractChunkListTimeseries }
    @assert(length(chunk_list_timeseries)>=1)
    @assert(haskey(line_list,:lambda_lo))
    @assert(haskey(line_list,:lambda_hi))
    idx_keep = trues(num_chunks(chunk_list_timeseries))
    for t in 1:length(chunk_list_timeseries)
        idx_bad_λ = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].λ)),1:num_chunks(chunk_list_timeseries))
        idx_bad_flux = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].flux)),1:num_chunks(chunk_list_timeseries))
        idx_bad_var = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].var)),1:num_chunks(chunk_list_timeseries))
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
        new_chunk_list_timeseries = [ChunkList(chunk_list_timeseries.chunk_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        return (chunk_timeseries=ChunkListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries), line_list=new_line_list)
    end
end

""" Create a range with equal spacing between points with end points set based on union of all chunks in timeseries.
    timeseries: ChunkListTimeseries
    chunk index:
    oversample_factor: (1)
"""
function make_grid_for_chunk(timeseries::ACLT, c::Integer; oversample_factor::Real = 1.0 ) where { ACLT<:AbstractChunkListTimeseries }
    num_obs = length(timeseries.chunk_list)
    @assert num_obs >= 1
    @assert 1<= c <= length(first(timeseries.chunk_list).data) # WARN: only tests first chunk_list
    λ_min = maximum(minimum(timeseries.chunk_list[t].data[c].λ) for t in 1:num_obs)
    λ_max = minimum(maximum(timeseries.chunk_list[t].data[c].λ) for t in 1:num_obs)
    Δλ_grid_obs = median((timeseries.chunk_list[t].data[c].λ[end]-
            timeseries.chunk_list[t].data[c].λ[1])/
              (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
    num_pixels_obs = mean(length(timeseries.chunk_list[t].data[c].λ) for t in 1:num_obs)
    num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
    Δλ_grid_gen = (λ_max-λ_min)/ (num_pixels_gen-1)
    range(λ_min,stop=λ_max,step=Δλ_grid_gen)
end



# Maniuplation Line lists
""" Return DataFrame with specifications for each chunk which will contain one or more lines.
    Input:  line_list a DataFrame with:
    -  lambda_lo, lambda_hi, lambda_mid, depth_vald
    Output:  DataFrame with
    - lambda_lo & lambda_hi: boundaries for chunk
    - line_λs & line_depths: arrays with info about each line
"""
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
function calc_normalization(chunk_list::ACL) where { ACL<:AbstractChunkList}
    total_flux = sum(sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end

""" Normalize spectrum, multiplying fluxes by scale_fac. """
function normalize_spectrum!(spectrum::ST, scale_fac::Real) where { ST<:AbstractSpectra }
    @assert 0 < scale_fac < Inf
    @assert !isnan(scale_fac^2)
    spectrum.flux .*= scale_fac
    spectrum.var .*= scale_fac^2
    return spectrum
end

""" Normalize each spectrum based on sum of fluxes in chunk_timeseries region of each spectrum. """
function normalize_spectra!(chunk_timeseries::ACLT, spectra::AS) where { ACLT<:AbstractChunkListTimeseries, ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    @assert length(chunk_timeseries) == length(spectra)
    for t in 1:length(chunk_timeseries)
        scale_fac = calc_normalization(chunk_timeseries.chunk_list[t])
        println("# t= ",t, " scale_fac= ", scale_fac)
        normalize_spectrum!(spectra[t], scale_fac)
    end
    return timeseries
end


# Unorganized code
""" Return true if all elements of array are equal to each other. """
@inline function allequal(x::AbstractArray{T,1}) where {T<:Real}
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end
