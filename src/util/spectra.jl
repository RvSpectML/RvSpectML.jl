"""
Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

"""
    apply_doppler_boost!(spectrum, doppler_factor) -> typeof(spectrum)
    apply_doppler_boost!(spectra, df) -> typeof(spectra)

Apply Doppler boost to spectra's λ's and update its metadata[:doppler_factor], so it will know how to undo the transform.
# Arguments:
* `spectrum::AbstractSpectra`: spectrum to be boosted
* `doppler_factor::Real`: boost factor (1 = noop)
or
* `spectra::AbstractArray{<:AbstractSpectra}`: spectra to be boosted
* `df::DataFrame`: provides `:drift` and `:ssb_rv` (in m/s) for calculating the Doppler boost for each spectrum

TODO: Improve documentation formatting.  This can serve as a template.
"""
function apply_doppler_boost! end

function apply_doppler_boost!(spectra::AS,doppler_factor::Real) where {AS<:AbstractSpectra}
    if doppler_factor == one(doppler_factor) return spectra end
    #println("# t= ",time, " doppler_factor= ",doppler_factor)
    spectra.λ .*= doppler_factor
    if hasproperty(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor
    end
    return spectra
end

function apply_doppler_boost!(spectra::AbstractArray{AS}, df::DataFrame ) where { AS<:AbstractSpectra }
    @assert size(spectra,1) == size(df,1)
    local doppler_factor = ones(size(spectra))
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :drift to apply."   end
    if  hasproperty(df,:drift)        doppler_factor .*= calc_doppler_factor.(df[!,:drift])          end
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :ssb_rv to apply."  end
    if  hasproperty(df,:ssb_rv)       doppler_factor   .*= calc_doppler_factor.(df[!,:ssb_rv])       end
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :diff_ext_rv to apply."  end
    if  hasproperty(df,:diff_ext_rv)  doppler_factor   .*= calc_doppler_factor.(df[!,:diff_ext_rv])  end
    map(x->apply_doppler_boost!(x[1],x[2]), zip(spectra,doppler_factor) );
end

"""  Calculate total SNR in (region of) spectra. """
function calc_snr(flux::AbstractArray{T1},var::AbstractArray{T2}) where {T1<:Real, T2<:Real}
    @assert size(flux) == size(var)
    sqrt(sum(flux./var))   # TODO: Generalize & Test when weights are far from equal
end

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
        # println("# t= ",t, " scale_fac= ", scale_fac)
        normalize_spectrum!(spectra[t], scale_fac)
    end
    return chunk_timeseries
end


function bin_consecutive_spectra(spectra::AbstractSpectralTimeSeriesCommonWavelengths, n::Integer)
  local num_in = size(spectra.flux,2)
  @assert num_in >= n
  local num_out = floor(Int,num_in//n)
  local flux_out = Array{eltype(spectra.flux),2}(undef,size(spectra.flux,1),num_out)
  local var_out  = Array{eltype(spectra.flux),2}(undef,size(spectra.flux,1),num_out)
  local time_idx_out  = Array{UnitRange,1}(undef,num_out)
  #idx_binned = map(i->1+(i-1)*obs_per_bin:obs_per_bin*i,1:num_obs_binned)

  local idx_start = 1
  for i in 1:num_out
    idx_stop = min(idx_start+n-1,num_in)
    flux_out[:,i] .= vec(sum(spectra.flux[:,idx_start:idx_stop],dims=2))
    var_out[:,i]  .= vec(sum((spectra.var[:,idx_start:idx_stop]),dims=2))
    time_idx_out[i] = idx_start:idx_stop
    idx_start = idx_stop + 1
  end
  metadata_out = MetadataT(:time_idx=>time_idx_out)
  return typeof(spectra)(spectra.λ,flux_out,var_out,spectra.chunk_map,spectra.inst,metadata_out)
end

function bin_times(spectra::AbstractSpectralTimeSeriesCommonWavelengths, times::AT, n::Integer) where { T<:Real, AT<:AbstractVector{T} }
    @assert haskey(spectra.metadata,:time_idx)
    map(idx->mean(times[idx]), spectra.metadata[:time_idx])
end

function bin_times(times::AT, n::Integer) where { T<:Real, AT<:AbstractVector{T} }
    local num_in = length(times)
    @assert num_in >= n
    local num_out = floor(Int,num_in//n)
    local times_out = Vector{eltype(times)}(undef,num_out)
    local time_idx_out  = Array{UnitRange,1}(undef,num_out)
    local idx_start = 1
    for i in 1:num_out
      idx_stop = min(idx_start+n-1,num_in)
      times_out[i] .= mean(times[idx_start:idx_stop])
      time_idx_out[i] = idx_start:idx_stop
      idx_start = idx_stop + 1
    end
    return times_out
end



#=
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
    if length(locs) == 0   return  missing end
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal_lo,goal_hi,lambda,var,Δ=Δ)
    if length(locs) == 0   return  missing end
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
    @assert hasproperty(line_list,:lambda_lo)
    @assert hasproperty(line_list,:lambda_hi)
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

function make_orders_into_chunks(spectra::AS, inst::AbstractInstrument;
        min_col::Integer = min_col_default(spectra.inst), # min_col_neid_default,
        max_col::Integer = max_col_default(spectra.inst), #max_col_neid_default ,
        orders_to_use = orders_to_use_default(spectra.inst) # 1:size(spectra.flux,2),
        ) where {AS<:AbstractSpectra }
    @assert eltype(orders_to_use) <: Integer
    @assert all( min_order(spectra.inst) .<= orders_to_use .<= max_order(spectra.inst) )
    @assert min_col >= min_pixel_in_order(spectra.inst)
    @assert max_col <= max_pixel_in_order(spectra.inst)
    pixels_to_use=fill(min_col:max_col,length(orders_to_use))
    make_orders_into_chunks(spectra,orders_to_use=orders_to_use, pixels_to_use=pixels_to_use)
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
    @assert(hasproperty(line_list,:lambda_lo))
    @assert(hasproperty(line_list,:lambda_hi))
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

function filter_bad_chunks(chunk_list_timeseries::ACLT; verbose::Bool = false) where { ACLT<:AbstractChunkListTimeseries }
    @assert(length(chunk_list_timeseries)>=1)
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
        if verbose   println("# Nothing to remove.")   end
        return chunk_list_timeseries
    else
        if verbose
            println("# Removing ", length(chunks_to_remove), " chunks due to NaNs.")
            map(c->println("# ",c,": ",line_list.lambda_lo[c]," - ",line_list.lambda_hi[c]),chunks_to_remove)
        end
        new_chunk_list_timeseries = [ChunkList(chunk_list_timeseries.chunk_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        return ChunkListTimeseries(chunk_list_timeseries.times[idx_keep],new_chunk_list_timeseries, inst=chunk_list_timeseries.inst, metadata=chunk_list_timeseries.metadata[idx_keep])
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
    # Create grid, so that chunks at all times include the grid's minimum and maximum wavelength.
    λ_min = maximum(minimum(timeseries.chunk_list[t].data[c].λ) for t in 1:num_obs)
    λ_max = minimum(maximum(timeseries.chunk_list[t].data[c].λ) for t in 1:num_obs)
    Δλ_grid_obs = median((timeseries.chunk_list[t].data[c].λ[end]-
            timeseries.chunk_list[t].data[c].λ[1])/
              (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
    num_pixels_obs = mean(length(timeseries.chunk_list[t].data[c].λ) for t in 1:num_obs)
    num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
    Δλ_grid_gen = (λ_max-λ_min)/ (num_pixels_gen-1)
    return range(λ_min,stop=λ_max,step=Δλ_grid_gen)
    # TODO: Tweak to make other code work with grids spaced eveningly in logλ
    #Δlnλ_grid_gen = log(λ_max/λ_min)/ (num_pixels_gen-1)
    #exp.(range(log(λ_min),stop=log(λ_max),step=Δlnλ_grid_gen))
end

=#
