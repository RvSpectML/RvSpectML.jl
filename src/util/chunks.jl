"""
Code for creating and manipulating chunklists (i.e., list of views into a spectrum).
For example, creating a list of views into the orders of a spectrum to be analyzed or
creating a list of views into chunks around entries in a line list.

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

""" make_orders_into_chunks

Return a ChunkList with a region of spectrum from each order in orders_to_use.
# Arguments
- spectra<:AbstractSpectra
- inst:  Instrument trait that provides default values
# Optional arguments
- orders_to_use: Range or Array (orders_to_use(inst))
- pixels_to_use: Array of Ranges (each from min_col to max_col)
or
- min_col: (min_col_default(inst,order)) and
- max_col: (max_col_default(inst,order))
"""
function make_orders_into_chunks
end

function make_orders_into_chunks(spectra::AS, inst::AbstractInstrument;
        orders_to_use = orders_to_use_default(spectra.inst) # 1:size(spectra.flux,2),
        ) where {AS<:AbstractSpectra }
    @assert eltype(orders_to_use) <: Integer
    @assert all( min_order(spectra.inst) .<= orders_to_use .<= max_order(spectra.inst) )
    pixels_to_use = map(ord->min_col_default(spectra.inst,ord):max_col_default(spectra.inst,ord),orders_to_use)
    make_orders_into_chunks(spectra,orders_to_use=orders_to_use, pixels_to_use=pixels_to_use)
end

function make_orders_into_chunks(spectra::AS;
        orders_to_use::Union{AR,AA1}, pixels_to_use::AAR ) where {
         AS<:AbstractSpectra, AR<:AbstractRange{Int64}, AA1<:AbstractArray{Int64,1}, AAR<:AbstractArray{AR,1} }
    @assert eltype(orders_to_use) <: Integer
    #@assert all( min_order(inst) .<= orders_to_use .<= max_order(inst) )
    #@assert minimum(pixels_to_use) >= min_pixel_in_order(inst)
    #@assert maximum(pixels_to_use) <= max_pixel_in_order(inst)
    ChunkList( map(order_idx->
                    ChunkOfSpectrum(spectra,(pixels=pixels_to_use[order_idx],order=orders_to_use[order_idx])),
                    1:length(orders_to_use) ))
end


""" make_grid_for_chunk
Create a range with equal spacing between points with end points set based on union of all chunks in timeseries.
# Arguments:
- timeseries: ChunkListTimeseries
- chunk index:
- oversample_factor: (1)
"""
function make_grid_for_chunk(timeseries::ACLT, c::Integer; oversample_factor::Real = 1.0, spacing::Symbol = :Linear, remove_rv_est::Bool = false ) where { ACLT<:AbstractChunkListTimeseries }
    num_obs = length(timeseries.chunk_list)
    @assert num_obs >= 1
    @assert 1<= c <= length(first(timeseries.chunk_list).data)
    @assert allequal(map(chunk->length(chunk.data),timeseries.chunk_list))
    @assert spacing == :Log || spacing == :Linear
    if spacing == :Log
        @warn "There's some issues with end points exceeding the bounds.  Round off error?  May cause bounds errors."
    end
    # Create grid, so that chunks at all times include the grid's minimum and maximum wavelength.
    if remove_rv_est   @assert haskey(first(timeseries.metadata),:rv_est)   end
    boost_factor = [ remove_rv_est ? calc_doppler_factor(timeseries.metadata[t][:rv_est]) : 1 for t in 1:num_obs ]
    λ_min = maximum(minimum(timeseries.chunk_list[t].data[c].λ)/boost_factor[t] for t in 1:num_obs)
    λ_max = minimum(maximum(timeseries.chunk_list[t].data[c].λ)/boost_factor[t] for t in 1:num_obs)
    Δλ_grid_obs = mean(log(timeseries.chunk_list[t].data[c].λ[end]/
                           timeseries.chunk_list[t].data[c].λ[1]   )/
                         (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
    num_pixels_obs = log(λ_max/λ_min)/Δλ_grid_obs
    num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
    if spacing == :Log
        Δlnλ_grid_obs = mean(log(timeseries.chunk_list[t].data[c].λ[end]/
                                 timeseries.chunk_list[t].data[c].λ[1]   )/
                                 (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
        num_pixels_obs = log(λ_max/λ_min)/Δlnλ_grid_obs
        num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
            Δlnλ_grid_gen = log(λ_max/λ_min)/ (num_pixels_gen-1)
        return exp.(range(log(λ_min),stop=log(λ_max),step=Δlnλ_grid_gen))
    elseif spacing == :Linear
        Δλ_grid_obs = mean((timeseries.chunk_list[t].data[c].λ[end]-timeseries.chunk_list[t].data[c].λ[1] )/
                            (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
        num_pixels_obs = (λ_max-λ_min)/Δλ_grid_obs
        num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
        Δλ_grid_gen = (λ_max-λ_min)/ (num_pixels_gen-1)
        return range(λ_min,stop=λ_max,step=Δλ_grid_gen)
    end
end

""" Return a ChunkList of best regions of spectrum with lines in line_line.
    line_list is a DataFrame containing :lambda_lo and :lambda_hi.
    Pads edges by Δ.
"""
function make_chunk_list(spectra::AS, line_list::DataFrame; rv_shift::Real = 0, Δ::Real=Δλoλ_edge_pad_default) where { AS<:AbstractSpectra }
    @assert hasproperty(line_list,:lambda_lo)
    @assert hasproperty(line_list,:lambda_hi)
    boost_factor = calc_doppler_factor(rv_shift)
    if rv_shift != 0
        @warn("I haven't tested this yet, especially the sign.")  # TODO
    end
    cl = ChunkList(map(row->ChunkOfSpectrum(spectra,find_line_best(row.lambda_lo*boost_factor,row.lambda_hi*boost_factor,spectra,Δ=Δ)), eachrow(line_list) ))

end

function make_chunk_list_timeseries(spectra::AS,chunk_list_df::DataFrame; rv_shift::Real = 0) where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    times = map(s->s.metadata[:bjd],spectra)
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    time_series_of_chunk_lists = map(spec->RvSpectML.make_chunk_list(spec,chunk_list_df, rv_shift=rv_shift),spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, time_series_of_chunk_lists, inst=first(spectra).inst, metadata=metadata )
end

function extract_chunk_list_timeseries_for_order(clt::AbstractChunkListTimeseries, order::Integer) # where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    @assert min_order(clt.inst) <= order <= max_order(clt.inst)
    num_chunks_to_search = num_chunks(clt)
    chunk_order = map(ch->first(clt.chunk_list).data[ch].λ.indices[2], 1:num_chunks_to_search )
    chunks_in_order = findall(chunk_order .== order)
    num_chunks_in_order = sum(chunks_in_order)
    if !(num_chunks_in_order >= 1)
        return nothing
    end
    new_chunk_list_timeseries = Vector{ChunkList}(undef, length(clt.times))
    for i in 1:length(clt.times)
        chunk_list = ChunkList(clt.chunk_list[i].data[chunks_in_order])
        new_chunk_list_timeseries[i] = chunk_list
    end
    output = ChunkListTimeseries(clt.times, new_chunk_list_timeseries, inst=clt.inst, metadata=clt.metadata )
    return output
end

#=
function make_order_list_timeseries(spectra::AS) #= , order_list::AOL ) =# where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} #=, CLT<:AbstractChunkList, AOL::AbstractArray{CLT,1} =# }
    times = map(s->s.metadata[:bjd],spectra)
    inst = first(spectra).inst
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    order_list = map( spec->RvSpectML.make_orders_into_chunks(spec,inst), spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, order_list, inst=first(spectra).inst, metadata=metadata )
end
=#

function make_order_list_timeseries(spectra::AS; orders_to_use = orders_to_use_default(spectra.inst) ) #= , order_list::AOL ) =# where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} #=, CLT<:AbstractChunkList, AOL::AbstractArray{CLT,1} =# }
    times = map(s->s.metadata[:bjd],spectra)
    inst = first(spectra).inst
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    order_list = map( spec->RvSpectML.make_orders_into_chunks(spec,inst, orders_to_use=orders_to_use), spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, order_list, inst=first(spectra).inst, metadata=metadata )
end


function make_chunk_list_expr(spectra::AS, range_list::DataFrame ) where { AS<:AbstractSpectra }
    @assert hasproperty(range_list,:lambda_lo)
    @assert hasproperty(range_list,:lambda_hi)
    cl = typeof(ChunkOfSpectrum(spectra,1,1:size(spectra.flux,1)))[]
    for row in eachrow(range_list)
        orders = find_orders_in_range(row.lambda_lo,row.lambda_hi, spectra.λ)
        if length(orders) == 1
            order = orders[1]
            pixels = find_cols_to_fit(spectra.λ[:,order],row.lambda_lo,row.lambda_hi)
            if 1 <= first(pixels) <= size(spectra.λ,1) && 1 <= last(pixels) <= size(spectra.λ,1) &&  # valid location
                calc_snr(spectra.flux[pixels,order],spectra.var[pixels,order]) > 0     # not just NaN's
                println("# Found one order (", order, ") for range ", row.lambda_lo, " - ", row.lambda_hi )
                push!(cl, ChunkOfSpectrum(spectra,order,pixels) )
            else
                println("# Found one order (", order, ") for range ", row.lambda_lo, " - ", row.lambda_hi, " but not valid or all NaNs.")
                continue
            end
        elseif length(orders) > 1
            pixels = map(ord->find_cols_to_fit(spectra.λ[:,ord],row.lambda_lo,row.lambda_hi), orders)
            scores = map( i->calc_snr(spectra.flux[pixels[i],orders[i]],spectra.var[pixels[i],orders[i]]), 1:length(orders) )
            idx_best = findmax(scores)[2]
            println("# Found ", length(orders), " orders for range ", row.lambda_lo, " - ", row.lambda_hi, " picking ", orders[idx_best])
            push!(cl, ChunkOfSpectrum(spectra,orders[idx_best],pixels[idx_best]) )
            idx_all = findall(scores)

        else # length(orders) == 0
            #println("# Found no orders for range ", row.lambda_lo, " - ", row.lambda_hi )

        end
    end
    return ChunkList(cl)
end

""" Return list of all orders that contain a pixel with wavelength lambda """
function find_orders_with_line(goal::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   for i in 1:5
       println("# i= ",i," order_min= ",order_min(i)," order_max= ",order_max(i), "   goal= ",goal)
   end
   flush(stdout)
   findall(i->order_min(i)<=goal<=order_max(i), 1:size(lambda,2) )
end

""" Return list of all orders that contain all pixels with wavelengths between goal_lo and goal_hi """
function find_orders_with_line(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal_lo && goal_hi<=order_max(i), 1:size(lambda,2) )
end

""" Return list of all orders that include any wavelengths between goal_lo and goal_hi """
function find_orders_in_range(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i-> (goal_lo<=order_min(i)<=goal_hi) || (goal_lo<=order_max(i)<=goal_hi), 1:size(lambda,2) )
end


# Find indicies for pixels around lines
const Δλoλ_fit_line_default = 5*(1.8*1000/speed_of_light_mps)
const Δλoλ_edge_pad_default = 0*(1.8*1000/speed_of_light_mps)

""" Return a range of columns indices with wavelengths within Δ of line_center """
function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_center::Real; Δ::Real = Δλoλ_fit_line_default) where T<:Real
    @assert Δ >= zero(Δ)
    first = findfirst(x->x>=line_center*(1-Δ),wavelengths)
    last = findlast(x->x<=line_center*(1+Δ),wavelengths)
    if isnothing(first) || isnothing(last)   return 0:0   end
    return first:last
end

""" Return a range of columns indices with wavelengths between line_lo and line_hi """
function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_lo::Real, line_hi::Real; Δ::Real = Δλoλ_edge_pad_default) where T<:Real
    @assert line_lo < line_hi
    first = findfirst(x->x>=line_lo*(1-Δ),wavelengths)
    last = findlast(x->x<=line_hi*(1+Δ),wavelengths)
    if isnothing(first) || isnothing(last)   return 0:0   end
    return first:last
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
    #locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_edge_pad_default, verbose::Bool = false) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal_lo < goal_hi <= lambda[end,end]
    orders = find_orders_with_line(goal_lo,goal_hi,lambda)
    #if ! (length(orders) >= 1) return end
    if verbose
        for i in 1:5
        println("# i= ",i," min(order)= ",minimum(lambda[:,i])," max(order)= ",maximum(lambda[:,i]), "   goal_lo= ",goal_lo, " goal_hi = ",goal_hi)
        end
    end
    flush(stdout)
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal_lo, goal_hi,Δ=Δ),order=o), orders)
    #locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal::Real,lambda::AbstractArray{T1,1},var::AbstractArray{T2,1}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real}
    @assert lambda[1] <= goal <= lambda[end]
    @assert size(lambda) == size(var)
    @assert Δ >= zero(Δ)
    locs = find_cols_to_fit(lambda,goal,Δ=Δ)
    locs_good_idx = findall(t->!any(isnan.(var[t])),locs)
    #locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,1},var::AbstractArray{T2,1}; Δ::Real = Δλoλ_edge_pad_default, verbose::Bool = false) where {T1<:Real, T2<:Real}
    @assert lambda[1] <= goal_lo < goal_hi <= lambda[end]
#=    if verbose
        for i in 1:5
        println("# i= ",i," min(order)= ",minimum(lambda[:,i])," max(order)= ",maximum(lambda[:,i]), "   goal_lo= ",goal_lo, " goal_hi = ",goal_hi)
        end
    end
    flush(stdout)
    =#
    locs = find_cols_to_fit(lambda,goal_lo, goal_hi,Δ=Δ)
    locs_good_idx = findall(t-> !any(isnan.(var[t])) ,locs)
    #locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
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
    if length(locs) == 0
        println("=>(",goal_lo, ", ",goal_hi, ")  Δ=",Δ)
        return  missing
    end
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal::Real,lambda::AbstractArray{T1,1},flux::AbstractArray{T2,1},var::AbstractArray{T3,1}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real, T3<:Real}
    cols = find_cols_to_fit(lambda,goal, Δ=Δ)
    @assert( ( first(cols)==0 && last(cols)==0)  || !any(isnan.(var[cols])) )
    return cols
    #=
    locs = findall_line(goal,lambda,var,Δ=Δ)
    if length(locs) == 0   return  missing end
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    return locs
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
    =#
end

function find_line_best(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,1},flux::AbstractArray{T2,1},var::AbstractArray{T3,1}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real, T3<:Real}
    cols = find_cols_to_fit(lambda,goal_lo, goal_hi, Δ=Δ)
    @assert( ( first(cols)==0 && last(cols)==0)  || !any(isnan.(var[cols])) )
    return cols
    #=
    locs = findall_line(goal_lo,goal_hi,lambda,var,Δ=Δ)
    if length(locs) == 0
        println("=>(",goal_lo, ", ",goal_hi, ")  Δ=",Δ)
        return  missing
    end
    return locs
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t],var[t]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
    =#
end

function find_line_best(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra}
    find_line_best(goal,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end

function find_line_best(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra}
    find_line_best(goal_lo,goal_hi,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end

# Manipulation data in chunks of spectrua

""" Return (chunk_timeseries, line_list) that have been trimmed of any chunks that are bad based on any spectra in the chunk_timeseries.
    chunk_timeseries: ChunkListTimeseries
    line_linst:  DataFrame w/ lambda_lo, lambda_hi
    verbose: print debugging info (false)
"""
#= Is there any reason to keep this version?
function filter_bad_chunks(chunk_list_timeseries::ACLT, line_list::DataFrame; verbose::Union{Int,Bool} = false) where { ACLT<:AbstractChunkListTimeseries }
    @assert(length(chunk_list_timeseries)>=1)
    @assert(hasproperty(line_list,:lambda_lo))
    @assert(hasproperty(line_list,:lambda_hi))
    inst = chunk_list_timeseries.inst
    idx_keep = trues(num_chunks(chunk_list_timeseries))
    for t in 1:length(chunk_list_timeseries)
        idx_bad_λ = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].λ)),1:num_chunks(chunk_list_timeseries))
        idx_bad_flux = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].flux)),1:num_chunks(chunk_list_timeseries))
        idx_bad_var = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].var)),1:num_chunks(chunk_list_timeseries))
        idx_not_sorted = findall(c->!issorted(chunk_list_timeseries.chunk_list[t].data[c].λ),1:num_chunks(chunk_list_timeseries))
        idx_keep[idx_bad_λ] .= false
        idx_keep[idx_bad_flux] .= false
        idx_keep[idx_bad_var] .= false
        idx_keep[idx_not_sorted] .= false
        if verbose && (length(idx_bad_λ)+length(idx_bad_flux)+length(idx_bad_var)+length(idx_not_sorted))
               flush(stdout)
               println("# Removing chunks", vcat(idx_bad_λ,idx_bad_flux,idx_bad_var), " at time ", t, " due to NaNs (",
                            length(idx_bad_λ),",",length(idx_bad_flux),",",length(idx_bad_var),") and ",
                            length(idx_not_sorted), " for not being sorted.")
        end
        #=
        # TODO: Move these kinds of checks to a traits/plan-based system
        if hasproperty(chunk_list_timeseries.metadata[t],:pixel_mask)
            idx_keep[.!chunk_list_timeseries.metadata[t][:pixel_mask]] .= false
        end
        if hasproperty(chunk_list_timeseries.metadata[t],:excalibur_mask)
            idx_keep[.!chunk_list_timeseries.metadata[t][:excalibur_mask]] .= false
        end
        =#
    end
    chunks_to_remove = findall(.!idx_keep)
    if length(chunks_to_remove) == 0
        println("# No lines to remove.")
        return (chunk_timeseries=chunk_list_timeseries, line_list=line_list)
    else
        println("# Removing ", length(chunks_to_remove), " chunks.")
        map(c->println("# ",c,": ",line_list.lambda_lo[c]," - ",line_list.lambda_hi[c]),chunks_to_remove)
        new_line_list = line_list[findall(idx_keep),:]
        new_chunk_list_timeseries = [ChunkList(chunk_list_timeseries.chunk_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        return (chunk_timeseries=ChunkListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries, inst=inst, metadata=chunk_list_timeseries.metadata), line_list=new_line_list)
    end
end
=#
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
            map(c->println("# ",c,": ",findall(.!idx_keep)))
        end
        new_chunk_list_timeseries = [ChunkList(chunk_list_timeseries.chunk_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        #return ChunkListTimeseries(chunk_list_timeseries.times[idx_keep],new_chunk_list_timeseries, inst=chunk_list_timeseries.inst, metadata=chunk_list_timeseries.metadata[idx_keep])
        return ChunkListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries, inst=chunk_list_timeseries.inst, metadata=chunk_list_timeseries.metadata)
    end
end

""" Find pixels included in a range of wavelengths """
function find_pixels_for_line_in_chunk( chunk::AbstractChunkOfSpectrum, λ_min::Real, λ_max::Real )# ; plan::LineFinderPlan = LineFinderPlan() )
  idx_lo = searchsortedfirst(chunk.λ, λ_min, by=x->x>=λ_min)
  idx_tmp = searchsortedlast(chunk.λ[idx_lo:end], λ_max, by=x->x<=λ_max, rev=true)
  idx_hi = idx_lo + idx_tmp - 1
  return idx_lo:idx_hi
end

function find_pixels_for_line_in_chunklist( chunk_list::AbstractChunkList, λ_min::Real, λ_max::Real; verbose::Bool = true )
  ch_idx_all = findall(c-> (λ_min <= minimum(chunk_list.data[c].λ)) && (maximum(chunk_list.data[c].λ) <= λ_max) ,1:length(chunk_list))
  println("Hello")
  #map(c->(chunk_idx=c, pixels=find_pixels_for_line_in_chunk(chunk_list.data[c], λ_min, λ_max) ), ch_idx)
  ch_idx = 0
  if length(ch_idx_all) > 1
    snr_of_chunks_with_line = map(c->RvSpectML.calc_snr(chunk_list.data[c].flux, chunk_list.data[c].var), ch_idx_all)
    ch_idx_to_keep = argmax(snr_of_chunks_with_line)
    ch_idx = ch_idx_all[ch_idx_to_keep]
    if verbose
      println(" Found λ=",λ_min,"-",λ_max," in chunks: ", ch_idx_all, " containing ", length.(ch_idx_all), " pixels.")
      println(" SNRs = ", snr_of_chunks_with_line)
      println(" Keeping chunk #",ch_idx)
    end
  elseif length(ch_idx_all) == 1
    ch_idx = first(ch_idx_all)
    if verbose
      println(" Found λ=",λ_min,"-",λ_max," in chunk: ", ch_idx, " containing ", length(ch_idx), " pixels.")
      snr_of_chunk_with_line = RvSpectML.calc_snr(chunk_list.data[ch_idx].flux, chunk_list.data[ch_idx].var)
      println(" SNRs = ", snr_of_chunk_with_line)
    end

  end
  if ch_idx == 0
    error("Didn't find λ = " *string(λ_min)*" - " *string(λ_max)* " in chunklist.")

  end
  return (chunk_idx=ch_idx, pixels=find_pixels_for_line_in_chunk(chunk_list.data[ch_idx], λ_min, λ_max) )
end

function find_pixels_for_line_in_chunklist( chunk_list::AbstractChunkList, λ_min::Real, λ_max::Real, chunk_id::Integer)
  return (chunk_idx=chunk_id, pixels=find_pixels_for_line_in_chunk(chunk_list.data[chunk_id], λ_min, λ_max) )
end



""" Calc normalization of chunk based on average flux in a ChunkOfSpectrum. """
function calc_normalization(chunk::AC) where { AC<:AbstractChunkOfSpectrum}
    total_flux = NaNMath.sum(Float64.(chunk.flux))
    num_pixels = length(flux)
    scale_fac = num_pixels / total_flux
end

""" Calc normalization of chunk based on average flux in a ChunkOfSpectrum using inverse variance weighting. """
function calc_normalization_var_weighted(chunk::AC) where { AC<:AbstractChunkOfSpectrum}
    sum_weighted_flux = NaNMath.sum(Float64.(chunk.flux) ./ chunk.var )
    sum_weights = NaNMath.sum(1.0 ./ chunk.var)
    scale_fac = sum_weights / sum_weighted_flux
end

""" Calc normalization of spectra based on average flux in a ChunkList. """
function calc_normalization(chunk_list::ACL) where { ACL<:AbstractChunkList}
    total_flux = NaNMath.sum(NaNMath.sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end

""" Calc normalization of spectra based on average flux in a ChunkList using inverse variance weighting. """
function calc_normalization_var_weighted(chunk_list::ACL) where { ACL<:AbstractChunkList}
    total_flux = NaNMath.sum(NaNMath.sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end
