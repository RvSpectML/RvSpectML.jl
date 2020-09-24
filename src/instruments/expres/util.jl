
""" Normalize spectrum based on continuum model from FITS file. """
function continuum_normalize_spectrum!(spectrum::ST) where { ST<:AbstractSpectra }
    @assert haskey(spectrum.metadata,:continuum)
    spectrum.flux ./= spectrum.metadata[:continuum] .* spectrum.metadata[:blaze]
    spectrum.var ./= (spectrum.metadata[:continuum].* spectrum.metadata[:blaze] ) .^2
    return spectrum
end



""" Normalize each spectrum based on continuum model from FITS files. """
function continuum_normalize_spectra!(spectra::AS) where { ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    for spectrum in spectra
        continuum_normalize_spectrum!(spectrum)
    end
    return spectra
end

function filter_line_list(df::DataFrame, inst::IT ; λmin::Real = default_filter_line_list_λmin, λmax::Real = default_filter_line_list_λmax ) where { # data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:EXPRES.AnyEXPRES } #, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   df |> @filter(λmin <= _.lambda <= λmax) |>
    #    @filter( _.lambda < 6000.0 ) |>                       # Avoid tellurics at redder wavelengths
    #    @filter( _.lambda >6157 || _.lambda < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame
end

function find_worst_telluric_in_each_chunk( clt::AbstractChunkListTimeseries, data::AbstractArray{AS,1} )  where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, IT<:AnyEXPRES, AS<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}  }
   num_lines = num_chunks(clt)
   num_obs = length(clt)
   min_telluric_model_one_obs = ones(num_lines, num_obs )
   min_telluric_model_all_obs  = ones(num_lines)
   for ch_idx in 1:num_chunks(clt)
       for t_idx in 1:num_obs
          view_indices = clt.chunk_list[t_idx].data[ch_idx].λ.indices
          cols = view_indices[1]
          order = view_indices[2]
          min_telluric_model_one_obs[ch_idx, t_idx] = minimum(view(data[t_idx].metadata[:tellurics], cols, order))
      end # times
      min_telluric_model_all_obs[ch_idx] = minimum(min_telluric_model_one_obs[ch_idx, :])
  end # lines
  return min_telluric_model_all_obs
end

default_min_Δv_clean = 8000.0

function find_ranges_with_tellurics_in_order(spectrum::ST, order::Integer; telluric_threshold::Real = 1, min_Δv_clean::Real = default_min_Δv_clean) where { ST<:AbstractSpectra }
    @assert haskey(spectrum.metadata,:tellurics)
    telluric_ranges = Vector{Tuple{Int64,Int64}}(undef,0)

    c_mps = 3e8
    all_wavelengths = range(first(spectrum.λ[:,order]),stop=last(spectrum.λ[:,order]),length=size(spectrum.λ,1) )
    tellurics = spectrum.metadata[:tellurics]
    start = findfirst(x->!isnan(x),view(tellurics,:,order) )
    stop = findlast(x->!isnan(x),view(tellurics,:,order) )
    if isnothing(start) || isnothing(stop)
        # println("# The entire order was NaNs!?!")
        return Vector{Tuple{eltype(spectrum.λ), eltype(spectrum.λ)} }(undef,0)
    end
    @assert stop >= start +1
    in_telluric = tellurics[start,order] < telluric_threshold
    idx_start_this_telluric = in_telluric ? start : 0
    idx_stop_this_telluric = 0
    for i in start+1:stop
        if in_telluric
            if !(tellurics[i,order] < telluric_threshold)  # exitted telluric
                idx_stop_this_telluric = i-1
                #println("# Order = " , order, "  Adding pixels = ", idx_start_this_telluric, " - ", idx_stop_this_telluric, " λ = ", spectrum.λ[idx_start_this_telluric,order], " - ", spectrum.λ[idx_stop_this_telluric,order] )
                push!(telluric_ranges, (idx_start_this_telluric,idx_stop_this_telluric) )
                #push!(telluric_ranges,(spectrum.λ[idx_start_this_telluric,order], spectrum.λ[idx_stop_this_telluric,order]) )
                in_telluric = false
            end
        else
            if (tellurics[i,order] < telluric_threshold)  # entered telluric
                idx_start_this_telluric = i
                if length(telluric_ranges) >= 1
                    idx_last_start = last(telluric_ranges)[1]
                    idx_last_stop = last(telluric_ranges)[2]
                    if spectrum.λ[idx_start_this_telluric,order] - spectrum.λ[idx_last_stop,order] <= min_Δv_clean/c_mps * spectrum.λ[idx_start_this_telluric,order]
                        idx_start_this_telluric = first(last(telluric_ranges))
                        pop!(telluric_ranges)
                    end
                end
                in_telluric = true
            end
        end
    end
    if in_telluric
        idx_stop_this_telluric = stop
        #println("# Order = " , order, "  Adding pixels = ", idx_start_this_telluric, " - ", idx_stop_this_telluric, " λ = ", spectrum.λ[idx_start_this_telluric,order], " - ", spectrum.λ[idx_stop_this_telluric,order] )
        push!(telluric_ranges,(idx_start_this_telluric,idx_stop_this_telluric) )
    end
    #lambda_range = map(r->(spectrum.λ[first(r),order], spectrum.λ[last(r),order] ), telluric_ranges)
    DataFrame(:lambda_lo=>map(r->spectrum.λ[first(r),order], telluric_ranges), :lambda_hi=>map(r->spectrum.λ[last(r),order], telluric_ranges) )
end


function merge_sorted_wavelength_ranges(list_of_ranges::Vector{Tuple{T,T}}; min_Δv_clean::Real = default_min_Δv_clean) where { T<:Real }
    c_mps = 3e8
    non_overlapping_ranges = typeof(list_of_ranges)(undef,0)
    last_range = list_of_ranges[1]
    for i in 2:length(list_of_ranges)
        this_range = list_of_ranges[i]
        #if first(list_of_ranges[i]) < last_range[i-1])
        #    println("# Need to merge ranges: ", list_of_ranges[i-1], "  and  ", list_of_ranges[i])
        #end
        if first(this_range) - last(last_range) > min_Δv_clean/c_mps*last(last_range)
            push!(non_overlapping_ranges, last_range)
            last_range = this_range
        else
            last_range = (last_range[1], this_range[2])
        end
    end
    return non_overlapping_ranges
end


function merge_sorted_wavelength_ranges(df::DataFrame; min_Δv_clean::Real = default_min_Δv_clean) where { T<:Real }
    @assert hasproperty(df,:lambda_lo)
    @assert hasproperty(df,:lambda_hi)
    c_mps = 3e8
    non_overlapping_ranges = DataFrame(:lambda_lo=>Float64[], :lambda_hi=>Float64[] )
    last_range = df[1,:]
    for i in 2:size(df,1)
        this_range = df[i,:]
        #if first(list_of_ranges[i]) < last_range[i-1])
        #    println("# Need to merge ranges: ", list_of_ranges[i-1], "  and  ", list_of_ranges[i])
        #end
        if this_range[:lambda_lo] - last_range[:lambda_hi] > min_Δv_clean/c_mps*last_range[:lambda_hi]
            push!(non_overlapping_ranges, last_range)
            last_range = this_range
        else
            last_range = Dict(:lambda_lo=>last_range[:lambda_lo], :lambda_hi=>this_range[:lambda_hi])
        end
    end
    return non_overlapping_ranges
end

function find_ranges_with_tellurics(spectrum::ST; min_order::Integer = 1, max_order::Integer = size(spectrum.λ,2), min_Δv_clean::Real = default_min_Δv_clean, telluric_threshold::Real = 1 ) where { ST<:AbstractSpectra }
    @assert haskey(spectrum.metadata,:tellurics)
    min_Δv_clean = 10000.0
    c_mps = 3e8

    list_of_ranges = Vector{Tuple{eltype(spectrum.λ),eltype(spectrum.λ)}}(undef,0)
    for ord in min_order:max_order
        ranges_for_order = find_ranges_with_tellurics_in_order(spectrum,ord, min_Δv_clean=min_Δv_clean, telluric_threshold=telluric_threshold)
        append!(list_of_ranges,ranges_for_order)
    end

    if length(list_of_ranges)<=1   return list_of_ranges   end
    sort!(list_of_ranges, by=x->x[1] )
    non_overlapping_ranges = merge_sorted_wavelength_ranges(list_of_ranges,min_Δv_clean=min_Δv_clean)
    return non_overlapping_ranges
end

#=
function find_ranges_with_tellurics(spectra::AS; min_order::Integer = 1, max_order::Integer = size(first(spectra).λ,2), min_Δv_clean::Real = default_min_Δv_clean, telluric_threshold::Real = 1 ) where { ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    list_of_ranges = mapreduce(s->find_ranges_with_tellurics(s, min_order=min_order, max_order=max_order, min_Δv_clean=min_Δv_clean, telluric_threshold=telluric_threshold), append!, spectra)
    #return list_of_ranges
    if size(list_of_ranges,1)<=1   return list_of_ranges   end
    #sort!(list_of_ranges, by=x->x[1] )
    sort!(list_of_ranges, :lambda_lo )
    non_overlapping_ranges = merge_sorted_wavelength_ranges(list_of_ranges, min_Δv_clean=min_Δv_clean)
    return non_overlapping_ranges
end
=#

function is_in_wavelength_range_list(list::AVT, λ::Real ) where { T<:Real, AVT<:AbstractVector{Tuple{T,T}} }
    idx =  searchsortedfirst(list, λ, lt=(x,y)->x[2]<y)
    return idx>length(list) || !(first(list[idx])<=λ<=last(list[idx])) ?  false : true
end

#=
function is_in_wavelength_range_list(list::AVT, λ::AA ) where { T1<:Real, AVT<:AbstractVector{Tuple{T1,T1}}, T2<:Real, AA<:AbstractArray{T2} }
    idx =  searchsortedfirst(list, λ, lt=(x,y)->x[2]<y)
    return idx>length(list) || !(first(list[idx])<=λ<=last(list[idx])) ?  false : true
end
=#

function is_in_wavelength_range_list(λ::Real; list::DataFrame  )
    @assert hasproperty(list, :lambda_lo)
    @assert hasproperty(list, :lambda_hi)
    idx =  searchsortedfirst(list[:,:lambda_hi], λ)
    return idx>size(list,1) || !(list[idx,:lambda_lo]<=λ<=list[idx,:lambda_hi]) ?  false : true
end

function make_ranges_without_tellurics(telluric_list::AVT) where { T<:Real, AVT<:AbstractVector{Tuple{T,T}} }
    chunk_list = map(i->( last(telluric_list[i]), first(telluric_list[i+1]) ), 1:(length(telluric_list)-1) )
    return chunk_list
end

#=
function make_ranges_without_tellurics(telluric_list::AVT; λ_start::Real, λ_stop::Real ) where { T<:Real, AVT<:AbstractVector{Tuple{T,T}} }
    chunk_list = vcat( (λ_start , telluric_list[1] ),
            map(i->( last(telluric_list[i]), first(telluric_list[i+1]) ), 1:(length(telluric_list)-1) ),
            ( telluric_list[end] , λ_stop) )
end

function make_ranges_without_tellurics(telluric_list::AVT; λ_start::Real ) where { T<:Real, AVT<:AbstractVector{Tuple{T,T}} }
    chunk_list = vcat( (λ_start , telluric_list[1] ),
            map(i->( last(telluric_list[i]), first(telluric_list[i+1]) ), 1:(length(telluric_list)-1) ) )
end

function make_ranges_without_tellurics(telluric_list::AVT; λ_stop::Real ) where { T<:Real, AVT<:AbstractVector{Tuple{T,T}} }
    chunk_list = vcat( map(i->( last(telluric_list[i]), first(telluric_list[i+1]) ), 1:(length(telluric_list)-1) ),
            ( telluric_list[end] , λ_stop) )
end
=#

function make_ranges_without_tellurics(df::DataFrame)
    DataFrame(:lambda_lo => df[1:(size(df,1)-1), :lambda_hi],  :lambda_hi => df[2:size(df,1), :lambda_lo] )
end

#=
function make_ranges_without_tellurics(df::DataFrame; λ_start::Real, λ_stop::Real )
    output = DataFrame(:lambda_lo=>[λ_start],:lambda_hi=>df[1, :lambda_lo])
    append!(output, DataFrame(:lambda_lo => df[1:(size(df,1)-1), :lambda_hi],  :lambda_hi => df[2:size(df,1), :lambda_lo] ) )
    append!(output, Dict(:lambda_lo => df[end, :lambda_hi],  :lambda_hi => λ_stop ) )
    return output
end

function make_ranges_without_tellurics(df::DataFrame; λ_start::Real )
    output = DataFrame(:lambda_lo=>[λ_start],:lambda_hi=>df[1, :lambda_lo])
    append!(output, DataFrame(:lambda_lo => df[1:(size(df,1)-1), :lambda_hi],  :lambda_hi => df[2:size(df,1), :lambda_lo] ) )
    return output
end

function make_ranges_without_tellurics(df::DataFrame; λ_stop::Real )
    output = DataFrame(:lambda_lo => df[1:(size(df,1)-1), :lambda_hi],  :lambda_hi => df[2:size(df,1), :lambda_lo]  )
    append!(output, Dict(:lambda_lo => df[end, :lambda_hi],  :lambda_hi => λ_stop ) )
    return output
end
=#

function break_chunk_into_chunks_without_tellurics(chunk::AbstractChunkOfSpectrum, tellurics::DataFrame)
    telluric_mask = is_in_wavelength_range_list.(chunk.λ, list=tellurics)
    new_chunks = Vector{typeof(chunk)}(undef,0)
    idx_start = findfirst(.!telluric_mask)
    idx_stop = findlast(.!telluric_mask)
    if isnothing(idx_start) || isnothing(idx_stop)
        return new_chunks
    end
    idx_lo = idx_start
    idx_hi = idx_start
    while idx_hi < idx_stop
        idx_lo = findfirst( view(telluric_mask,idx_hi:idx_stop) )
        if isnothing(idx_lo)
            break
        end
        idx_lo += idx_lo-1
        idx_hi = findfirst( view(telluric_mask,idx_lo:idx_stop) )
        if isnothing(idx_hi)
            idx_hi = idx_stop
        else
            idx_hi += idx_lo-2
        end
        if idx_hi-idx_lo > 10 # TODO Make min_chunk_length
            push!(new_chunks,ChunkOfSpectrum())
        end
    end
    return new_chunks
end
