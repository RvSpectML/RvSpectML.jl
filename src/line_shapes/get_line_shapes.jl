
"""
Code to measure line bisectors and widths as a function of depth.
Author: Alex Wise
Created: September 2020
"""

""" get_line_shapes(line_list, clt, depths, v_width)
For each spectrum in clt (chunklist timerseries), returns data frames for
line bisector and width as a function of fractional line depth given in frac_depths
note: frac_depths must be sorted from largest to smallest
"""
function get_line_shapes(line_list::A1, clt::AbstractChunkListTimeseries;
                    frac_depths::A2 = 0.99:-0.01:0.05, v_width::Real=15e3,
                    RV_estimate::Real=-5000.0) where { T1<:Real, A1<:AbstractArray{T1,1},
                    T2<:Real, A2<:AbstractArray{T2,1}}
    ll = copy(line_list)
    n_lines = length(ll)
    n_spectra = length(clt.chunk_list)
    n_chunks = length(clt.chunk_list[1].data) # this assumes the number of chunks in every spectrum is the same
    n_depths = length(frac_depths)
    #allocate memory for line shapes
    bisectors = zeros(n_lines,n_depths,n_spectra)
    widths = zeros(n_lines,n_depths,n_spectra)
    depths = zeros(n_lines,n_spectra)
    #shift each entry in line list to RV estimate of star
    ll .*= calc_doppler_factor(RV_estimate)
    #find center of each chunk so that we can assign each line in line list to a chunk
    chunk_centers = zeros(n_chunks)
    for i in 1:n_chunks
        chunk_centers[i] = (first(clt.chunk_list[1].data[i].λ) + last(clt.chunk_list[1].data[i].λ)) / 2.0
    end
    #pick which chunk to measure each line shape in
    line_chunk_index = zeros(Int64,n_lines)
    for i in 1:n_lines
        line_chunk_index[i] = findmin(abs.((chunk_centers .- ll[i]) ./ chunk_centers))[2]
    end
    #measure the bisectors and widths
    valid_lines = ones(Bool,n_lines)
    for k in 1:n_spectra
        #println("k = ",k)
        for i in 1:n_lines
            if valid_lines[i]
                #find pixels within v_width of ll[i]
                #println("i = ",i)
                line_high = ll[i] * calc_doppler_factor(v_width)
                line_low = ll[i] * calc_doppler_factor(-v_width)
                index_low = searchsortedfirst(clt.chunk_list[k].data[line_chunk_index[i]].λ,line_low) - 1
                index_high = searchsortedfirst(clt.chunk_list[k].data[line_chunk_index[i]].λ,line_high)
                if (index_low < 1) | (index_high > length(clt.chunk_list[k].data[line_chunk_index[i]].λ))
                    println("line number ",i," is not contained in chunk list ",k,". setting bisector and width functions to zeros for this line in all chunk lists.")
                    valid_lines[i] = false
                    bisectors[i,:,:] .= 0.0
                    widths[i,:,:] .= 0.0
                else
                    #get λs and fluxes for these pixels
                    local_λ = view(clt.chunk_list[k].data[line_chunk_index[i]].λ,index_low:index_high)
                    local_flux = view(clt.chunk_list[k].data[line_chunk_index[i]].flux,index_low:index_high)
                    depths[i,k] = 1.0 - minimum(local_flux[ Int(floor(length(local_flux) / 4)) : Int(ceil(length(local_flux) * 3 / 4)) ]) / maximum(local_flux)
                    for j in 1:n_depths
                        b = RvSpectMLBase.calc_line_bisector_at_frac_depth(local_λ, local_flux, frac_depth=frac_depths[j])
                        w = RvSpectMLBase.calc_line_width(local_λ, local_flux, frac_depth=frac_depths[j])
                        if isnan(b)
                            bisectors[i,j,k] = bisectors[i,j-1,k]
                        else
                            bisectors[i,j,k] = b
                        end
                        if isnan(w)
                            widths[i,j,k] = widths[i,j-1,k] #TODO: figure out what to do here to make more sense
                        else
                            widths[i,j,k] = w
                        end
                    end
                end
            end
        end
    end
    return bisectors, widths, depths
end


""" get_line_shapes(line_list, chunk_list, depths, v_width)
For a single spectrum in chunk_list, returns data frames for
line bisector and width as a function of fractional line depth given in depths
note: frac_depths must be sorted from largest to smallest
"""
function get_line_shapes(line_list::A1, chunk_list::AbstractChunkList;
         frac_depths::A2 = 0.99:-0.01:0.05, v_width::Real=15e3,
          RV_estimate::Real=-5000.0) where { T1<:Real, A1<:AbstractArray{T1,1},
          T2<:Real, A2<:AbstractArray{T2,1}}
    ll = copy(line_list)
    n_lines = length(ll)
    n_chunks = length(clt.chunk_list[1].data) # this assumes the number of chunks in every spectrum is the same
    n_depths = length(frac_depths)
    #allocate memory for line shapes
    bisectors = zeros(n_lines,n_depths)
    widths = zeros(n_lines,n_depths)
    depths = zeros(n_lines)
    #shift each entry in line list to RV estimate of star
    ll .*= calc_doppler_factor(RV_estimate)
    #find center of each chunk so that we can assign each line in line list to a chunk
    chunk_centers = zeros(n_chunks)
    for i in 1:n_chunks
        chunk_centers[i] = (first(chunk_list.data[i].λ) + last(chunk_list.data[i].λ)) / 2.0
    end
    #pick which chunk to measure each line shape in
    line_chunk_index = zeros(n_lines)
    for i in 1:n_lines
        line_chunk_index[i] = findmin(abs.((chunk_centers .- ll[i]) ./ chunk_centers))[2]
    end
    #measure the bisectors and widths
    for i in 1:n_lines
        #find pixels within v_width of ll[i]
        line_high = ll[i] * calc_doppler_factor(v_width)
        line_low = ll[i] * calc_doppler_factor(-v_width)
        index_low = searchsortedfirst(chunk_list.data[line_chunk_index[i]].λ,line_low) - 1
        index_high = searchsortedfirst(chunk_list.data[line_chunk_index[i]].λ,line_high)
        #get λs and fluxes for these pixels
        local_λ = view(chunk_list.data[line_chunk_index[i]].λ,index_low:index_high)
        local_flux = view(chunk_list.data[line_chunk_index[i]].flux,index_low:index_high)
        depths[i] = 1.0 - minimum(local_flux[ Int(floor(length(local_flux) / 4)) : Int(ceil(length(local_flux) * 3 / 4)) ]) / maximum(local_flux)
        for j in 1:n_depths
            b = RvSpectMLBase.calc_line_bisector_at_frac_depth(local_λ, local_flux, frac_depth=frac_depths[j])
            w = RvSpectMLBase.calc_line_width(local_λ, local_flux, frac_depth=frac_depths[j])
            # in the case of NaNs, use the previous depth's value (this is why depths must be sorted decreasing)
            if isnan(b)
                bisectors[i,j] = bisectors[i,j-1]
            else
                bisectors[i,j] = b
            end
            if isnan(w)
                widths[i,j] = widths[i,j-1]
            else
                widths[i,j] = w
            end
        end
    end
    return bisectors, widths, depths
end


