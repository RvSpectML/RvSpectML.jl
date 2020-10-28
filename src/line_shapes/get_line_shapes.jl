
"""
Code to measure line bisectors and widths as a function of depth.
Author: Alex Wise
Created: September 2020
"""

""" get_line_shapes(line_list, clt, frac_depths, v_width, RV_estimate, verbose)
For each spectrum in clt, returns arrays for
line bisector and width as a function of fractional line depth given in frac_depths.
Note that frac_depths must be sorted from largest to smallest.
Also returns arrays for line depth, presently a simple minimum over the line's pixels, and valid lines,
a boolean array for which lines' shape measurements were calculated successfully.
v_width is the width of the window (in m/s) around each line's center to use.
RV_estimate (in m/s) is a rough estimate of the star's RV.
Note: RV_estimate shifts the line centers in line_list to the RV of star, so use RV_estimate=0 if line_list is already aligned to star.
"""
function get_line_shapes(line_list::A1, clt::AbstractChunkListTimeseries;
                    frac_depths::A2 = 0.99:-0.01:0.05, v_width::Real=15e3,
                    RV_estimate::Real=-5000.0, verbose::Bool = true) where { T1<:Real, A1<:AbstractArray{T1,1},
                    T2<:Real, A2<:AbstractArray{T2,1}}

    @assert issorted(frac_depths, rev=true) #make sure frac_depths are sorted largest to smallest
    n_lines = length(line_list)
    n_spectra = length(clt.chunk_list)
    n_depths = length(frac_depths)
    #allocate memory for line shapes
    bisectors = zeros(n_depths,n_spectra,n_lines)
    widths = zeros(n_depths,n_spectra,n_lines)
    depths = zeros(n_spectra,n_lines)
    valid_lines = ones(Bool,n_lines)
    #find center of each chunk so that we can assign each line in line list to a chunk
    chunk_centers = zeros(length(clt.chunk_list[1])) #note: this assumes all spectra have the same number of chunks
    calc_chunk_centers!(chunk_centers, clt.chunk_list[1]) # TODO: change chunk centers to average over all chunks in timeseries
    #measure the bisectors and widths
    for k in 1:n_spectra
        line_shapes = get_line_shapes(line_list, clt.chunk_list[k], chunk_centers=chunk_centers,
        frac_depths=frac_depths, v_width=v_width, RV_estimate=RV_estimate, verbose=verbose)
        bisectors[:,k,:] = line_shapes[:bisectors]
        widths[:,k,:] = line_shapes[:widths]
        depths[k,:] = line_shapes[:depths]
        valid_lines .= valid_lines .& line_shapes[:valid_lines]
    end
    #return bisectors, widths, depths, and valid indices
    (bisectors = bisectors, widths=widths, depths=depths, valid_lines=valid_lines)
end


""" get_line_shapes(line_list, chunk_list, chunk_centers, frac_depths, v_width, RV_estimate, verbose)
For a single spectrum in chunk_list, returns arrays for
line bisector and width as a function of fractional line depth given in frac_depths.
Note that frac_depths must be sorted from largest to smallest.
Also returns arrays for line depth, presently a simple minimum over the line's pixels, and valid lines,
a boolean array for which lines' shape measurements were calculated successfully.
v_width is the width of the window (in m/s) around each line's center to use.
RV_estimate (in m/s) is a rough estimate of the star's RV.
Note: RV_estimate shifts the line centers in line_list to the RV of star, so use RV_estimate=0 if line_list is already aligned to star.
"""
function get_line_shapes(line_list::A1, chunk_list::AbstractChunkList;
        chunk_centers::A2 = calc_chunk_centers!(zeros(length(chunk_list)),chunk_list),
        frac_depths::A3 = 0.99:-0.01:0.05, v_width::Real=15e3,
          RV_estimate::Real=-5000.0, verbose::Bool = true) where { T1<:Real, A1<:AbstractArray{T1,1},
          T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1}}
    @assert issorted(frac_depths, rev=true) #make sure frac_depths are sorted largest to smallest
    ll = deepcopy(line_list)
    #shift each entry in line list to RV estimate of star
    ll .*= calc_doppler_factor(RV_estimate)
    #get array dimensions
    n_lines = length(ll)
    n_depths = length(frac_depths)
    #allocate memory for line shapes - could be optimized by referencing memory in get_line_shapes(liine_list,clt)
    bisectors = zeros(n_depths,n_lines)
    widths = zeros(n_depths,n_lines)
    depths = zeros(n_lines)
    valid_lines = ones(Bool,n_lines)
    #pick which chunk to measure each line shape in
    line_chunk_index = zeros(Int64,n_lines)
    for i in 1:n_lines
        line_chunk_index[i] = findmin(abs.((chunk_centers .- ll[i]) ./ chunk_centers))[2] # TODO: change line_chunk_index to pick higher SNR of valid chunks
    end
    #measure the bisectors and widths
    #TOTO optional: might want to make this into a mutating function that operates on a single line, and is called here.
    for i in 1:n_lines
        #find pixels within v_width of ll[i]
        line_high = ll[i] * calc_doppler_factor(v_width)
        line_low = ll[i] * calc_doppler_factor(-v_width)
        index_low = searchsortedfirst(chunk_list.data[line_chunk_index[i]].λ,line_low) - 1
        index_high = searchsortedfirst(chunk_list.data[line_chunk_index[i]].λ,line_high)
        if (index_low < 1) | (index_high > length(chunk_list.data[line_chunk_index[i]].λ))
            if verbose
                println("Warning: line number ",i," is not contained in chunk_list.")
            end
            valid_lines[i] = false
        else
            #get λs and fluxes for these pixels
            local_λ = view(chunk_list.data[line_chunk_index[i]].λ,index_low:index_high)
            local_flux = view(chunk_list.data[line_chunk_index[i]].flux,index_low:index_high)
            #calculate the depth by taking the minimum of the middle 50% of the line's pixels #TODO: add make this fraction an argument as in calc_line_width()
            depths[i] = 1.0 - minimum(local_flux[ Int(floor(length(local_flux) / 4)) : Int(ceil(length(local_flux) * 3 / 4)) ]) / maximum(local_flux)
            for j in 1:n_depths
                #calculate the line widths and bisectors # TODO opt: could merge calc_line_bis with calc_line width for performance
                b = RvSpectMLBase.calc_line_bisector_at_frac_depth(local_λ, local_flux, frac_depth=frac_depths[j])
                w = RvSpectMLBase.calc_line_width(local_λ, local_flux, frac_depth=frac_depths[j])
                @assert j != 1 || (!isnan(b) && !isnan(w)) #line minumum is wrong or line doesn't exist
                # in the case of NaNs, use the previous depth's value (this is why depths must be sorted decreasing)
                if isnan(b)
                    bisectors[j,i] = bisectors[j-1,i]
                else
                    bisectors[j,i] = b
                end
                if isnan(w)
                    widths[j,i] = widths[j-1,i] #TODO: figure out what to do here to make more sense
                else
                    widths[j,i] = w
                end #end "if isnan(w)"
            end #end "for j in 1:n_depths"
        end #end if line is not in chunk_list
    end #end "for i in 1:n_lines"
    #return bisectors, widths, depths, and valid indices
    (bisectors = bisectors, widths=widths, depths=depths, valid_lines=valid_lines)
end


function calc_chunk_centers!(chunk_centers::A1, chunk_list::AbstractChunkList) where { T1<:Real, A1<:AbstractArray{T1,1}}
    for i in 1:length(chunk_centers)
        chunk_centers[i] = (first(chunk_list.data[i].λ) + last(chunk_list.data[i].λ)) / 2.0
    end
end

