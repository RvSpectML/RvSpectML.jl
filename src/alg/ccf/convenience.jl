"""
    Code for convenience functions for calculating CCFs on RvSpectML's types
Author: Eric Ford
Created: August 2019
"""

"""  calc_ccf_chunk( chunk, line_list::ALL )
Convenience function to compute CCF for one chunk of spectrum.
# Inputs:
- chunk
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
CCF for one chunk of spectrum, evaluated using mask_shape and plan
"""
function calc_ccf_chunk(chunk::AbstractChuckOfSpectrum,
                            #line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan::PlanT = BasicCCFPlan() ) where {
                                    PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                                    # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  v_grid = calc_ccf_v_grid(plan)
  ccf_out = zeros(size(v_grid))
  #ccf_1D!(ccf_out, chunk.λ, chunk.flux, line_list, mask_shape=mask_shape, plan=plan)
  ccf_1D!(ccf_out, chunk.λ, chunk.flux, plan)
  # TODO Change once ready use generic CCF mask shape (currently don't understand why normalization of output CCF differs)
  #ccf_1D_expr!(ccf_out, chunk.λ, chunk.flux, plan)
  return ccf_out
end

"""  calc_ccf_chunk_expr( chunk, line_list::ALL )
Convenience function to compute CCF for one chunk of spectrum.
Uses experimental version of project_mask_expr!.
Need to understand why difference before merging this in.
# Inputs:
- chunk
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
CCF for one chunk of spectrum, evaluated using mask_shape and plan
"""
function calc_ccf_chunk_expr(chunk::AbstractChuckOfSpectrum,
                            #line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan::PlanT = BasicCCFPlan() ) where {
                                    PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                                    # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  v_grid = calc_ccf_v_grid(plan)
  ccf_out = zeros(size(v_grid))
  #ccf_1D!(ccf_out, chunk.λ, chunk.flux, line_list, mask_shape=mask_shape, plan=plan)
  #ccf_1D!(ccf_out, chunk.λ, chunk.flux, plan)
  # TODO Change once ready use generic CCF mask shape (currently don't understand why normalization of output CCF differs)
  ccf_1D_expr!(ccf_out, chunk.λ, chunk.flux, plan)
  return ccf_out
end

"""  calc_ccf_chunklist ( chunklist, line_list )
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using mask_shape and plan.
"""
function calc_ccf_chunklist(chunk_list::AbstractChunkList,
                            #line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
                                plan_for_chunk::AbstractVector{PlanT}  ) where {
                                    #ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                                            PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  @assert length(chunk_list) == length(plan_for_chunk)
  #mapreduce(chunk->calc_ccf_chunk(chunk, line_list,mask_shape=mask_shape, plan=plan), +, chunk_list.data)
  mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid]), +, 1:length(chunk_list.data) )
end

"""  calc_ccf_chunklist_timeseries( chunklist_timeseries, line_list )
Convenience function to compute CCF for a timeseries of spectra, each with a chunklist.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
- verbose (false)
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the ccf_plan.
Note that the ccf_plan provided is used as a template for creating a custom ccf_plan for each chunk that 
    only includes lines that reliably appear in that order for all spectra in the chunklist_timeseries.
"""
function calc_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                                plan::PlanT = BasicCCFPlan(); verbose::Bool = false ) where {
                                    PlanT<:AbstractCCFPlan }

  num_lines = length(plan.line_list.λ)
  plan_for_chunk = Vector{BasicCCFPlan}(undef,num_chunks(clt))
  for chid in 1:num_chunks(clt)
      λmin = maximum(map(obsid->first(clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
      λmax = minimum(map(obsid->last( clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
      upper_edge_of_mask_for_line_at_λmin = λ_max(plan.mask_shape,λmin)
      lower_edge_of_mask_for_line_at_λmax = λ_min(plan.mask_shape,λmax)
      start_line_idx = searchsortedfirst(plan.line_list.λ,upper_edge_of_mask_for_line_at_λmin)
      if verbose
          flush(stdout)
          println("extrema(plan.line_list.λ) = ",extrema(plan.line_list.λ) )
          println("upper_edge_of_mask_for_line_at_λmin = ", upper_edge_of_mask_for_line_at_λmin, " LOWer_edge_of_mask_for_line_at_λmax = ", lower_edge_of_mask_for_line_at_λmax)
      end
      stop_line_idx  = num_lines+1 - searchsortedfirst(view(plan.line_list.λ,num_lines:-1:1),lower_edge_of_mask_for_line_at_λmax,rev=true)
      if (1 <= start_line_idx <= length(plan.line_list.λ)) && (1 <= stop_line_idx <= length(plan.line_list.λ))
          if verbose
               println("start_line_idx = ", start_line_idx, " λ= ", plan.line_list.λ[start_line_idx])
               println("stop_line_idx = ", stop_line_idx, " λ= ", plan.line_list.λ[stop_line_idx])
               println(" Using ", length(start_line_idx:stop_line_idx), " lines for chunk ", chid)
           end
          line_list_for_chunk = BasicLineList(view(plan.line_list.λ,start_line_idx:stop_line_idx), view(plan.line_list.weight,start_line_idx:stop_line_idx) )
      else  # No lines in this chunk!
          line_list_for_chunk = EmptyBasicLineList()
      end
      plan_for_chunk[chid] = BasicCCFPlan( line_list=line_list_for_chunk, midpoint=plan.v_center, step=plan.v_step, max=plan.v_max, mask_shape=plan.mask_shape )
  end
  @threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk),hcat, 1:length(clt) )
end

"""  calc_order_ccfs_chunklist ( chunklist_timeseries, line_list )
Convenience function to compute separate CCFs for each chunk in a spectrum.
CCF evaluated using mask_shape and plan.
# Inputs:
- chunklist_timeseries
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccfs_chunklist(chunk_list::AbstractChunkList, # line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
    plan::PlanT = BasicCCFPlan() ) where {
                        PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
        # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
    #mapreduce(chunk->calc_ccf_chunk(chunk, line_list, mask_shape=mask_shape, plan=plan), hcat, chunk_list.data)
    mapreduce(chunk->calc_ccf_chunk(chunk, plan), hcat, chunk_list.data)
end

"""  calc_order_ccf_chunklist_timeseries( chunklist_timeseries, line_list )
Convenience function to compute separate CCFs for each chunk of each spectrum in a timeseries.
CCF is evaluated using mask_shape and plan.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
- line_list
# Optional Arguments:
- mask_shape (TopHatCCFMask())
- ccf_plan (BasicCCFPlan())
# Return:
A 3-d array containing the CCF at each (velocity, chunk, spectrum)
"""
function calc_order_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                plan::PlanT = BasicCCFPlan() ) where {
        #                    line_list::ALL; mask_shape::ACMS = TopHatCCFMask(),
    #                            plan::AbstractCCFPlan = BasicCCFPlan() ) where {
#                                    ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                    PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
  nvs = length(calc_ccf_v_grid(plan))
  norders = length(clt.chunk_list[1].data)
  nobs =  length(clt.chunk_list)
  order_ccfs = zeros(nvs, norders, nobs)
  Threads.@threads for i in 1:nobs
      #order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], line_list,  mask_shape=mask_shape, plan=plan)
      order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], plan)
  end
  return order_ccfs
end

#=
"""
   find_bin_edges( pixel_centers )

Internal function used by project_mask!.
"""
function find_bin_edges_opt(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = Array{eltype(fws),1}(undef,length(fws)+2)
    fwse[2:end-2] = (fws[2:end] + fws[1:end-1]) ./ 2.0
    fwse[1] = 2.0 * fws[1] - fwse[2]
    fwse[end-1] = 2.0 * fws[end] - fwse[end-2]
    fwse[end] = zero(eltype(A))
    return fwse
end


function find_bin_edges_compare(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = (fws[2:end] + fws[1:end-1]) ./ 2.0
    le = 2.0 * fws[1] - fwse[1]
    re = 2.0 * fws[end] - fwse[end]
    fwsle = vcat(le, fwse)
    fwsre = vcat(fwse, re)
    return fwsle, fwsre
end



function project_mask_compare!(projection::A2, λs::A1, plan::PlanT = BasicCCFPlan() ) where {
    #mask_in::ALL; shift_factor::Real=one(T1),
            #mask_shape::ACMS = TopHatCCFMask() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2},
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                #plan::PlanT = BasicCCFPlan() ) where {
                #ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }

#function project_mask_compare(fws::A1, mask::A2; shift_factor::Real=one(T1)) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    # find bin edges
    fwsle, fwsre = find_bin_edges_compare(λs)

    # allocate memory for line_list projection
    nm = length(mask_in)
    nx = size(fwsle, 1)
    ny = size(fwsle, 2)
    projection = zeros(nx, ny)
    #flush(stdout);  println("size(projection) = ", size(projection), " size(fws) = ", size(fws), " size(mask) = ", size(mask), " shift_factr = ", shift_factor)

    # shift the mask
    mask_shifted = zeros(length(mask_in),3)
    mask_shifted[:,1] = plan.mask_shape.λ_lo  .* shift_factor  # view(mask,:,1) .* shift_factor
    mask_shifted[:,2] = plan.mask_shape.λ_hi  .* shift_factor  # view(mask,:,2) .* shift_factor
    mask_shifted[:,3] = plan.mask_shape.weight                 # view(mask,:,3)
    local mask = mask_shifted

    # set indexing variables
    p = 1
    m = 1
    on_mask = false

    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if fwsre[p] > mask[m,1]
                if fwsre[p] > mask[m,2]
                    projection[p] += mask[m,3] * (mask[m,2] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    m += 1
                else
                    projection[p] += mask[m,3] * (fwsre[p] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    on_mask = true
                    p+=1
                end
            else
                p+=1
            end
        else
            if fwsre[p] > mask[m,2]
                projection[p] += mask[m,3] * (mask[m,2] - fwsle[p]) / (fwsre[p] - fwsle[p])
                on_mask = false
                m += 1
            else
                projection[p] += mask[m,3]
                p += 1
            end
        end
        if p>length(projection) break end
    end
    return projection
end
=#
