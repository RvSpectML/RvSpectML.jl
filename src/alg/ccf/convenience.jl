"""
    Code for convenience functions for calculating CCFs on RvSpectML's types
Author: Eric Ford
Created: August 2020
"""

using Statistics
using ThreadedIterables

"""  calc_ccf_chunk( chunk, ccf_plan )
Convenience function to compute CCF for one chunk of spectrum.
# Inputs:
- chunk
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
# Return:
CCF for one chunk of spectrum, evaluated using mask_shape and line list from ccf plan
"""
function calc_ccf_chunk(chunk::AbstractChunkOfSpectrum, plan::PlanT = BasicCCFPlan()
                 ; var::AbstractVector{T} = chunk.var,
                 assume_sorted::Bool = false, use_pixel_vars::Bool = false ) where { T<:Real, PlanT<:AbstractCCFPlan }
  @assert assume_sorted || issorted( plan.line_list.λ )
  len_v_grid = calc_length_ccf_v_grid(plan)
  ccf_out = zeros(len_v_grid)
  if use_pixel_vars
      ccf_1D!(ccf_out, chunk.λ, chunk.flux, var, plan; assume_sorted=true)
  else
      ccf_1D!(ccf_out, chunk.λ, chunk.flux, plan; assume_sorted=true)
  end
  return ccf_out
end

"""  calc_ccf_chunklist ( chunklist, ccf_plans )
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_chunklist(chunk_list::AbstractChunkList,
                                plan_for_chunk::AbstractVector{PlanT};
                                assume_sorted::Bool = false, use_pixel_vars::Bool = false   ) where {
                                            PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid],
                    assume_sorted=assume_sorted, use_pixel_vars=use_pixel_vars), +, 1:length(chunk_list.data) )
end


"""  calc_ccf_chunklist ( chunklist, var_list, ccf_plans )
Convenience function to compute CCF based on a spectrum's chunklist.
# Inputs:
- chunklist
- var_list: vector of variance vectors for each chunk
- ccf_palns: vector of ccf plans (one for each chunk)
# Optional Arguments:
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_chunklist(chunk_list::AbstractChunkList, var_list::AbstractVector{A1},
                                plan_for_chunk::AbstractVector{PlanT};
                                assume_sorted::Bool = false, use_pixel_vars::Bool = false   ) where {
                                            T1<:Real, A1<:AbstractVector{T1}, PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  @assert length(chunk_list) == length(var_list)
  @assert all(map(chid->length(chunk_list.data[chid].var) == length(var_list[chid]),1:length(chunk_list.data)))
  #var_weight = mapreduce(chid->median(chunk_list.data[chid].var), +, 1:length(chunk_list.data) ) / mapreduce(chid->median(var_list[chid]), +, 1:length(chunk_list.data) )
  var_weight = 1.0 #mapreduce(chid->mean(chunk_list.data[chid].var), +, 1:length(chunk_list.data) ) / mapreduce(chid->mean(var_list[chid]), +, 1:length(chunk_list.data) )
  #println("# var_weight*N_obs = ", var_weight .* length(chunk_list))
  mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid], var=var_list[chid] .* var_weight,
                    assume_sorted=assume_sorted, use_pixel_vars=use_pixel_vars), +, 1:length(chunk_list.data) )
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
                                plan::PlanT = BasicCCFPlan(); verbose::Bool = false,
                                use_pixel_vars::Bool = false  ) where { PlanT<:AbstractCCFPlan }

  @assert issorted( plan.line_list.λ )
  num_lines = length(plan.line_list.λ)
  plan_for_chunk = Vector{BasicCCFPlan}(undef,num_chunks(clt))
  if use_pixel_vars
      var_list = Vector{Vector{Float64}}(undef,num_chunks(clt))
  end
  for chid in 1:num_chunks(clt)
      # find the maximum lower wavelength for the chunk, and the minumum upper wavelength, over all observations
      λmin = maximum(map(obsid->first(clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
      λmax = minimum(map(obsid->last( clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
      # extend λmin/λmax by the velocity range over which we don't want the mask to change
      λmin  = λmin / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(-plan.v_range_no_mask_change))
      λmax  = λmax / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(plan.v_range_no_mask_change))
      # extend λmin/λmax by the mask width
      upper_edge_of_mask_for_line_at_λmin = λ_max(plan.mask_shape,λmin)
      lower_edge_of_mask_for_line_at_λmax = λ_min(plan.mask_shape,λmax)
      # find the first and last mask entries to use in each chunk
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
      #create a plan for this chunk that only includes the mask entries we want for this chunk
      plan_for_chunk[chid] = BasicCCFPlan( line_list=line_list_for_chunk, midpoint=plan.v_center, step=plan.v_step, max=plan.v_max, mask_shape=plan.mask_shape )
      if use_pixel_vars
          # Warning: Orders have consistent length over time, but chunks may not.  So currently only works with full orders
          var_list[chid] = mapreduce(obsid->clt.chunk_list[obsid].data[chid].var,+, 1:length(clt) )
          var_list[chid] ./= length(clt)
          #var_list[chid] = vec(median(mapreduce(obsid->clt.chunk_list[obsid].data[chid].var,hcat, 1:length(clt) ),dims=2))
            #return var_list
      end
  end
  if use_pixel_vars
      return @threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], var_list, plan_for_chunk,assume_sorted=true, use_pixel_vars=true),hcat, 1:length(clt) )
  else
      return @threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk,assume_sorted=true, use_pixel_vars=false),hcat, 1:length(clt) )
  end
end


"""  calc_order_ccfs_chunklist ( chunklist_timeseries, line_list )
Convenience function to compute separate CCFs for each chunk (typically an order) in a spectrum.
CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
# Inputs:
- chunklist_timeseries
- vector of ccf plans (one for each chunk)
# Return:
A 2-d array containing the CCF at each (velocity, chunk)
"""
function calc_order_ccfs_chunklist(chunk_list::AbstractChunkList,
    plan_for_chunk::AbstractVector{PlanT} = BasicCCFPlan(); assume_sorted::Bool = false, use_pixel_vars::Bool = false  ) where {
                        PlanT<:AbstractCCFPlan }
    @assert length(chunk_list) == length(plan_for_chunk)
    @assert assume_sorted || issorted( first(plan_for_chunk).line_list.λ )
    mapreduce(chid->calc_ccf_chunk(chunk_list.data[chid], plan_for_chunk[chid], assume_sorted=assume_sorted, use_pixel_vars=use_pixel_vars), hcat, 1:length(chunk_list.data) )
end

"""  calc_order_ccf_chunklist_timeseries( chunklist_timeseries, ccf_plan )
Convenience function to compute separate CCFs for each chunk (typically an order) of each spectrum in a timeseries.
    CCF is evaluated using line list and mask_shape provided by the ccf plan for each chunk.
Uses multiple threads if avaliable.
# Inputs:
- chunklist_timeseries
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
# Return:
A 3-d array containing the CCF at each (velocity, chunk, spectrum)
Note that the ccf_plan provided is used as a template for creating a custom ccf_plan for each chunk that
    only includes lines that reliably appear in that order for all spectra in the chunklist_timeseries.
"""
function calc_order_ccf_chunklist_timeseries(clt::AbstractChunkListTimeseries,
                plan::PlanT = BasicCCFPlan(); verbose::Bool = false, use_pixel_vars::Bool = false  ) where {
                    PlanT<:AbstractCCFPlan }
  @assert issorted( plan.line_list.λ )
  #nvs = length(calc_ccf_v_grid(plan))
  nvs = calc_length_ccf_v_grid(plan)
  norders = length(clt.chunk_list[1].data)
  nobs =  length(clt.chunk_list)
  order_ccfs = zeros(nvs, norders, nobs)
  num_lines = length(plan.line_list.λ)
  plan_for_chunk = Vector{BasicCCFPlan}(undef,num_chunks(clt))
  for chid in 1:num_chunks(clt)
      # find the maximum lower wavelength for the chunk, and the minumum upper wavelength, over all observations
      λmin = maximum(map(obsid->first(clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
      λmax = minimum(map(obsid->last( clt.chunk_list[obsid].data[chid].λ), 1:length(clt) ))
      # extend λmin/λmax by the velocity range over which we don't want the mask to change
      λmin  = λmin / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(-plan.v_range_no_mask_change))
      λmax  = λmax / (calc_doppler_factor(plan.v_center)*calc_doppler_factor(plan.v_range_no_mask_change))
      # extend λmin/λmax by the mask width
      upper_edge_of_mask_for_line_at_λmin = λ_max(plan.mask_shape,λmin)
      lower_edge_of_mask_for_line_at_λmax = λ_min(plan.mask_shape,λmax)
      # find the first and last mask entries to use in each chunk
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
      #create a plan for this chunk that only includes the mask entries we want for this chunk
      plan_for_chunk[chid] = BasicCCFPlan( line_list=line_list_for_chunk, midpoint=plan.v_center, step=plan.v_step, max=plan.v_max, mask_shape=plan.mask_shape )
  end
  #@threaded mapreduce(obsid->calc_ccf_chunklist(clt.chunk_list[obsid], plan_for_chunk),hcat, 1:length(clt) )
  #Threads.@threads for i in 1:nobs
    #  order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], plan)
  #end
  Threads.@threads for i in 1:nobs
    order_ccfs[:,:,i] .= calc_order_ccfs_chunklist(clt.chunk_list[i], plan_for_chunk, assume_sorted=true, use_pixel_vars=use_pixel_vars)
  end

  return order_ccfs
end
