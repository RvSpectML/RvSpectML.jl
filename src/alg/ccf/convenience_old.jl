"""
    Code for convenience functions for calculating CCFs on RvSpectML's types
Author: Eric Ford
Created: August 2020
"""

using ThreadedIterables

"""  calc_ccf_chunk_old( chunk, ccf_plan )
Convenience function to compute CCF for one chunk of spectrum.
Uses old version of project_mask_old! that is hardwired for a tophat mask.
# Inputs:
- chunk
# Optional Arguments:
- ccf_plan (BasicCCFPlan())
# Return:
CCF for one chunk of spectrum, evaluated using mask_shape and line list from ccf plan
"""
function calc_ccf_chunk_old(chunk::AbstractChuckOfSpectrum, plan::PlanT = BasicCCFPlan()
             ; assume_sorted::Bool = false ) where { PlanT<:AbstractCCFPlan, T1<:Real, T2<:Real }
  @assert assume_sorted || issorted( plan.line_list.λ )
  v_grid = calc_ccf_v_grid(plan)
  ccf_out = zeros(size(v_grid))
  ccf_1D_old!(ccf_out, chunk.λ, chunk.flux, plan; assume_sorted=true)
  return ccf_out
end

"""  calc_ccf_chunklist_old ( chunklist, ccf_plans )
Convenience function to compute CCF based on a spectrum's chunklist.
Uses old version of project_mask_old! that is hardwired for a tophat mask.
# Inputs:
- chunklist
- vector of ccf plans (one for each chunk)
# Optional Arguments:
# Return:
CCF summed over all chunks in a spectrum's chunklist, evaluated using the
line list and mask_shape from the ccf plan for each chunk.
"""
function calc_ccf_chunklist_old(chunk_list::AbstractChunkList,
                                plan_for_chunk::AbstractVector{PlanT};
                                assume_sorted::Bool = false ) where {
                                            PlanT<:AbstractCCFPlan }
  @assert length(chunk_list) == length(plan_for_chunk)
  mapreduce(chid->calc_ccf_chunk_old(chunk_list.data[chid], plan_for_chunk[chid], assume_sorted=assume_sorted), +, 1:length(chunk_list.data) )
end

"""  calc_ccf_chunklist_timeseries_old ( chunklist_timeseries, line_list )
Convenience function to compute CCF for a timeseries of spectra, each with a chunklist.
Uses old version of project_mask_old! that is hardwired for a tophat mask.
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
function calc_ccf_chunklist_timeseries_old(clt::AbstractChunkListTimeseries,
                                plan::PlanT = BasicCCFPlan(); verbose::Bool = false ) where {
                                    PlanT<:AbstractCCFPlan }

  @assert issorted( plan.line_list.λ )
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
  @threaded mapreduce(obsid->calc_ccf_chunklist_old(clt.chunk_list[obsid], plan_for_chunk, assume_sorted=true),hcat, 1:length(clt) )
end
