"""
   find_worst_telluric_in_each_chunk( chunk_list_timeseries, array_of_spectra )
Returns a vector of the worst telluric (i.e., smallest value) within each chunk at any observations.

Warning: Assumes a tellurics value in metadata for each spectra, such as is provided by EXPRES.
"""
function find_worst_telluric_in_each_chunk( clt::AbstractChunkListTimeseries, data::AbstractArray{AS,1} )  where {AS<:AbstractSpectra}
   @assert typeof(first(data).inst) <: EXPRES.AnyEXPRES
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
