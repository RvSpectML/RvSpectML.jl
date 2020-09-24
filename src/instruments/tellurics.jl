"""
   find_worst_telluric_in_each_chunk( chunk_list_timeseries, array_of_spectra )
Returns a vector of the worst telluric (i.e., smallest value) within each chunk at any observations.

Warning: Assumes a tellurics value in metadata for each spectra, such as is provided by EXPRES.
"""
function find_worst_telluric_in_each_chunk end




"""
   make_clean_line_list_from_tellurics_expres(line_list, expres_data; Δv_to_avoid_tellurics )
Returns a new line list that excludes lines with telluric contamination.
Inputs:
- line_list:  Dataframe containing field lambda
- expres_data:  Array of spectra
- Δv_to_avoid_tellurics:  in m/s
Outputs:
- line_list_without_tellurics:   DataFrame with fields, lambda, weight, lambda_lo, and lambda_hi.
Warning: Currently, assumes a tellurics value in metadata for each spectra, such as is provided by EXPRES.
"""
function make_clean_line_list_from_tellurics_expres(line_list::DataFrame, expres_data::DT; Δv_to_avoid_tellurics::Real = default_Δv_to_avoid_tellurics,
            v_center_to_avoid_tellurics::Real = 0.0
               ) where { T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:EXPRES.AnyEXPRES, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1} }
   @assert 1000 <= Δv_to_avoid_tellurics <= 50000
   line_list_to_search_for_tellurics = copy(line_list)
   # TODO: Check sign convention for v_center_to_avoid_tellurics
   line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda./calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
   line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda.*calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
   chunk_list_timeseries = RvSpectML.make_chunk_list_timeseries(expres_data,line_list_to_search_for_tellurics)
   line_list_to_search_for_tellurics.min_telluric_model_all_obs = EXPRES.find_worst_telluric_in_each_chunk( chunk_list_timeseries, expres_data)
   line_list_no_tellurics_df = line_list_to_search_for_tellurics |> @filter(_.min_telluric_model_all_obs == 1.0) |> DataFrame
end
