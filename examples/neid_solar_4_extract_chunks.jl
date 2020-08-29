# Run code for previous steps
include("neid_solar_1_read.jl")
# Don't import Plots if make_plots set to false
make_plots = isdefined(Main,:make_plots) ? make_plots : true
 if make_plots
   using Plots
 end

order_list_timeseries = RvSpectML.make_order_list_timeseries(solar_data)

order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries)

RvSpectML.normalize_spectra!(order_list_timeseries,solar_data);

if make_plots
   order_idx = 20:22
   plt = RvSpectML.plot_spectrum_chunks(order_list_timeseries, order_idx)
end

using DataFrames, CSV, Query
 using Statistics

vald_filename = joinpath(ancilary_data_path,"VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas")
vald_df = RvSpectML.read_mask_vald(vald_filename)

lambda_range_with_data = (min = maximum(d->minimum(d.λ),solar_data), max = minimum(d->maximum(d.λ),solar_data) )
line_list_df = vald_df |>
   @filter(lambda_range_with_data.min <= _.lambda_lo ) |>
   @filter( _.lambda_hi < lambda_range_with_data.max) |>
   @filter( _.lambda_hi < 6000.0 ) |>
#   @filter( _.lambda_lo >6157 || _.lambda_hi < 6155  ) |>   # Avoid "line" w/ large variability
   DataFrame

find_overlapping_chunks(line_list_df)

chunk_list_df = RvSpectML.merge_lines(line_list_df)
@assert find_overlapping_chunks(chunk_list_df) == nothing

size(solar_data)
chunk_list_df = chunk_list_df[1:end,:]

# TODO:  FIX:  Problem with ESPRESSO mask having lots of lines getting meged into humongous chunks that don't fit into the order!
chunk_list_timeseries = RvSpectML.make_chunk_list_timeseries(solar_data,chunk_list_df)

# Check that no NaN's included
(chunk_list_timeseries, chunk_list_df) = RvSpectML.filter_bad_chunks(chunk_list_timeseries,chunk_list_df)
println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )

using Plots
if make_plots
   #order_idx = 20:22
   # plt = RvSpectML.plot_spectrum_chunks(order_list_timeseries, order_idx)
   chunk_idx = 23
   for c in chunk_idx
       t = 1
       #if(sum(chunk_list_df.line_depths[c])<0.25) continue end
       λ_mid = sqrt(chunk_list_df.lambda_hi[c]*chunk_list_df.lambda_lo[c])
       println("c= ",c , " λs= ",chunk_list_df.lambda[c]," depths= ",chunk_list_df.depth[c])
       flush(stdout)
   end
   plt = RvSpectML.plot_spectrum_chunks(chunk_list_timeseries, chunk_idx, plt=plot(), color=:black)
   display(plt)

end
