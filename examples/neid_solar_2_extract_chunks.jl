# Run code for previous steps
include("neid_solar_1_read.jl")
# Don't import Plots if make_plots set to false
make_plots = isdefined(Main,:make_plots) ? make_plots : true
 if make_plots
   using Plots
 end

order_list_timeseries = RvSpectML.make_order_list_timeseries(solar_data) # ,chunk_list_df)

order_list_timeseries = RvSpectML.filter_bad_chunks(order_list_timeseries)

RvSpectML.normalize_spectra!(order_list_timeseries,solar_data);

if make_plots
   order_idx = 20:22
   plt = RvSpectML.plot_spectrum_chunks(order_list_timeseries, order_idx)
end

using DataFrames, CSV, Query
 using Statistics

vald_filename = joinpath(ancilary_data_path,"VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas")
lambda_range_with_data = (min = maximum(d->minimum(d.λ),solar_data), max = minimum(d->maximum(d.λ),solar_data) )
vald_df = RvSpectML.read_mask_vald(vald_filename)
 line_list_df = vald_df |>
   @filter(lambda_range_with_data.min <= _.lambda_lo ) |>
   @filter( _.lambda_hi < lambda_range_with_data.max) |>
#   @filter( _.lambda_lo >6157 || _.lambda_hi < 6155  ) |>   # Avoid "line" w/ large variability
   DataFrame

# TODO: Once find good values, move this code to adjust width of chunks into read_mask_vald
chunk_size_factor = 3       # TODO: Figure out what value to use
 max_vald_line_offset = 0.0       # km/s
 line_width = RvSpectML.predict_line_width(5780,v_rot=1.8) # # km/s
 Δλoλ_fit_line = (max_vald_line_offset+chunk_size_factor*line_width)*1000/RvSpectML.speed_of_light_mps
 println("# Δλ/λ = ",Δλoλ_fit_line)
 line_list_df.lambda_hi .= line_list_df.lambda*(1 + Δλoλ_fit_line)
 line_list_df.lambda_lo .= line_list_df.lambda/(1 + Δλoλ_fit_line)

chunk_list_df = line_list_df

find_overlapping_chunks(line_list_df)

chunk_list_df = RvSpectML.merge_lines(line_list_df)
 @assert find_overlapping_chunks(chunk_list_df) == nothing

chunk_list_timeseries = RvSpectML.make_chunk_list_timeseries(solar_data,chunk_list_df)
println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )
# Check that no NaN's included
(chunk_list_timeseries, chunk_list_df) = RvSpectML.filter_bad_chunks(chunk_list_timeseries,chunk_list_df)
println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )

using Plots
if make_plots
   #order_idx = 20:22
   # plt = RvSpectML.plot_spectrum_chunks(order_list_timeseries, order_idx)
   chunk_idx = 13:22
   RvSpectML.plot_spectrum_chunks(chunk_list_timeseries, chunk_idx, plt=plt, color=:black)
end

#= Goofing around to see how well chunks line up
chunk_idx = 13:20
xmin = minimum(chunk_list_df.lambda_lo[chunk_idx])
xmax = maximum(chunk_list_df.lambda_hi[chunk_idx])
plt = plot(legend=:none)
#xlims!(xmin,xmax)
for c in chunk_idx
    t = 1
    #if(sum(chunk_list_df.line_depths[c])<0.25) continue end
    λ_mid = sqrt(chunk_list_df.lambda_hi[c]*chunk_list_df.lambda_lo[c])
    println("c= ",c , " λs= ",chunk_list_df.line_λs[c]," depths= ",chunk_list_df.line_depths[c])
    #println("  λlo= ",chunk_list_df.lambda_lo[c]," λhi= ",chunk_list_df.lambda_hi[c], " Δλ= ",chunk_list_df.lambda_hi[c]-chunk_list_df.lambda_lo[c])
    plot!(plt,chunk_list_timeseries.chunk_list[t].data[c].λ.-λ_mid,chunk_list_timeseries.chunk_list[t].data[c].flux)
end
#plot!(plt,solar_data[1].λ,solar_data[1].flux)
#xlims!(4560,4565)
display(plt)
=#
