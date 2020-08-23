using RvSpectML

# 2458745.1296134139 to 2458745.4284185418. Out of those, I usually also require the JD to be less than 2458745.283
data_path = "20190918"   # TODO:  Update with path to data on your local machine
num_spectra_to_use = 20
bjd_first_good = 2458745.1296134139
bjd_last_good = 2458745.283

df_files = make_manifest(data_path)
#=
dir_filelist = readdir(data_path,join=true)
idx_spectra = map(fn->occursin(r"^neid\w+\.fits$", last(split(fn,'/')) ),dir_filelist)
spectra_filelist = dir_filelist[idx_spectra]

df_files = DataFrame(Filename = String[], target = String[], bjd = Float64[], ssbz=Float64[] )
map(fn->add_metadata_from_files_with_neid_solar_data!(df_files,fn),spectra_filelist)
=#

df_files_use = df_files |>
  @filter( _.target == "Solar" ) |>
  @filter(bjd_first_good <= _.bjd < bjd_last_good) |>
  @take(num_spectra_to_use) |>
  DataFrame

@time solar_data = map(read_neid_solar_data,eachrow(df_files_use))



if true
    drift_corrections = CSV.read("SolarRV20190918_JD_SciRV_CalRV.txt", header=["bjd", "sci_drift", "cal_drift"]);

    drift_interp = LinearInterpolation(drift_corrections[:bjd],drift_corrections[:cal_drift])
    df_files_use[!,:drift] = drift_interp.(df_files_use.bjd)
    if ! "doppler_factor" in names(df_files_use)
        df_files_use[!,:doppler_factor] = calc_doppler_factor.(df_files_use[:drift])
    else
        df_files_use[!,:doppler_factor] .*= calc_doppler_factor.(df_files_use[:drift])
    end

    @time map(x->apply_doppler_factor!(x[1],x[2]), zip(solar_data,df_files_use.doppler_factor));

    # TODO: WARNING: NEED TO CHECK SIGN OF DRIFT CORRECTION!!!
#=
function apply_drift_correction!(spectra::S,time::Real,drift_cor::DataFrame; rv_field::Symbol=:cal_drift) where {S<:AbstractSpectra}
    @assert any(isequal(:bjd),propertynames(drift_cor))
    @assert any(isequal(rv_field),propertynames(drift_cor))
    #@assert haskey(spectra.metadata,:doppler_factor)
    idx = findmin(abs.(time .- drift_cor.bjd))[2]
    doppler_factor = calc_doppler_factor(drift_cor[rv_field][idx])
    # println("# t= ",time, " drift_cor= ",drift_cor[rv_field][idx], " df= ",doppler_factor)
    spectra.λ .*= doppler_factor
    if haskey(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor
    end
    return spectra
end
@time map(x->apply_drift_correction!(x[1],x[2],drift_corrections), zip(solar_data,df_files_use.bjd));
drift_corrections[:cal_drift].-mean(drift_corrections[:cal_drift])
=#

end

if true
    ssb_corrections = CSV.read("SolarTelescope2019-09-18_inclGravRedshiftAirMassAltitude.csv", header=["bjd","rv_ssb"], select=[1,2], types=types=[Float64,Float64], datarow=2, silencewarnings=true);
@assert any(isequal(:bjd),propertynames(ssb_corrections))
@assert any(isequal(:rv_ssb),propertynames(ssb_corrections))
ssb_interp = LinearInterpolation(ssb_corrections.bjd, ssb_corrections.rv_ssb)

function apply_ssb_correction!(spectra::AS,time::Real,ssb_cor::AI) where {AS<:AbstractSpectra,AI<:AbstractInterpolation}
    doppler_factor = calc_doppler_factor(ssb_cor(time))
    #println("# t= ",time, " ssb_cor= ",ssb_cor(time), " df= ",doppler_factor)
    spectra.λ .*= doppler_factor
    if haskey(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = doppler_factor
    end
    return spectra
end
@time map(x->apply_ssb_correction!(x[1],x[2],ssb_interp), zip(solar_data,df_files_use.bjd));
ssb_interp.(df_files_use.bjd)
end
