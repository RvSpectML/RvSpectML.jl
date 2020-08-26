using Pkg
Pkg.activate(".")

using Revise
using RvSpectML
using DataFrames, Query

# 2458745.1296134139 to 2458745.4284185418. Out of those, I usually also require the JD to be less than 2458745.283
solar_data_path = "20190918"   # TODO:  Update with path to data on your local machine
ancilary_data_path = "."
 if isfile(joinpath(pkgdir(RvSpectML),"examples","data_paths.jl"))
    include(joinpath(pkgdir(RvSpectML),"examples","data_paths.jl"))
 end
 solar_data_path

num_spectra_to_use = 20
 bjd_first_good = 2458745.1296134139
 bjd_last_good = 2458745.283
 df_files = NEID.make_manifest(solar_data_path)

df_files_use = df_files |>
  @filter( _.target == "Solar" ) |>
  @filter(bjd_first_good <= _.bjd < bjd_last_good) |>
#  @take(num_spectra_to_use) |>
  DataFrame

solar_data = map(NEID.read_solar_data,eachrow(df_files_use))

NEID.read_drift_corrections!(joinpath(ancilary_data_path,"SolarRV20190918_JD_SciRV_CalRV.txt"), df_files_use)

NEID.read_barycentric_corrections!(joinpath(ancilary_data_path,"SolarTelescope2019-09-18_inclGravRedshiftAirMassAltitude.csv"), df_files_use)

RvSpectML.apply_doppler_boost!(solar_data,df_files_use)
