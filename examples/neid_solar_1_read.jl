using Pkg
Pkg.activate(".")
verbose = true
using Revise
if verbose   println("# Loading RvSpecML")    end
using RvSpectML
if verbose   println("# Loading other packages")    end
using DataFrames, Query

# TODO: USER:  Either the paths that specify where datafiles are stored here or in examples/data_paths.jl
if verbose   println("# Creating manifest of files to process.")    end
solar_data_path = "20190918"
ancilary_solar_data_path = "."
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

if verbose println("# Reading in FITS files.")  end
@time solar_data = map(NEID.read_solar_data,eachrow(df_files_use))

if verbose println("# Applying wavelength corrections.")  end
NEID.read_drift_corrections!(joinpath(ancilary_solar_data_path,"SolarRV20190918_JD_SciRV_CalRV.txt"), df_files_use)

NEID.read_barycentric_corrections!(joinpath(ancilary_solar_data_path,"SolarTelescope2019-09-18_inclGravRedshiftAirMassAltitude.csv"), df_files_use)

NEID.read_differential_extinctions!(joinpath(ancilary_solar_data_path,"20190918_diff_ex_full_fixed.txt"), df_files_use)

apply_doppler_boost!(solar_data,df_files_use)
