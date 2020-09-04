using Pkg
Pkg.activate(".")

verbose = true
using Revise
if verbose   println("# Loading RvSpecML")    end
using RvSpectML
if verbose   println("# Loading other packages")    end
using DataFrames, Query

expres_data_path = "."
 target_subdir = "101501"

 # TODO: USER:  Either the paths that specify where datafiles are stored here or in examples/data_paths.jl
 if verbose   println("# Creating manifest of files to process.")    end
 if isdir("examples")
    cd("examples")
 end
 if isfile("data_paths.jl")
      include("data_paths.jl")
 end
 expres_data_path

num_spectra_to_use = 20
 #bjd_first_good = 2458745.1296134139
 #bjd_last_good = 2458745.283
 df_files = EXPRES.make_manifest(joinpath(expres_data_path,target_subdir))
 df_files_use = df_files |>
   @filter( _.target == "101501" ) |>
   # @take(num_spectra_to_use) |>
   DataFrame

if verbose println("# Reading in FITS files.")  end
@time expres_data = map(EXPRES.read_data,eachrow(df_files_use))
