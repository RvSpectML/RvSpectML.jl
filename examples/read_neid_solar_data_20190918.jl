using Pkg|
 Pkg.activate(".")

verbose = true
if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML

# USER: You must create a data_paths.jl file in one of the default_paths_to_search listed below. It need only contain one line:
# solar_data_path = "/home/eford/Data/SolarSpectra/NEID_solar/"
target_subdir = "20190918"   # USER: Replace with directory of your choice
 fits_target_str = "Solar"
 output_dir = "examples/output"
 default_paths_to_search = [pwd(),"examples",joinpath(pkgdir(RvSpectML),"examples"), "/gpfs/group/ebf11/default/ebf11/neid_solar"]
 # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
 pipeline = PipelinePlan()

RvSpectML.Pipeline.reset_all_needs!(pipeline)
 if need_to(pipeline,:read_spectra)
   if verbose println("# Finding what data files are avaliable.")  end
   df_files = make_manifest(target_subdir, NEID )

   if verbose println("# Reading in customized parameters from param.jl.")  end
   eval(code_to_include_param_jl())

   if verbose println("# Reading in ", size(df_files,1), " FITS files.")  end
   @time all_spectra = map(NEID.read_solar_data,eachrow(df_files_use))
   dont_need_to!(pipeline,:read_spectra)

   if verbose println("# Applying wavelength corrections.")  end
   NEID.read_drift_corrections!(joinpath(ancilary_solar_data_path,"SolarRV20190918_JD_SciRV_CalRV.txt"), df_files_use)
   NEID.read_barycentric_corrections!(joinpath(ancilary_solar_data_path,"SolarTelescope2019-09-18_inclGravRedshiftAirMassAltitude.csv"), df_files_use)
   NEID.read_differential_extinctions!(joinpath(ancilary_solar_data_path,"20190918_diff_ex_full_fixed.txt"), df_files_use)
   apply_doppler_boost!(all_spectra,df_files_use)
   all_spectra
 end
