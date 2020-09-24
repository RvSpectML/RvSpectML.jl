verbose = true
 if verbose && !isdefined(Main,:RvSpectML)   println("# Loading RvSpecML")    end
 using RvSpectML

# USER: You must create a data_paths.jl file in one of the default_paths_to_search listed below. It need only contain one line:
# expres_data_path = "/path/to/EXPRES/data/not/including/target_subdir"
target_subdir = "101501"   # USER: Replace with directory of your choice
 fits_target_str = "101501"
 output_dir = "examples/output/"
 default_paths_to_search = [pwd(),"examples",joinpath(pkgdir(RvSpectML),"examples"),"/gpfs/group/ebf11/default/ebf11/expres/inputs"]
 # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
 pipeline_plan = PipelinePlan()
 dont_make_plot!(pipeline_plan, :movie)

reset_all_needs!(pipeline_plan)
 if need_to(pipeline_plan,:read_spectra)
   if verbose println("# Finding what data files are avaliable.")  end
   df_files = make_manifest(target_subdir, EXPRES )

   if verbose println("# Reading in customized parameters from param.jl.")  end
   eval(code_to_include_param_jl())

   if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = map(EXPRES.read_data,eachrow(df_files_use))
   RvSpectML.discard_blaze(all_spectra)
   RvSpectML.discard_continuum(all_spectra)
   dont_need_to!(pipeline_plan,:read_spectra)
   #all_spectra
 end
