"""
Convenience code to deal with finding files, reading parameters, etc.

Author: Eric Ford
Created: Sept 2020
"""

const default_paths_to_search = [pwd(),"examples",joinpath(pkgdir(RvSpectML),"examples"),"/gpfs/group/ebf11/default/ebf11/expres/inputs"]

""" make_manifest(target_subdir::String, Inst::Module; [opts] )
Returns a dataframe containing a list of files to be read and some metadata (e.g., observation times)

# Optional arguements
- max_spectra_to_use (default_max_spectra_to_use)
- paths_to_search (default_paths_to_search)
- verbose = true

Warning:  Malicious users could insert arbitrary code into data_paths.jl.  Don't be a malicous user.
"""
function make_manifest(target_subdir::String, Inst::Module; paths_to_search::Union{String,AbstractVector{String}} = default_paths_to_search,
      # max_spectra_to_use::Integer = 1000, target_fits::String = target_subdir,
      verbose::Bool = true)

   if verbose   println("# Looking for data paths and config files in default_paths_to_search.")    end
   idx_path = findfirst(isfile,map(d->joinpath(d,"data_paths.jl"),paths_to_search))
   if isnothing(idx_path)
      @error("Can't find data_paths.jl in any of: ", paths_to_search)
   end
   data_paths_jl = paths_to_search[idx_path]
   println("# Found ", data_paths_jl)
   flush(stdout)
   if !isdir(data_paths_jl)
         @error("Can't access instrument's base data directory ", data_paths_jl, ".")
   end

   if isfile(joinpath(data_paths_jl,"data_paths.jl"))
       include(joinpath(pwd(),data_paths_jl,"data_paths.jl"))
   end

   if Inst == RvSpectML.EXPRES
      #@assert isdefined(Main,:expres_data_path)
      data_path = joinpath(expres_data_path,target_subdir)
   elseif Inst == RvSpectML.NEID
      #  @assert isdefined(Main,:solar_data_path)
      data_path = joinpath(solar_data_path,target_subdir)
   else
      @error "Specified instrument didn't match a module name, ", string(Inst), "."
   end
   if !isdir(data_path)
      @error "Can't access target data directory ", data_path, "."
   end

   if verbose  println("# Creating manifest of files to process.")    end
   df_files = Inst.make_manifest(data_path)
   if size(df_files,1) < 1
      @error("Did not find any files in ", data_path ,".")
   end
   #=
   df_files = df_files |> @filter( _.target == target_fits ) |> DataFrame
   if size(df_files,1) < 1
      @error("Did not find any files with target field matching ", target_fits,".")
   end
   if size(df_files,1) > max_spectra_to_use
      df_files = df_files |> @take(max_spectra_to_use) |> DataFrame
   end
   =#
   return df_files
end

"""   code_toread_param_jl( path_to_search )

Returns a Code object.  After `res = code_toread_param_jl( path_to_search )`,
execute `eval(res)` to actually include the param.jl file.
This is useful since it allows variables to be placed into caller's namespace.

Warning:  Malicious users could insert arbitrary code into param.jl.  Don't be a malicous user.
"""
function code_to_include_param_jl(paths_to_search::Union{String,AbstractVector{String}} = default_paths_to_search; filename::String = "param.jl", verbose::Bool = true)
   if verbose   println("# Looking for param.jl file to set configuration parameters.")    end
   idx_path = findfirst(isfile,map(d->joinpath(d,filename),default_paths_to_search))
   if isnothing(idx_path)
      @error(" Couldn't find $filename in $default_paths_to_search")
      return Expr()
   end
   data_path = default_paths_to_search[idx_path]

   if isfile(joinpath(data_path,filename))
      println("# Reading parameter values from ", filename)
      code_to_include_param = quote
         include(joinpath(pwd(),$data_path,$filename))
      end
   else
      println("# Did not locate ", filename, " in ",data_path, ".  Returning code to set default values for lots of things.")
      # Set defaults in case not in param.jl
      code_to_include_param = quote
         df_files_use = df_files
         espresso_filename = joinpath(pkgdir(RvSpectML),"data","masks","G2.espresso.mas")
         ccf_mid_velocity = 0
         tophap_ccf_mask_scale_factor=1.6
      end
   end
   return code_to_include_param
end
