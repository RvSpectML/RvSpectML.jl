global linelist_for_ccf_filename = "G2.espresso.mas"
global ccf_mid_velocity = 0
global tophap_ccf_mask_scale_factor=1.6

global max_spectra_to_use = 20
global fits_target_str 
if fits_target_str == "Solar"
   global ancilary_solar_data_path = "/home/eford/Data/SolarSpectra/NEID_solar/"
   global bjd_first_good = 2458745.1296134139
   global bjd_last_good = 2458745.283
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      @take(max_spectra_to_use) |>
      DataFrame
elseif fits_target_str == "101501"
   global ccf_mid_velocity = -5e3
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @take(max_spectra_to_use) |>
      DataFrame
end
