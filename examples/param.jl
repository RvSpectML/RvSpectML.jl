linelist_for_ccf_filename = "G2.espresso.mas"
ccf_mid_velocity = 0
tophap_ccf_mask_scale_factor=1.6

max_spectra_to_use = 200
if fits_target_str == "Solar"
   ancilary_solar_data_path = "/home/eford/Data/SolarSpectra/NEID_solar/"
   bjd_first_good = 2458745.1296134139
   bjd_last_good = 2458745.283
   df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      @take(max_spectra_to_use) |>
      DataFrame
elseif fits_target_str == "101501"
   ccf_mid_velocity = -5e3
   df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @take(max_spectra_to_use) |>
      DataFrame
end
