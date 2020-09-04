
fits_target_str = "101501"
num_max_spectra_to_use = 10

df_files_use = df_files |>
   @filter( _.target == fits_target_str ) |>
   @take(num_max_spectra_to_use) |>
   DataFrame

linelist_for_ccf_filename = "G2.espresso.mas"

ccf_mid_velocity = -5e3
tophap_ccf_mask_scale_factor=1.6
write_ccf_to_csv = false
write_ccf_to_csv = true
