# Default values shared across instruments
const default_line_width_kmps = predict_intrinsic_stellar_line_width(5780,v_rot=1.8)  # km/s
default_chunk_size_factor = 3       # For default_calc_chunk_width TODO: Figure out what value to use.  Ask Alex
default_calc_chunk_width = ChunkWidthFixedΔlnλ(default_chunk_size_factor*default_line_width_kmps*1000/speed_of_light_mps)  # Used by read_mask_espresso and read_vald for assigning lambda_lo and lambda_hi
default_min_chunk_Δv = 20           # km/s  for ChunkWidthFixedΔlnλ

default_Δv_to_avoid_tellurics = 15000.0  # m/s

 # TODO: OPT: Make above const once settle on good values
