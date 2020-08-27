# Default values shared across instruments
default_chunk_size_factor = 3       # TODO: Figure out what value to use.  Ask Alex
default_min_chunk_Δv = 20          # km/s
default_line_width_kmps = predict_line_width(5780,v_rot=1.8)  # km/s
default_calc_chunk_width = ChunkWidthFixedΔlnλ(default_chunk_size_factor*default_line_width_kmps*1000/speed_of_light_mps)
