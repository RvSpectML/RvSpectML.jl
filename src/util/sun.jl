
""" get_obs_loc(obs::Symbol)
 Returns a Dict with long & lat (degrees) and elevation (km)
 Warning: Currently only has info for :WIYN and :HARPSN.
"""
function get_obs_loc(obs::Symbol)
	valid_obs = [:WIYN, :HARPSN]
	@assert obs ∈ valid_obs
	if obs == :HARPSN
		loc = Dict("lon"=> -17.88905, "lat"=> 28.754, "elevation"=> 2.3872)
	elseif obs == :WIYN
		loc = Dict("lon"=> -111.600562, "lat"=> 31.958092, "elevation"=> 2.091)
	else
		@error("Don't have coordinates for obs = " * string(obs))
	end
	return loc
end

function calc_solar_alt(time::Real; obs::Symbol)
	loc = get_obs_loc(obs)
  #=
	obs = AstroLib.observatories["kpno"]
	lat = obs.latitude
	long = obs.longitude
	=#
  #=
	 Values for specifically for WIYN
  long = -111.600562
  lat =  31.958092
  alt =  2091.0
	=#
  (ra_sun, dec_sun, ) = sunpos(time)
  lst = ct2lst(loc["lon"],time)
  (alt, az) = hadec2altaz(lst-ra_sun/360*24, dec_sun,loc["lat"])
  return alt
end
