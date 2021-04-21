function calc_solar_alt(time::Real; loc::Symbol)
  #=
	obs = AstroLib.observatories["kpno"]
	lat = obs.latitude
	long = obs.longitude
	=#
  # Values for specifically for WIYN
  long = -111.600562
  lat =  31.958092
  alt =  2091.0
  (ra_sun, dec_sun, ) = sunpos(time)
  lst = ct2lst(long,time)
  (alt, az) = hadec2altaz(lst-ra_sun/360*24, dec_sun,lat)
  return alt
end
