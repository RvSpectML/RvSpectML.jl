function calc_solar_alt(time::Real)
  obs = AstroLib.observatories["kpno"]
   (ra_sun, dec_sun, ) = sunpos(time)
   lst = ct2lst(obs.longitude,time)
   (alt, az) = hadec2altaz(lst-ra_sun/360*24, dec_sun,obs.latitude)
   return alt
end
