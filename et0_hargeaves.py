# ------------------------------------------------------------------
# Compute Reference Evapotranspiration ET0 on a daily time step
# Michel Le Page, 2011
#
# When solar radiation data, relative humidity data and/or wind speed data are missing, they should be estimated using the procedures presented in this section. As an alternative, ETo can be estimated using the Hargreaves ETo equation
#
# Hargreaves, G.L., Hargreaves, G.H., and Riley, J.P. 1985. Agricultural benefits for Senegal River Basin. J. Irrigation and Drainage Engr., ASCE 111:113-124.
#
# day1: julian day from 0 to 366
# alt : altitude of the meteo station (meters)
# tmoy: average temperature of the day (Celsius Degrees)
# tmin: minimum temperature of the day (Celsius Degrees)
# tmax: maximum temperature of the day (Celsius Degrees)
#
# ------------------------------------------------------------------

def et0_hargreaves(day1,lat,tmoy,tmin,tmax):
	if (tmin==tmax): return(-1)
	try:
		conv_rad=lat * 3.1416 / 180
		const_sol=0.082
		p=3.142
		dr= 1+0.033*cos(2*pi*day1/365)
		d= 0.409*sin((2*pi*day1/365)-1.39)
		ws= acos(-tan(conv_rad)*tan(d))
		ra= (24*60/pi)*const_sol*dr*(ws*sin(conv_rad)*sin(d)+cos(conv_rad)*cos(d)*sin(ws))
		ra=ra/2.45
		et0=0.0023*(tmoy+17.8)*sqrt(tmax-tmin)*ra
		return(et0)
	except:
		return(-1)