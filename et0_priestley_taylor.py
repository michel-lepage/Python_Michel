# ------------------------------------------------------------------
# Compute Reference Evapotranspiration ET0 with Priestley Taylor method on a daily time step
# Michel Le Page, 2011

# The Priestley-Taylor method (Priestley and Taylor, 1972) for the calculation of daily ET0 (mm d-1) replaces the aerodynamic term of Penman-Monteith equation by a dimensionless empirical multiplier (a, Priestley-Taylor coefficient)
#
# day1: julian day from 0 to 366
# alt : altitude of the meteo station (meters)
# lat : latitude of the meteo station (decimal degrees)
# tmoy: average temperature of the day (Celsius Degrees)
# tmin: minimum temperature of the day (Celsius Degrees)
# tmax: maximum temperature of the day (Celsius Degrees)
# hrmin: minimum Relative Humidity of the day (%)
# hrmax: maximum Relative Humidity of the day (%)
# rs: Solar Radiation averaged for the day (W/m2)
# 
# ------------------------------------------------------------------

from math import *

def et0_priestley_taylor(day1,alt,lat,tmoy,tmin,tmax,hrmin,hrmax,rs):
	if ((hrmax-hrmin)<.1 or tmin==tmax): return(-1)
	try: 
		lambda=2.45
		conv_rad=lat * 3.1416 / 180
		pt_coeff=1.26
		const_sol=0.082
		g=0.063
		gflux=0
		p=3.142
		const_stef=4.9*pow(10,-9)

		rs_mj= rs*24*3600*0.000001
		dr= 1+0.033*cos(2*p*day1/365)
		d= 0.409*sin((2*p*day1/365)-1.39)
		ws= acos(-tan(conv_rad)*tan(d))
		ra= (24*60/p)*const_sol*dr*(ws*sin(conv_rad)*sin(d)+cos(conv_rad)*cos(d)*sin(ws))
		rso= ra*(0.75+0.00002*alt)
		ea= (hrmin*0.6108*e(17.27*tmax/(tmax+237.3))+hrmax*0.6108*e(17.27*tmin/(tmin+237.3)))/(2*100)
		f= 1.35*(rs_mj/rso)-0.35
		rnl= const_stef*0.5*(pow((tmin+273),4)+pow((tmax+273),4))*(0.34-0.14*sqrt(ea))*f
		rn= (1-0.23)*rs_mj-rnl
		delta= 4098*0.6108*e(17.27*tmoy /(tmoy +237.3))/pow((tmoy +237.3),2)

		et0=pt_coeff*(delta/(delta+g))*(rn/lambda)
		return(et0)
	except:
		return(-1)