# ------------------------------------------------------------------
# Compute Reference Evapotranspiration ET0 on a daily time step
# Michel Le Page, 2011
#
# The FAO-56 Penman-Monteith equation (Allen et al., 1998) estimates the reference crop evapotranspiration, ET0 (mm d-1)
#
# Monteith, J.L., 1965. Evaporation and Environment. 19th Symposia of the Society for Experimental Biology, University Press, Cambridge, 19:205-234.
#
# day1: julian day from 0 to 366
# alt : altitude of the meteo station (meters)
# hmes: height of measurement (meters)
# lat : latitude of the meteo station (decimal degrees)
# tmoy: average temperature of the day (Celsius Degrees)
# tmin: minimum temperature of the day (Celsius Degrees)
# tmax: maximum temperature of the day (Celsius Degrees)
# vv: Average Wind Speed of the day (m/s)
# hrmoy: average Relative Humidity of the day (%)
# hrmin: minimum Relative Humidity of the day (%)
# hrmax: maximum Relative Humidity of the day (%)
# rs: Solar Radiation averaged for the day (W/m2)
# 
# ------------------------------------------------------------------

def et0_pm(day1,alt,hmes,lat,tmoy,tmin,tmax,vv,hrmoy,hrmin,hrmax,rs):
  if ((hrmax-hrmin)<.1 or tmin==tmax): return(-1)
  try: 
    conv_rad=lat * pi / 180.
    # -- constantes --
    const_sol=0.082
    const_stef=4.9*pow(10,-9)
    g=0.063
 
    u2= vv * 4.87 / log(67.8*hmes-5.42)
    rs_mj= rs*24*3600*0.000001
    dr= 1+0.033*cos(2*pi*day1/365)
    d= 0.409*sin((2*pi*day1/365)-1.39)
    ws= acos(-tan(conv_rad)*tan(d))
    ra= (24*60/pi)*const_sol*dr*(ws*sin(conv_rad)*sin(d)+cos(conv_rad)*cos(d)*sin(ws))
    rso= ra*(0.75+0.00002*alt)
    f= 1.35*(rs_mj/rso)-0.35
    es= (0.6108*exp(17.27*tmin/(tmin+237.3))+0.6108*exp(17.27*tmax/(tmax+237.3)))/2
    ea= (hrmin*0.6108*exp(17.27*tmax/(tmax+237.3))+hrmax*0.6108*exp(17.27*tmin/(tmin+237.3)))/(2*100)
    rnl= const_stef*0.5*(pow((tmin+273),4)+pow((tmax+273),4))*(0.34-0.14*sqrt(ea))*f
    rn=(1-0.23)*rs_mj-rnl
    
    gflux=0
    delta= 4098*0.6108*exp(17.27*tmoy /(tmoy +237.3))/pow((tmoy +237.3),2)
    et0= (0.408*delta*(rn)+(900*g/(tmin+273))*u2*(es-ea))/(delta+g*(1+0.34*u2))
    return(et0)
  except:
    return(-1)