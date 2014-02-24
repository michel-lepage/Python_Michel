#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-

# Michel Le Page, 20/04/2013
# Winstral, Spatial Snow Modeling of Wind-Redistributed Snow Using Terrain-Based Parameters
# Journal of Hydrometeorology, 2002
# seules les eq 1 et 2 sont calculées.
# Input : un DEM
# Output: average of n maximum upwind slopes relative to seasonally averaged winds (TIF Float16)
#
# Python 3.3


from math import *
import numpy as np
from osgeo import gdal
from osgeo.gdalconst import *
import sys


if __name__ == "__main__":

        # ======= DEBUT CONFIG ==============
    # Fichiers
    infile="testDEM.tif"
    outfile1="testDEM_winstral.tif"
    
    resol=30.               #resolution en metres du DEM d'entree

    azimuth_vent=90.        # azimuth du vent dominant en degrés
    az1=azimuth_vent-15.    #azimuth dominant - un certain angle (-15° par défaut)
    az2=azimuth_vent+15.    #azimuth dominant + un certain angle (+15° par défaut)
    inc=5.                  # incrément de calcul des droites de Vent (5° par défaut)
    dmax=500                # distance dmax

    # ======= FIN CONFIG ==============

    nbvect=int(((az2-az1)/inc)+1)   # Nombre de vecteurs= nombre de directions de calcul (fig3, p528)
    longvect=int(dmax/resol)        # Longueur des vecteurs en fonction du param dmax.

    # lecture du fichier d'entree
    inDs = gdal.Open(infile)
    cols = inDs.RasterXSize
    rows = inDs.RasterYSize
    bands =  inDs.RasterCount
    driver = inDs.GetDriver()
    print("Cols:", cols, "Ligs:", rows, "Bands:", bands,driver)
    band=inDs.GetRasterBand(1)
    mnt=band.ReadAsArray()
    mnt=np.float32(mnt)

    # on passe tout en radians
    azimuth_vent=azimuth_vent/180*pi 
    az1=az1/180.*pi
    az2=az2/180.*pi
    inc=inc/180.*pi

    # vect[][0]: les x, vect[][1]=mles y, vect[][2]=les dists
    vects=np.zeros((nbvect,3,longvect),np.float32)

    # Calcul des coordonnées pixel des droites.

    for n in range(0,nbvect):
        az=az1+(n*inc)
        print("Azimuth=",az/pi*180)
        x=dmax*cos(az)
        y=dmax*sin(az)
        pixx=x/dmax
        pixy=y/dmax
        print(x,y,pixx,pixy)

        for xx in range(1,longvect+1):
            vects[n][0][xx-1]=round(xx*pixx)
            vects[n][1][xx-1]=round(xx*pixy)
            vects[n][2][xx-1]=sqrt(pow(float(vects[n][0][xx-1])*resol,2)+pow(float(vects[n][1][xx-1])*resol,2))
            print(vects[n][0][xx-1],vects[n][1][xx-1],vects[n][2][xx-1])

    maxx=max(0,int(np.max(vects[:,0,:]))+1)
    maxy=max(0,int(np.max(vects[:,1,:]))+1)
    minx=max(0,(int(np.min(vects[:,0,:]))-1)*-1)
    miny=max(0,(int(np.min(vects[:,1,:]))-1)*-1)
    print(minx,miny,maxx,maxy)
    
    # matrice de sortie
    mntout=np.zeros((rows,cols),np.float32)

    # stockage temporaire des pentes
    sx=np.zeros(nbvect)
    
    #calculs pour chaque pixel
    hash10=range(0,rows,int(rows/10))
    hash100=range(0,rows,int(rows/100))
    maxyy=rows-maxy-1
    maxxx=cols-maxx-1
    for j in range(miny,maxyy,1):
        if j in hash10: print(hash10.index(j)*10,'%',end="")
        if j in hash100: print('.',end="")
        sys.stdout.flush()
        sx=sx*0.
        for i in range(minx,maxxx,1):
            z=mnt[j,i]
            # calcule les sx pour les 7 (nbvect) droites.
            for n in range(0,nbvect):
                alt=mnt[np.int32(vects[n][1]+j),np.int32(vects[n][0]+i)] # numérateur de eq1 p528
                sx[n]=np.max(np.tan((alt-z)/vects[n][2])) # eq 1 p528
            mntout[j,i]=np.average(sx) # eq 2 p 529 (ou comment écrire quelque chose de simple de manière compliquée!)
    
    mntout=np.arctan(mntout)/pi*180
    #ecriture de mnt out
    outDs=driver.Create(outfile1,cols,rows,bands,GDT_Float32)
    outDs.SetGeoTransform(inDs.GetGeoTransform())
    outDs.SetProjection(inDs.GetProjection())
    bando=outDs.GetRasterBand(1)
    bando.WriteArray(mntout)
    outDs=None

    # done!        print "done!"


