#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-

# Michel Le Page, 08/04/2013
#
# This code calculate some simple stats from an image file (in this case a DEM) related to a shapefile
# input: a rasterfile (repert and infile) and a shapefile (shapefile) with its ID column (attribut)
# output: a text file with, id,number,average,stdev,min,max
#


import os
from math import *
import numpy as np
from osgeo import gdal
import csv
from osgeo.gdalconst import *
from osgeo import gdal, gdal_array
import sys

if __name__ == '__main__':

    gdaldem="C:\\OSGeo4W\\bin\\gdaldem"
    gdal_rasterize="C:\\OSGeo4W\\bin\\gdal_rasterize"

    repert="D:\\DEM\\500m_SRTM (CGIAR_2008) Global et Maroc"
    infile="srtm_500m_maroc.tif"

    shapefile="d:\\shp\\communes_maroc_geo.shp"    
    attribut="COMMUNE_ID"
    layer="communes_maroc_geo"

    os.chdir(repert)

    inDs = gdal.Open(infile)
    cols = inDs.RasterXSize
    rows = inDs.RasterYSize
    bands =  inDs.RasterCount
    driver = inDs.GetDriver()
    geo=inDs.GetGeoTransform()
    print geo
    proj=inDs.GetProjection()
    print proj
    band=inDs.GetRasterBand(1)
    maxi=band.GetMaximum()
    mini=band.GetMinimum()
    inDs=0
    

    # sauvergarder la projection
    #f = open('proj.txt', 'w')
    #f.write(proj)
    #f.close()
    # rasteriser le fichier shape
    commande=gdal_rasterize+" -a "+attribut+" -ts "+str(cols)+" "+str(rows)
    commande+=" -l "+layer+" "+shapefile+" tmp.tif -ot Int16"
    #commande+=" -a_srs proj.txt " 
    commande+=" -te "+str(geo[0])+" "+str(geo[3]+geo[5]*rows)+" "+str(geo[0]+geo[1]*cols)+" "+str(geo[3])

    print commande

    #os.system(commande)

    inDs = gdal.Open(infile)
    band=inDs.GetRasterBand(1)
    mnt=band.ReadAsArray()

    inDs2 = gdal.Open("tmp.tif")
    driver2 = inDs2.GetDriver()
    band2=inDs2.GetRasterBand(1)
    shp=band2.ReadAsArray()
    maxishp=shp.max()
    minishp=shp.min()   
    print minishp,maxishp
    
    f = open('stats.txt', 'w')
    out= "i,n,moy,std,min,max\n"
    f.write(out)
    
    aa=np.unique(shp)

    for i in aa:
        a=np.where(shp==i)
        if mnt[a].size>0:
            moy=mnt[a].mean()
            std=mnt[a].std()
            min=mnt[a].min()
            max=mnt[a].max()
            out= str(i)+","+str(mnt[a].size)+","+str(moy)+","+str(std)+","+str(min)+","+str(max)+"\n"
            f.write(out)
        if (i%100)==0:print(i)
            
    f.close()

