# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:29:18 2013

@author: Michel Le Page

GRACE Total Water Content data was downloaded from http://www.esrl.noaa.gov/psd/data/gridded/data.gpcc.html
Three different solutions (GFZ, CSR, JPL) are given into netcdf files at a 1 degree resolution (360x180 pixels) for the period 2003 to present, at a monthly timestep. The matrix is oriented from 90S à 90N and 0 to 360 degrees

This code extract a sub-area and a sub-period, and save it to a TIF file

INPUT: infiles (netcdf GRACE files)
PARAMS:     year_begin, year_end, lat_sud, lat_nord, lon_ouest, lon_est,
OUTPUT: Three TIF files
"""

from math import *
import numpy as np
from osgeo import gdal
from osgeo.gdalconst import *
import sys
import scipy.ndimage


if __name__ == "__main__":

    # params des fichiers d'entree

    os.chdir("C:\\Users\\michel\\Documents\\GRACE\\med_dataset\\GRACE\\netcdf")

    infile_sf="CLM4.RL05.SCALE_FACTOR.DS.G200KM.nc"    
    inDs_sf=gdal.Open('NETCDF:"'+infile_sf+'":SCALE_FACTOR')
    cols = inDs_sf.RasterXSize
    rows = inDs_sf.RasterYSize
    bands = inDs_sf.RasterCount    
    bandi = inDs_sf.GetRasterBand(1)
    scale_factor=bandi.ReadAsArray(0,0, cols, rows)    
    inDs_sf=None
    
    infiles=["GRACE.CSR.LAND.RL05.DS.G200KM.nc","GRACE.GFZ.LAND.RL05.DS.G200KM.nc","GRACE.JPL.LAND.RL05.DS.G200KM.nc"]
    outfiles=["..\\GRACE.CSR.LAND.RL05.DS.G200KM_MEDIT_2003_2012.tif","..\\GRACE.GFZ.LAND.RL05.DS.G200KM_MEDIT_2003_2012.tif","..\\GRACE.JPL.LAND.RL05.DS.G200KM_MEDIT_2003_2012.tif"]
    rows=180
    cols=360
    resolution=1

    # facteur du resampling
    
    resize_factor=1
    
    # params de la zone à extraire
    year_begin=2003
    year_end=2013        
    lat_sud=24.
    lat_nord=48.
    lon_ouest=-15
    lon_est=40.

    # attention la matrice GRACE va de 90S à 90N et de 0 à 360°

    x0_ouest=int((360+lon_ouest)/resolution)
    x0_est=cols
    
    x1_ouest=0
    x1_est=int(lon_est/resolution)    

    y_sud=int((90-lat_sud)/resolution)
    y_nord=int((90-lat_nord)/resolution)

    rows_out=int((y_sud-y_nord)*resize_factor)
    cols_out=int(((x0_est-x0_ouest)+(x1_est-x1_ouest))*resize_factor)
    nbband_out=int((year_end-year_begin)*12)

    
    for ifile in range(0,3):
        infile=infiles[ifile]
        outfile=outfiles[ifile]
        print(infile)
        
        inDs=gdal.Open('NETCDF:"'+infile+'":lwe_thickness')
    
        cols = inDs.RasterXSize
        rows = inDs.RasterYSize
        bands =  inDs.RasterCount
    
        driver = gdal.GetDriverByName("GTiff")
        #outDs=driver.Create("GPCC_precip.mon.total.v6.tif",37+98,58,bands,GDT_Float32)
        outDs=driver.Create(outfile,cols_out,rows_out,nbband_out,GDT_Float32)    
        #outDs.SetGeoTransform(inDs.GetGeoTransform())
        #outDs.SetProjection(inDs.GetProjection())

        band_begin=(year_begin-2003)*12
        band_end=(year_end-2003)*12
    
        for band in range(band_begin,band_end):
            bandi=inDs.GetRasterBand(band+1)
            raster=bandi.ReadAsArray(0,0, cols, rows)
            
            raster=raster*scale_factor
    
            rasterout=np.zeros((rows_out/resize_factor,cols_out/resize_factor),np.float32)
         
            for j in range(y_nord,y_sud):
                for i in range(x0_ouest,x0_est):
                    rasterout[j-y_nord,i-x0_ouest]=raster[j,i]
                for i in range(x1_ouest,x1_est):
                    rasterout[j-y_nord,i+(x0_est-x0_ouest)]=raster[j,i]         
    
            bando=outDs.GetRasterBand(band-band_begin+1)
    
            if resize_factor > 1:
                raster_resample=scipy.ndimage.zoom(rasterout,resize_factor,order=3)
                bando.WriteArray( raster_resample )    
            if resize_factor < 1:
                raster_resample=block_reduce(rasterout, block_size=(1./resize_factor, 1./resize_factor), func=np.mean)
                bando.WriteArray( raster_resample )                   
            if resize_factor == 1:
                bando.WriteArray( rasterout ) 
            

        outDs=None
     
     