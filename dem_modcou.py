#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-

# Michel Le Page, 06/04/2013
#
# Create a unregular matrix from a DEM file for Modcou model
# Input: A DEM file with a pixel size wihch is multiple of 1000 meters (eg 500, 250, 125, 100...)
# Outputs: Two filse. One is an index file that can be used within a GIS app' to be vectorize. The second one if the average height of each new grid.

from math import *
import numpy as np
#import scipy as sp
from osgeo import gdal
import csv
from osgeo.gdalconst import *


if __name__ == "__main__":

	# Fichiers
	infile="D:\DEM\GMTED2010\Gmted_250m_maroc.tif"
	outfile1="D:\DEM\GMTED2010\Gmted_250m_maroc_surfexind.tif"
	outfile2="D:\DEM\GMTED2010\Gmted_250m_maroc_surfexval.tif"
	
	resol=250 #resolution en metres du DEM d'entree, attention, doit etre un multiple de 1000m (500, 250, 125, 100...)
	
	pasmax=8000/resol # pas en pixel pour 8km
	pasmin=1000/resol # pas en pixel pour 1km
	stdevmin=30 # min stdev pour regrouper des pixels
	resols=(pasmax/8,pasmax/4,pasmax/2,pasmax) # résolutions intermediaires à 1, 2, 4 et 8 km, attention l'ordre du + petit au + grand est important

	inDs = gdal.Open(infile)
	cols = inDs.RasterXSize
	rows = inDs.RasterYSize
	bands =  inDs.RasterCount
	driver = inDs.GetDriver()

	print "Cols:", cols, "Ligs:", rows, "Bands:", bands,driver


	band=inDs.GetRasterBand(1)
	  
	mnt=band.ReadAsArray()
	
	matout=np.zeros([rows,cols], np.uint32)
	mntout=np.zeros([rows/pasmin,cols/pasmin], np.float32)
	
	val=1L

	for i in range(0,rows,pasmax):
		if i%256==0: print i
		i1=i/pasmin
		for j in range(0,cols,pasmax):
			j1=j/pasmin
			# pour les 4 résolutions intermédiaires, on teste de la plus petite à la plus grande resol
			for k in resols:
				k1=k/pasmin
				for ii in range(0,pasmax,k):
					ii1=ii/pasmin
					for jj in range(0,pasmax, k):
						stdev=np.std(mnt[i+ii:i+ii+k,j+jj:j+jj+k])
						if stdev<stdevmin or (k==pasmin and stdev>=stdevmin):
							#print "ok:", (i*pasmax)+ii,(i*pasmax)+ii+k,(j*pasmax)+jj,(j*pasmax)+jj+k, val
							matout[i+ii:i+ii+k,j+jj:j+jj+k]=val
							mean=np.mean(mnt[i+ii:i+ii+k,j+jj:j+jj+k])
							mntout[i1+ii1:i1+ii1+k1,j1+jj/pasmin:j1+jj/pasmin+k1]=mean
							val=val+1
						#else:
                                                        #        print i,j,stdev, mnt[i+ii:i+ii+k,j+jj:j+jj+k]
                                                        
		
	outDs=driver.Create(outfile1,cols,rows,bands,GDT_UInt32)
	outDs.SetGeoTransform(inDs.GetGeoTransform())
	outDs.SetProjection(inDs.GetProjection())
	bando=outDs.GetRasterBand(1)
	bando.WriteArray(matout)
	outDs=None

	outDs=driver.Create(outfile2,cols/pasmin,rows/pasmin,bands,GDT_Float32)
	geo=inDs.GetGeoTransform()
	# changement de la résolution de sortie à 1km
	lgeo=list(geo)
	lgeo[1]*=pasmin
	lgeo[5]*=pasmin
	tgeo=tuple(lgeo)
	outDs.SetGeoTransform(tgeo)
	outDs.SetProjection(inDs.GetProjection())
	bando=outDs.GetRasterBand(1)
	bando.WriteArray(mntout)
	outDs=None
	
	inDs=None


