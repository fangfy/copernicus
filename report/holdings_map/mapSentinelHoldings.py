#!/usr/bin/env python
"""
Generate a raster image showing where the AusCopHub server holds Sentinel-2 imagery. 
Uses a list of all of the AusCopHub XML files directly from the server, and extracts the
coordinates of the footprint polygons. Using just the vertices in those, create a
raster image which has 1 in every pixel whas has any footprint vertices in it. Note that it
does not count them, just records that at least one has occurred. 

"""
from __future__ import print_function, division

import sys
import os
import argparse
import fnmatch

import numpy
from osgeo import ogr, gdal, osr

from auscophub import auscophubmeta

def getCmdargs():
    """
    Get commandline arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("--outfile", default="sentinelholdings.img",
        help="Name of output image file (default=%(default)s)")
    p.add_argument("--searchdir", default="/g/data/fj7/Copernicus/Sentinel-2",
        help="Directory from which to search for XML files (default=%(default)s)")
    cmdargs = p.parse_args()
    
    if cmdargs.outfile is None:
        print("\n\nRequire --outfile argument\n\n")
        p.print_help()
        sys.exit()
    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    # Walk through the subdirectories of search dir, and get a list of all "*.xml" files.     
    xmlfilelist = getXmlList(cmdargs)
    if len(xmlfilelist) == 0:
        print("No XML files found")
        sys.exit()
    
    # Make a list of every vertex of every footprint. 
    pointsList = []
    for xmlfile in xmlfilelist:
        info = auscophubmeta.AusCopHubMeta(filename=xmlfile)
        geom = ogr.Geometry(wkt=str(info.footprintWkt))
        coords = eval(geom.ExportToJson())['coordinates']
        coords = coords[0]
        pointsList.extend(coords)
    
    # A blank raster image array
    img = numpy.zeros((180, 360), dtype=numpy.uint8)
    
    # Round the vertices to whole degrees
    pointsArr = numpy.array(pointsList).round().astype(numpy.int32)
    
    # Turn the whole degree values into row/col coords in the array
    row = (179 - (pointsArr[:, 1] + 90).clip(0, 179))
    col = (pointsArr[:, 0] + 180).clip(0, 359)
    
    # For every point, set the corresponding pixel to 1
    img[row, col] = 1
    
    # Write the raster array into an image file. 
    drvr = gdal.GetDriverByName('HFA')
    if os.path.exists(cmdargs.outfile):
        drvr.Delete(cmdargs.outfile)
        
    ds = drvr.Create(cmdargs.outfile, 360, 180, 1, gdal.GDT_Byte, ['COMPRESS=YES'])
    band = ds.GetRasterBand(1)
    band.WriteArray(img)
    band.SetNoDataValue(0)
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)
    ds.SetProjection(sr.ExportToWkt())
    ds.SetGeoTransform((-180, 1, 0, 90, 0, -1))
    del ds


def getXmlList(cmdargs):
    """
    Walk the directory structure to find all files matching "*.xml". 
    """
    xmlfilelist = []
    for (dirname, subdirs, filenames) in os.walk(cmdargs.searchdir):
        for filename in filenames:
            if fnmatch.fnmatch(filename, "*.xml"):
                if '.png' in filename: continue
                fullXmlFilename = os.path.join(dirname, filename)
                xmlfilelist.append(fullXmlFilename)
    return xmlfilelist


if __name__ == "__main__":
    main()

