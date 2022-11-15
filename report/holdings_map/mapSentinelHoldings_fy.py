#!/usr/bin/env python
"""
Generate a daily time-lapse showing where the AusCopHub server holds Sentinel-1/2/3 imagery. 
Uses a list of all of the AusCopHub XML files directly from the server, and extracts the
coordinates of the footprint polygons. Counting the center locations, create
raster images at 1 degree resolution. These images are then compiled into a time-lapse. 

Adapted from Neil Flood's mapSentinelHoldings.py.

"""
from __future__ import print_function, division

import sys
import os
import argparse
import fnmatch
import datetime

import numpy
from osgeo import ogr, gdal, osr

from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import colors

from images2gif import writeGif
from PIL import Image

from auscophub import auscophubmeta, geomutils


def getCmdargs():
    """
    Get commandline arguments
    """
    p = argparse.ArgumentParser()
    #p.add_argument("--outfile", default="sentinelholdings.png",
    #    help="Name of output image file (default=%(default)s)")
    p.add_argument("--collection", default="Sentinel-2",
                   help="Which collection to search for XML files (default=%(default)s)")
    p.add_argument("--product", default="L1C",
                   help="Which product type to search for XML files (default=%(default)s)")
    p.add_argument("--startdate", default=(datetime.datetime.today()-datetime.timedelta(90)).strftime('%Y%m%d'),
                   help="Start date of search (default=%(default)s)")
    p.add_argument("--enddate", default=(datetime.datetime.today()-datetime.timedelta(1)).strftime('%Y%m%d'),
                   help="End date of search (default=%(default)s)")
    p.add_argument("--remake", action='store_true',
                   help="Remake a png even if already exists")
    p.add_argument("--outdir", default="/short/fj7/fxy120/maps/",
                   help="Directory for output maps (default=%(default)s)")
    cmdargs = p.parse_args()
    
    #if cmdargs.outfile is None:
    #    print("\n\nRequire --outfile argument\n\n")
    #    p.print_help()
    #    sys.exit()
    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()
    if not os.path.exists(cmdargs.outdir):
        os.mkdir(cmdargs.outdir)

    startdate = datetime.datetime.strptime(cmdargs.startdate,'%Y%m%d')
    enddate = datetime.datetime.strptime(cmdargs.enddate,'%Y%m%d')
    
    fullxmllist=[]
    pnglist=[]
    #loopthrough every day
    d=startdate
    while d<=enddate:
        print (d)
        #reset xmlfilelist
        xmlfilelist=[]
        #define search pattern
        searchdir, pattern = formSearchDir(cmdargs.collection,cmdargs.product,d)
        #print (searchdir,pattern)
        #find xml files
        if len(searchdir)>0:
            xmlfilelist=getXmlList(searchdir, pattern)
        #plot for this date
        if len(xmlfilelist)>0:
            fullxmllist+=xmlfilelist
            #title='%s_'%cmdargs.collection+datetime.datetime.strftime(d,'%Y%m%d')
            #if cmdargs.remake or (not os.path.exists(cmdargs.outdir+'/'+title+'.png')):
            #    print ("making png")
            #    pngname=plotXml(xmlfilelist, title,outputdir=cmdargs.outdir)
            #if os.path.exists(cmdargs.outdir+'/'+title+'.png'):
            #    pnglist.append(cmdargs.outdir+'/'+title+'.png')
        d+=datetime.timedelta(1)
        
    #finally make the full coverage
    title='%s_'%cmdargs.collection+datetime.datetime.strftime(startdate,'%Y%m%d')+'_'+datetime.datetime.strftime(enddate,'%Y%m%d')
    pngname=plotXml(fullxmllist, title, outputdir=cmdargs.outdir)
    if pngname: pnglist.append(pngname)

    #now make a movie
    if len(pnglist)>1:
        images=[Image.open(fn) for fn in pnglist]
        writeGif(cmdargs.outdir+'/'+title+'.gif',images,duration=1.,repeat=False)

    
def plotXml(xmlfilelist, title, outputdir='.'):
    # Make a list of every vertex of every footprint. 
    pointsList = []
    for xmlfile in xmlfilelist:
        info = auscophubmeta.AusCopHubMeta(filename=xmlfile)
        geom = ogr.Geometry(wkt=str(info.footprintWkt))
        # find centroid
        prefEpsg = geomutils. findSensibleProjection(geom)
        centroidXY = geomutils.findCentroid(geom, prefEpsg)
        pointsList.extend([centroidXY])
        ## find vertices
        #coords = eval(geom.ExportToJson())['coordinates']
        #coords = coords[0]
        #pointsList.extend(coords)
    
    if len(pointsList)==0: return None
    # Round the vertices/centroids to whole degrees
    pointsArr = numpy.array(pointsList).round().astype(numpy.int32)
    try:
        west=numpy.where(pointsArr[:,0]<0)[0]
        pointsArr[west,:]=pointsArr[west,:]+[360.,0]
    except:
        pass

    #find out boundry of lat and lon
    minlat=pointsArr[:,1].min()
    if minlat<-90: minlat=-90.
    maxlat=pointsArr[:,1].max()
    if maxlat>90: maxlat=90.
    minlon=pointsArr[:,0].min()
    maxlon=pointsArr[:,0].max()
    #make them the same for all
    #minlon=20
    #maxlon=220
    #minlat=-80 
    #maxlat=50
    maxlon_wrap=maxlon-360.

    #make image grid mesh
    imgridx, imgridy = numpy.meshgrid(numpy.arange(minlon,maxlon+1),numpy.arange(minlat,maxlat+1))
    # A blank raster image array
    img = numpy.zeros(imgridx.shape, dtype=numpy.uint8)
    # Turn the whole degree values into row/col coords in the array
    row = (pointsArr[:, 1] - minlat)
    col = (pointsArr[:, 0] - minlon) #.clip(0, 359)

    pts=zip(col,row)
    upts=list(set(pts))
    for upt in upts:
        img[upt[1],upt[0]]=pts.count(upt)
    
    #wrap imgridx
    westgrid=numpy.where(imgridx>180.)
    imgridx[westgrid]=imgridx[westgrid]-360.
        
    # For every point, set the corresponding pixel to 1
    #for count in range(len(row)):
    #    img[row[count], col[count]] +=1
    #img=img[minlon:maxlon+1,minlat:maxlat+1]
    
    # plot with matplotlib basemap
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # base map with countries
    m = Basemap(projection='aea',
                llcrnrlon=minlon, 
                llcrnrlat = minlat, 
                urcrnrlon = maxlon_wrap, 
                urcrnrlat = maxlat,
                resolution='l',lon_0=130,lat_0=-15, lat_1=-20, lat_2=-10,lat_ts=-15)
    #m.drawcountries()
    m.drawcoastlines()
    m.drawparallels(numpy.arange(-80.,81.,20.))
    m.drawmeridians(numpy.arange(-180.,181.,20.))

    # obtain map grid
    x,y = m(imgridx,imgridy)
    vmax=img.max()
    if vmax<10: vmax=10
    im = m.pcolormesh(x,y, img, cmap='jet', norm=colors.LogNorm(vmin=1,vmax=vmax))
    cb = m.colorbar(im, "right", size="5%", pad='2%')
    ax.set_title(title.replace('_',' '))
    fig.savefig(outputdir+'/'+title+'.png')
    return outputdir+'/'+title+'.png'

def formSearchDir(collection,product,d):
    """
    return searchdir and pattern based on collection name
    S1 SLC: /g/data/fj7/Copernicus/Sentinel-1/C-SAR/SLC/yyyy/yyyy-mm/*/*_????????T??????_yyyymmddT??????_*.xml
    S2 MSI: /g/data/fj7/Copernicus/Sentinel-2/MSI/L1C/yyyy/yyyy-mm/*/*_yyyymmddT??????_*_????????T??????.xml
    S3 OLCI: /g/data/fj7/Copernicus/Sentinel-3/OLCI/OL_1_EFR___/yyyy/yyyy-mm/yyyy-mm-dd/*_????????T??????_????????T??????_yyyymmddT??????_*.xml
    """
    searchdir=''
    pattern=''
    if collection == 'Sentinel-1':
        searchdir='/g/data/fj7/Copernicus/Sentinel-1/C-SAR/%s/%s'%(product, d.strftime('%Y/%Y-%m/'))
        pattern=datetime.datetime.strftime(d,'*_????????T??????_%Y%m%dT??????_*.xml')
    if collection == 'Sentinel-2':
        searchdir='/g/data/fj7/Copernicus/Sentinel-2/MSI/%s/%s'%(product,d.strftime('%Y/%Y-%m/'))
        pattern=datetime.datetime.strftime(d,'*_%Y%m%dT??????_*_????????T??????.xml')
    if collection == 'Sentinel-3':
        if 'OL' in product: instrument='OLCI'
        if 'SL' in product: instrument='SLSTR'
        if 'SR' in product: instrument='SRAL'
        if len(product)<11: product+='_'*(11-len(product))
        searchdir='/g/data/fj7/Copernicus/Sentinel-3/%s/%s/%s'%(instrument, product, d.strftime('%Y/%Y-%m/%Y-%m-%d/'))
        pattern=datetime.datetime.strftime(d,'*_????????T??????_????????T??????_%Y%m%dT??????_*.xml')
    return searchdir, pattern
                                           

def getXmlList(searchdir,pattern):
    """
    Warlk the searchdir, find filename matching pattern
    """
    xmlfilelist = []
    for (dirname, subdirs, filenames) in os.walk(searchdir):
        for filename in filenames:
            if fnmatch.fnmatch(filename, pattern):
                if '.png' in filename: continue
                fullXmlFilename = os.path.join(dirname, filename)
                xmlfilelist.append(fullXmlFilename)
    return xmlfilelist


if __name__ == "__main__":
    main()

