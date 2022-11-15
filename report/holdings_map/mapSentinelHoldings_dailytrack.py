#!/usr/bin/env python
"""
Generate heatmap for Sentinel-1/2/3 acquisitions.

Adapted from Neil Flood's mapSentinelHoldings.py and Damien Ayers's datacube map.

"""
from __future__ import print_function, division

import sys
import os
import argparse
import fnmatch
import pickle

import datetime
import numpy as np

from functools import partial
import pyproj
import shapely.ops
from shapely.wkt import loads
from xml.etree import ElementTree

from affine import Affine
from rasterio import features

from mpl_toolkits.basemap import Basemap 
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import cartopy.crs as ccrs


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
    Walk the searchdir, find filename matching pattern
    """
    xmlfilelist = []
    for (dirname, subdirs, filenames) in os.walk(searchdir):
        for filename in filenames:
            if fnmatch.fnmatch(filename, pattern):
                if '.png' in filename: continue
                fullXmlFilename = os.path.join(dirname, filename)
                xmlfilelist.append(fullXmlFilename)
    return xmlfilelist


def polygon_from_xml(metafile, valid_mode='IW', valid_pols=['VH','VV']):
    """
    load wkt footprint from a xml file
    """
    root = ElementTree.parse(metafile)
    footprint = root.findall('ESA_TILEOUTLINE_FOOTPRINT_WKT')[0].text
    if 'MULTIPOLYGON' in footprint: return None # no Multipolygon
    mode = root.findall('MODE')[0].get('value')
    if mode!=valid_mode: return None # only IW mode
    pols = root.findall('POLARISATION')[0].get('values').split(',')
    for pol in valid_pols:
        if not pol in pols: return None #required pols
    if len(pols)!=len(valid_pols): return None
    try:
        poly=loads(footprint.strip())
    except:
        print("Problem with footprint:",footprint.strip(), metafile)
        return None
    if not poly.is_valid: poly=poly.buffer(0)
    if not poly.is_valid:
        print("not valid?:",footprint.strip(), metafile)
        return None
    #cross dateline?
    bounds=poly.bounds
    if abs(bounds[2]-bounds[0])>180:
        line=loads("LINESTRING(-180 90, -180 -90)")
        newpoly=shapely.ops.split(poly,line)
        for p in newpoly:
            b=p.bounds
            if abs(bounds[2]-bounds[0])>180:
                return None
        return [p for p in newpoly]
    return [poly]


def shape_for_list(xmlList, tolerance=10,transform=None,**kwargs):
    """
    Combine shapes from list of xml files
    """
    ps=[]
    for mf in xmlList:
        p=polygon_from_xml(mf,**kwargs)
        if p: ps.extend(p)
    if len(ps)>0:
        polys = shapely.ops.unary_union(ps)
        polys = polys.simplify(tolerance=tolerance)
        if transform:
            project = partial(
                pyproj.transform,
                pyproj.Proj(init='epsg:4326'),
                #pyproj.Proj(init='epsg:3577')
                pyproj.Proj(init=transform)
                )
            polys=shapely.ops.transform(project, polys)
        return polys
    else:
        return None

def polys_for_days(collection,product,start,end,transform=False):
    polys=[]
    d=start
    while d<=end:
        searchdir, pattern = formSearchDir(collection,product,d)
        if len(searchdir)>0:
            xmlList=getXmlList(searchdir, pattern)
            if len(xmlList)>0:
                track=shape_for_list(xmlList,transform=transform)
                if track: polys.append(track)
        d+=datetime.timedelta(1)
    return polys


def affine_from_extent(extent, width=1024):
    x_min, x_max, y_min, y_max = extent
    x_res = 1.*(x_max - x_min) / width
    height = int((y_max - y_min) / x_res)
    transform = Affine.translation(x_min, y_min) * Affine.scale(x_res, x_res)
    return transform, height


def rasterize_polys_to_heatmap(polys, transform, height, width):
    output = np.zeros((height, width), dtype='int16')
    tmp_image = np.zeros((height, width), dtype='int16')
    for poly in polys:
        if not hasattr(poly,'__iter__'):poly=[poly]
        features.rasterize(poly, (height, width), out=tmp_image, transform=transform, dtype='int16')
        output += tmp_image
    return output


def heatmap_for_days(collection,product,start,end, transform, height, width, projtransform=False,
                     **kwargs):
    output = np.zeros((height, width), dtype='int16')
    d=start
    while d<=end:
        searchdir, pattern = formSearchDir(collection,product,d)
        if len(searchdir)>0:
            xmlList=getXmlList(searchdir, pattern)
            if len(xmlList)>0:
                track=shape_for_list(xmlList,transform=projtransform, **kwargs)
                if track:
                    if type(track) is not list: track=[track]
                    tmp_image=rasterize_polys_to_heatmap(track,transform,height,width)
                    output += 1*(tmp_image>0)                    
        d+=datetime.timedelta(1)
    return output


def create_heatmap_plot(count_data, projection, extent, title):
    # Create figure and axes
    fig = plt.figure(figsize=(14, 10))
    ax = plt.axes(projection=projection)

    # Setup some output
    plt.title(title,fontsize=32,verticalalignment='bottom')
    ax.background_patch.set_visible(False)
    ax.outline_patch.set_visible(False)
    ax.coastlines()
    #ax.gridlines()

    # Display our raster count
    cmap=plt.cm.jet
    cmap.set_under(color='white')
    img_ax = ax.imshow(count_data, cmap=cmap, vmin=1, extent=extent, transform=projection)

    # Add a colorbar
    cbar = plt.colorbar(img_ax, ax=ax)    
    fig.savefig('%s.png'%('_'.join(title.split())))


def main(collection, product, start, end, ausonly=True, heatmappickle=None, title=None, **kwargs):

    if ausonly:
        #Define Australian Extents and a suitable Affine Translation for the whole contintent   
        proj='epsg:3577'
        projection = ccrs.epsg(3577)
        AOI_EXTENT=-2.69e6, 2.58e6, -5.1e6, -1.e6
    else:
        proj='esri:102030'
        projection = ccrs.LambertConformal(central_longitude=125, central_latitude=-15,standard_parallels=(7,-32))
        AOI_EXTENT = -1e7, 1e7, -7.5e6, 6e6 
        
    width = 1024
    transform, height = affine_from_extent(AOI_EXTENT,width=width)
    
    heatmap=None
    if heatmappickle:
        if os.path.exists(heatmappickle):
            heatmap=pickle.load(open(heatmappickle,'rb'))
    if heatmap is None:
        heatmap = heatmap_for_days(collection,product,start,end, transform, height, width, 
                                   projtransform=proj,**kwargs)
            
    if heatmappickle:
        if not os.path.exists(heatmappickle):
            pickle.dump(heatmap,open(heatmappickle,'wb'))
    
    if title is None:
        title='{} {} {} to {} {}'.format(collection, product, 
                                         start.strftime('%Y-%m-%d'),end.strftime('%Y-%m-%d'),
                                         proj.split(':')[-1])
    create_heatmap_plot(heatmap, projection, AOI_EXTENT, title)

if __name__ == "__main__":
    #find data and extract polygons
    collection='Sentinel-1'
    product='SLC'
    mode='IW'
    for year in [2015,2016,2017]:
        for pols in [['VH','VV']]:#,['HV','HH'],['VV'],['HH'],['HV'],['VH']]:
            title='%s %s %s %s %d'%(collection,product,mode,'+'.join(pols), year)
            heatmappickle='_'.join(title.split(' '))+'.pkl'
            start=datetime.datetime(year,1,1)
            end=datetime.datetime(year,12,31)
            main(collection, product, start, end, ausonly=True, valid_mode=mode,valid_pols=pols, 
                 heatmappickle=heatmappickle,title=title)
