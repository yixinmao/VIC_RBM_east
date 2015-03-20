#!/usr/bin/env python

# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 1992 - 2012
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2015/1/28
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

# -*- coding: iso-8859-1 -*-
#-------------------------------------------------------------------------------
''' This function file contains many functions for converting rasters or gridded
datasets to polygon shapefiles, and for computing the relationship between polygons
in two shapefiles. Other functions are used for projecting shapefiles and writing
polygon correspondence files to netCDF for use in regridding operations.'''
#-------------------------------------------------------------------------------

'''
 Name:         GDAL_OGR_Functions.py
 Author:       Kevin Sampson
               Associate Scientist
               National Center for Atmospheric Research
               ksampson@ucar.edu

 Created:      10/23/2013
 Modified:     01/28/2015
'''

import os, sys, ogr, osr, gdal, shapely, time, numpy
from gdalconst import *
from shapely import speedups
from shapely.wkb import loads
from shapely.ops import cascaded_union
from netCDF4 import Dataset

# Setup directories to import ogr2ogr.py
pwd = os.path.dirname(__file__)                                                 # Store current working directory
sys.path.append(pwd)                                                            # Append current directory to the python path
#os.chdir(r'E:\Projects\Clark\Grid2Basin')                                      # Set the working directory so that we can import ogr2ogr.py
import ogr2ogr

if speedups.available == True:
    speedups.enable()

def featureIsValid(shgeom):
    '''Return True if this feature passes validity tests.'''
    if not shgeom.is_valid:
        return "Geometry is not valid."
    if shgeom.is_empty:
        return "Geometry is empty."
    return None

def getfieldinfo(field_defn, fieldname):
    '''Get information about field type for buildng the output NetCDF file later'''
    if field_defn.GetType() == ogr.OFTInteger:
        fieldtype = 'integer'                                                   #print "%d" % feat.GetFieldAsInteger(i)
    elif field_defn.GetType() == ogr.OFTReal:
        fieldtype = 'real'                                                      #print "%.3f" % feat.GetFieldAsDouble(i)
        print "field type: OFTReal not currently supported in output NetCDF file."
    elif field_defn.GetType() == ogr.OFTString:
        fieldtype = 'string'                                                    #print "%s" % feat.GetFieldAsString(i)
    else:
        fieldtype = 'string'                                                    #print "%s" % feat.GetFieldAsString(i)
    print "    Field Type for field '%s': %s (%s)" %(fieldname, field_defn.GetType(), fieldtype)
    return fieldtype

def checkfield(layer, fieldname, string1):
    '''Check for existence of provided fieldnames'''
    layerDefinition = layer.GetLayerDefn()
    fieldslist = []
    for i in range(layerDefinition.GetFieldCount()):
        fieldslist.append(layerDefinition.GetFieldDefn(i).GetName())
    if fieldname in fieldslist:
        i = fieldslist.index(fieldname)
        field_defn = layerDefinition.GetFieldDefn(i)
    else:
        print " Field %s not found in input %s. Terminating..." %(fieldname, string1)
        raise SystemExit
    return field_defn, fieldslist

def polygonize_raster(inraster, outputFile, valfield=True):
    '''Function takes a directory and a GDAL-compatible raster and creates a polygon
    shapefile with one rectanglular polygon feature per raster cell. Each feature
    will include information about the i,j index of the grid cell, the latitude and
    longitude of the grid cell centroid (WGS84), and the raster cell data value.'''

    gdal.AllRegister()
    dataset = gdal.Open(inraster, GA_ReadOnly)                                  # Opening the file with GDAL, with read only acces

    # Getting raster dataset information
    print 'Input Raster Size: %s x %s x %s' %(dataset.RasterXSize,dataset.RasterYSize,dataset.RasterCount)
    print 'Projection of input raster: %s' %dataset.GetProjection()
    geotransform = dataset.GetGeoTransform()
    x0 = geotransform[0]                                                        # upper left corner of upper left pixel x
    y0 = geotransform[3]                                                        # upper left corner of upper left pixel y
    pwidth = geotransform[1]                                                    # pixel width, if north is up.
    pheight = geotransform[5]                                                   # pixel height is negative because it's measured from top, if north is up.
    if not geotransform is None:
        print 'Origin (x,y): %s,%s' %(x0,y0)
        print 'Pixel Size (x,y): %s,%s' %(pwidth,pheight)

    # Get raster values as array (must have Numpy 1.8.1 or greater)
    band = dataset.GetRasterBand(1)
    ndv = float(band.GetNoDataValue())
    DataType = band.DataType
    print 'NoData Value: %s' %ndv
    print 'Raster DataType: %s' %gdal.GetDataTypeName(DataType)
    data = band.ReadAsArray(0, 0, dataset.RasterXSize, dataset.RasterYSize).astype(numpy.float)

    # Now convert it to a shapefile with OGR
    driver = ogr.GetDriverByName('Esri Shapefile') # reading from ORG
    datasource = driver.CreateDataSource(outputFile)
    #driver = ogr.GetDriverByName('Memory')                                     # For in-memory vector file
    #datasource = driver.CreateDataSource('out')                                # For in-memory vector file
    if datasource is None:
        print "Creation of output file failed.\n"
        raise SystemExit

    # Create the SpatialReference
    coordinateSystem = osr.SpatialReference()
    coordinateSystem.ImportFromWkt(dataset.GetProjection())                     # Use projection from input raster

    # Create a new layer on the data source with the indicated name, coordinate system, geometry type
    layer = datasource.CreateLayer(outputFile, coordinateSystem, geom_type=ogr.wkbPolygon)
    if layer is None:
        print "Layer creation failed.\n"
        raise SystemExit

    # Create a new field on a layer. Add one attribute
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('i_index', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('j_index', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('lon_cen', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('lat_cen', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('CELLVALUE', ogr.OFTReal))

    # Fetch the schema information for this layer
    LayerDef = layer.GetLayerDefn()

    # Set up transformation in case input is not in WGS84
    point_ref = ogr.osr.SpatialReference()                                      # Get blank spatial reference object
    point_ref.ImportFromEPSG(4326)                                              # WGS84
    ctran=ogr.osr.CoordinateTransformation(coordinateSystem, point_ref)         # Create transformation for converting to WGS84

    # For each cell, create a rectangle
    i = 1
    counter = 0
    for x in range(dataset.RasterXSize):
        j = 1
        for y in range(dataset.RasterYSize):

            # Control flow in order to avoid creating polygons where pixel value is nodata
            if valfield == True and data[y,x] == ndv:
                j += 1
                continue

            # Calculating the polygon's coordinates that frame the raster image
            x00 = x0 + (pwidth * x)
            y00 = y0 - (abs(pheight) * y)
            x1 = x00 + pwidth
            y1 = y00
            x2 = x00
            y2 = y00 - abs(pheight)
            x3 = x1
            y3 = y2

            #create polygon object:
            myRing = ogr.Geometry(type=ogr.wkbLinearRing)
            myRing.AddPoint(x00, y00)
            myRing.AddPoint(x1, y1)
            myRing.AddPoint(x3, y3)
            myRing.AddPoint(x2, y2)
            myRing.AddPoint(x00, y00)                                           #close ring
            geometry = ogr.Geometry(type=ogr.wkbPolygon)
            geometry.AddGeometry(myRing)

            # create point geometry for coordinate system tranformation
            pt = geometry.Centroid()
            pt.Transform(ctran)

            # Create a new feature (attribute and geometry)
            feature = ogr.Feature(LayerDef)
            feature.SetField('id', counter)
            feature.SetField('i_index', i)
            feature.SetField('j_index', j)
            feature.SetField('lon_cen', pt.GetX())
            feature.SetField('lat_cen', pt.GetY())
            feature.SetField('CELLVALUE', data[y,x])                            # Add in the raster value

            # Make a geometry from Shapely object
            feature.SetGeometry(geometry)
            layer.CreateFeature(feature)

            #flush memory
            feature = geometry = myRing = pt = None  # destroy these
            counter += 1
            j += 1
        i += 1

    # Save and close everything
    datasource = layer = None
    return outputFile

def project_to_input(infeatures, outputFile, to_project):
    '''This function projects the input features (to_project) in ESRI Shapefile
    format and uses ogr2ogr.py (if present in the same directory as this script)
    to the coordinate system of 'infeatures' shapefile.'''

    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(infeatures, 0)                                            # Opening the file with GDAL, with read only acces
    lyr = shp.GetLayer()
    spatialref = lyr.GetSpatialRef().ExportToWkt()
    if shp is None:
        print "Open failed.\n"
    ogr2ogr.main(["","-f", "ESRI Shapefile", "-t_srs", spatialref, outputFile, to_project])   #note: main is expecting sys.argv, where the first argument is the script name, so the argument indices in the array need to be offset by 1
    shp = lyr = spatialref = None
    return outputFile

def spatial_weights(inshp1, inshp2, fieldname1, fieldname2, gridflag=0):
    '''This function takes two input shapefiles (both must be in the same coordinate
    reference sysetm) and performs polygon-by-polygon analysis to compute the
    individual weight of each grid cell using OGR routines.  If you wish to use a
    fieldname with non-unique values in inshp1, then it should be the first input
    (inshp1 and fieldname1), and not the second.'''

    tic = time.time()
    spatialweights = {}                                                         # Initiate dictionaries

    # Open the basin shapefile file with OGR, with read only access
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp1 = driver.Open(inshp1, 0)                                               # 0 means read-only. 1 means writeable.
    if shp1 is None:
        print "Open failed.\n"
        raise SystemExit
    layer1 = shp1.GetLayer()
    extent1 = layer1.GetExtent()

    # Open the second shapefile with OGR, with read only access
    shp2 = driver.Open(inshp2, 0)                                               # 0 means read-only. 1 means writeable.
    if shp2 is None:
        print "Open failed.\n"
        raise SystemExit
    layer2 = shp2.GetLayer()

    # Create a Polygon from the extent tuple
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(extent1[0], extent1[2])
    ring.AddPoint(extent1[1], extent1[2])
    ring.AddPoint(extent1[1], extent1[3])
    ring.AddPoint(extent1[0], extent1[3])
    ring.AddPoint(extent1[0], extent1[2])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    layer2.SetSpatialFilter(poly)
    poly = None                                                                 # Eliminate the polygon part

    # Set up transformation in case input is not in WGS84
    coordinateSystem = osr.SpatialReference()
    coordinateSystem.ImportFromWkt(layer2.GetSpatialRef().ExportToWkt())        # Use projection from input shapefile2
    point_ref = osr.SpatialReference()
    geoSR = osr.SpatialReference()
    point_ref.ImportFromEPSG(4326)                                              # WGS84
    ctran = osr.CoordinateTransformation(coordinateSystem, point_ref)           # Create transformation for converting to WGS84

    # Check for existence of provided fieldnames
    field_defn1, fieldslist1 = checkfield(layer1, fieldname1, 'shapefile1')
    field_defn2, fieldslist2 = checkfield(layer2, fieldname2, 'shapefile2')

    # Get information about field types for buildng the output NetCDF file later
    fieldtype1 = getfieldinfo(field_defn1, fieldname1)
    fieldtype2 = getfieldinfo(field_defn2, fieldname2)

    # Load entire polygon shapefile into shapely to create a spatial filter
    if gridflag == 1:
        polys2 = [(feat.GetField(fieldname2), loads(feat.GetGeometryRef().ExportToWkb()), feat.GetField('lon_cen'), feat.GetField('lat_cen'), feat.GetField('i_index'), feat.GetField('j_index')) for feat in layer2]
    elif gridflag == 0:
        polys2 = [(feat.GetField(fieldname2), loads(feat.GetGeometryRef().ExportToWkb())) for feat in layer2]
    layer2.ResetReading()
    if len(polys2) > 1:                                                         # Use cascated union geometry to limit results using spatial filter on layer1
        shgeom = shapely.ops.cascaded_union([poly[1] for poly in polys2])
    else:
        shgeom = polys2[0][1]
    geom = ogr.CreateGeometryFromWkb(shgeom.wkb)                                # Create OGR feature from Shapely geometry
    print "    Cascaded union completed on shapefile2. Setting spatial filter..."

    # Get list of unique field values in input field for shapefile1
    layer1.SetSpatialFilter(geom)
    trials = layer1.GetFeatureCount()
    valuelist = []
    feature = layer1.GetNextFeature()
    while feature is not None:
        valuelist.append(feature.GetField(fieldname1))
        feature.Destroy()
        feature = layer1.GetNextFeature()
    uniques = list(set([x for x in valuelist if valuelist.count(x) == 1]))      # Find all values of fieldname1 that occur exactly once
    nonuniques = list(set([x for x in valuelist if valuelist.count(x) > 1]))    # Find all values of fieldname1 that occur more than once
    del valuelist
    layer1.ResetReading()
    print '    Finished gathering %s unique fieldnames for shapefile1...' %(len(uniques)+len(nonuniques))

    # Initiate polygon geometry list
    polys1 = []

    # Check to see if the number of features is the same as the number of unique field values
    if len(uniques) == trials:                                                  # Fast loop - all unique
        print '    Field %s in shapefile1 contains unique values...' %fieldname1
        feature = layer1.GetNextFeature()
        while feature is not None:
            value = feature.GetField(fieldname1)
            feat = loads(feature.GetGeometryRef().ExportToWkb())
            polys1.append((value, feat))
            feature.Destroy()
            feature = layer1.GetNextFeature()
        layer1.ResetReading()
    else:                                                                       # Slow loop - resolving non-unique values
        print '    Field %s in shapefile1 contains non-unique values...' %fieldname1
        for value in nonuniques:                                                # loop through the non-unique values of input features to create cascaded unions
            layer1.SetAttributeFilter("%s = '%s'" %(fieldname1, value))
            polys = [loads(feat.GetGeometryRef().ExportToWkb()) for feat in layer1]
            layer1.ResetReading()
            polys1.append((value, shapely.ops.cascaded_union(polys)))
            print "        Adding polygon %s to polygons list and performing cascaded union." %value

        # loop through the unique values of input features to create geometry list
        feature = layer1.GetNextFeature()
        while feature is not None:
            value = feature.GetField(fieldname1)
            if value in uniques:
                feat = loads(feature.GetGeometryRef().ExportToWkb())
                polys1.append((value, feat))
                print "        Adding polygon %s to polygons list." %value
            feature.Destroy()
            feature = layer1.GetNextFeature()
        layer1.ResetReading()

    # Loop through the basin polygons
    invalidcount = 0
    for polygon in polys1:
        print "        Gathering correspondence information for polygon %s" %polygon[0]
        shgeom = polygon[1]
        polygon_area = shgeom.area

        try:
            # Attempt to find all overlapping pixels
            if gridflag == 1:                                                   # Read centroid info from attirbute table instead of from geometry
                pixelweights = [(poly[0], (poly[1].intersection(shgeom).area / polygon_area), poly[2], poly[3], poly[4], poly[5]) for poly in polys2 if shgeom.intersects(poly[1])==True]
            elif gridflag == 0:
                pixelweights = [(poly[0], (poly[1].intersection(shgeom).area / polygon_area), ctran.TransformPoint(poly[1].centroid.x, poly[1].centroid.y)[0], ctran.TransformPoint(poly[1].centroid.x, poly[1].centroid.y)[1]) for poly in polys2 if shgeom.intersects(poly[1])==True]
            spatialweights[polygon[0]] = pixelweights
        except:
            invalidcount += 1
            pass

    print "    Invalid Count: %s" %invalidcount
    print "    %s trials completed in %s seconds" %(trials, time.time()-tic)
    return spatialweights, fieldtype1, fieldtype2

def create_ncfile(spatialweights, outputfile, fieldtype1, fieldtype2, gridflag=0):
    '''This function creates an netcdf file with two dimensions ('polyid' and 'gridoverlaps')
    which describe each basin and the grid cells that intersect it.  This leaves
    'whitespace' or blank index values and takes up unecessary disk space.'''

    # Create netcdf file for this simulation
    rootgrp = Dataset(outputfile, 'w', format='NETCDF4')

    # Get length of largest grid to basin overlap
    overlaps = max([len(x) for x in spatialweights.values()])

    # Create dimensions and set other attribute information
    dim1 = 'polyid'
    dim2 = 'gridoverlaps'
    dim = rootgrp.createDimension(dim1, len(spatialweights))
    overlap = rootgrp.createDimension(dim2, overlaps)

    # Handle the data type of the polygon identifier
    if fieldtype1 == 'integer':
        ids = rootgrp.createVariable(dim1,'i4',(dim1))                          # Coordinate Variable (32-bit signed integer)
    elif fieldtype1 == 'string':
        ids = rootgrp.createVariable(dim1,str,(dim1))                           # Coordinate Variable (string type character)

    if fieldtype2 == 'integer':
        ids2 = rootgrp.createVariable('intersector', 'i4', (dim1,dim2))         # (32-bit signed integer)
    elif fieldtype2 == 'string':
        ids2 = rootgrp.createVariable('intersector', str, (dim1,dim2))          # (string type character)

    # Create fixed-length variables
    overlapers = rootgrp.createVariable('overlap','i4',(dim2))                  # Coordinate Variable (32-bit signed integer)
    weights = rootgrp.createVariable('weight', 'f8', (dim1,dim2))               # (64-bit floating point)
    lats = rootgrp.createVariable('latitude', 'f8',(dim1,dim2))                 # (64-bit floating point)
    lons = rootgrp.createVariable('longitude', 'f8',(dim1,dim2))                # (64-bit floating point)
    overlapps = rootgrp.createVariable('overlaps','i4',(dim1))                  # (32-bit signed integer)

    if gridflag == 1:
        iindex = rootgrp.createVariable('i_index', 'i4', (dim1,dim2))           # (32-bit signed integer)
        jindex = rootgrp.createVariable('j_index', 'i4', (dim1,dim2))           # (32-bit signed integer)
        iindex.long_name = 'Index in the x dimension of the raster grid (starting with 0,0 in UL corner)'
        jindex.long_name = 'Index in the y dimension of the raster grid (starting with 0,0 in UL corner)'

    # Set variable descriptions
    weights.long_name = 'fraction of polygon(polyid) intersected by polygon identified by poly2'
    ids.long_name = 'ID of polygon'
    overlapers.long_name = 'Overlap number (1-n)'
    lats.long_name = 'Centroid latitude of intersecting polygon in degrees north in WGS84 EPSG:4326'
    lons.long_name = 'Centroid longitude of intersecting polygon in degrees east in WGS84 EPSG:4326'
    overlapps.long_name = 'Number of intersecting polygons'
    ids2.long_name = 'ID of the polygon that intersetcs polyid'

    # Set units
    lats.units = 'degrees north'
    lons.units = 'degrees east'

    # Fill in global attributes
    rootgrp.history = 'Created %s' %time.ctime()

    # Start filling in elements
    if fieldtype1 == 'integer':
        ids[:] = numpy.array(spatialweights.keys())                         # Ths method works for int-type netcdf variable
    elif fieldtype1 == 'string':
        ids[:] = numpy.array(spatialweights.keys(), dtype=numpy.object)     # This method works for a string-type netcdf variable

    overlapers[:] = numpy.array(range(overlaps))
    overlapps[:] = numpy.array([len(x) for x in spatialweights.values()])
    for i in range(len(spatialweights.keys())):                                 # For each basin
        indexval = ids[i]                                                       # Set id for the basin
        weightslist = [weight for weight in spatialweights[indexval]]           # Generate list of gridcell IDs and spatial weights
        for x in weightslist:                                                   # For each grid cell in each basin
            ids2[i, weightslist.index(x)] = x[0]
            weights[i, weightslist.index(x)] = x[1]
            lons[i, weightslist.index(x)] = x[2]
            lats[i, weightslist.index(x)] = x[3]
            if gridflag == 1:
                iindex[i, weightslist.index(x)] = x[4]
                jindex[i, weightslist.index(x)] = x[5]

    # Close file
    rootgrp.close()

def create_ragged_ncfile(spatialweights, outputfile, fieldtype1, fieldtype2, gridflag=0):
    '''This function creates a ragged-array netCDF file, with variable length dimensions
    for each basin. '''

    # Create netcdf file for this simulation
    rootgrp = Dataset(outputfile, 'w', format='NETCDF4')

    # Create variable length (ragged) data type
    vlen_i = rootgrp.createVLType(numpy.int32, 'vlen_int')
    vlen_f = rootgrp.createVLType(numpy.float64, 'vlen_float')

    # Create dimensions and set other attribute information
    dim1 = 'polyid'
    dim = rootgrp.createDimension(dim1, len(spatialweights))

    # Create variable-length variables
    weights = rootgrp.createVariable('weight', vlen_f, (dim1))                  # (64-bit floating point)
    lats = rootgrp.createVariable('latitude', vlen_f, (dim1))                   # (64-bit floating point)
    lons = rootgrp.createVariable('longitude', vlen_f, (dim1))                  # (64-bit floating point)
    if fieldtype2 == 'integer':
        ids2 = rootgrp.createVariable('intersector', vlen_i, (dim1))            # Coordinate Variable (32-bit signed integer)
    elif fieldtype2 == 'string':
        ids2 = rootgrp.createVariable('intersector', vlen_f, (dim1))            # Coordinate Variable (string type character)

    if gridflag == 1:
        iindex = rootgrp.createVariable('i_index', vlen_i, (dim1))              # (32-bit signed integer)
        jindex = rootgrp.createVariable('j_index', vlen_i, (dim1))              # (32-bit signed integer)
        iindex.long_name = 'Index in the x dimension of the raster grid (starting with 0,0 in UL corner)'
        jindex.long_name = 'Index in the y dimension of the raster grid (starting with 0,0 in UL corner)'

    # Handle the data type of the polygon identifier
    if fieldtype1 == 'integer':
        ids = rootgrp.createVariable(dim1,'i4',(dim1))                          # Coordinate Variable (32-bit signed integer)
    elif fieldtype1 == 'string':
        ids = rootgrp.createVariable(dim1,str,(dim1))                           # Coordinate Variable (string type character)

    # Create fixed-length variables
    overlaps = rootgrp.createVariable('overlaps','i4',(dim1))                   # 32-bit signed integer

    # Set variable descriptions
    weights.long_name = 'fraction of polygon(polyid) intersected by polygon identified by poly2'
    ids.long_name = 'ID of polygon'
    overlaps.long_name = 'Number of intersecting polygons'
    lats.long_name = 'Centroid latitude of intersecting polygon in degrees north in WGS84 EPSG:4326'
    lons.long_name = 'Centroid longitude of intersecting polygon in degrees east in WGS84 EPSG:4326'
    ids2.long_name = 'ID of the polygon that intersetcs polyid'

    # Set units
    lats.units = 'degrees north'
    lons.units = 'degrees east'

    # Fill in global attributes
    rootgrp.history = 'Created %s' %time.ctime()

    # Start filling in elements
    if fieldtype1 == 'integer':
        ids[:] = numpy.array(spatialweights.keys())                         # Ths method works for int-type netcdf variable
    elif fieldtype1 == 'string':
        ids[:] = numpy.array(spatialweights.keys(), dtype=numpy.object)     # This method works for a string-type netcdf variable

    #ids[:] = numpy.array([str(x) for x in spatialweights.keys()], dtype=numpy.object)   # This method will convert to string-type netcdf variable from int
    overlaps[:] = numpy.array([len(x) for x in spatialweights.values()])
    ids2[:] = numpy.array([numpy.array([x[0] for x in spatialweights[spatialweights.keys()[i]]]) for i in range(len(spatialweights.keys()))], dtype=numpy.object)
    weights[:] = numpy.array([numpy.array([weight[1] for weight in spatialweights[spatialweights.keys()[i]]]) for i in range(len(spatialweights.keys()))], dtype=numpy.object)
    lons[:] = numpy.array([numpy.array([x[2] for x in spatialweights[spatialweights.keys()[i]]]) for i in range(len(spatialweights.keys()))], dtype=numpy.object)
    lats[:] = numpy.array([numpy.array([x[3] for x in spatialweights[spatialweights.keys()[i]]]) for i in range(len(spatialweights.keys()))], dtype=numpy.object)

    if gridflag == 1:
        iindex[:] = numpy.array([numpy.array([x[4] for x in spatialweights[spatialweights.keys()[i]]]) for i in range(len(spatialweights.keys()))], dtype=numpy.object)
        jindex[:] = numpy.array([numpy.array([x[5] for x in spatialweights[spatialweights.keys()[i]]]) for i in range(len(spatialweights.keys()))], dtype=numpy.object)

    # Close file
    rootgrp.close()

if __name__ == '__main__':
    print "Cannot run as __main__ !"
    raise SystemExit