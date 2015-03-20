#!/usr/bin/env python

# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 1992 - 2012
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2015/1/28
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

'''
 Name:         poly2poly.py
 Author:       Kevin Sampson
               Associate Scientist
               National Center for Atmospheric Research
               ksampson@ucar.edu

 Created:      10/23/2013
 Modified:     01/28/2015
'''

# Import modules
import os, sys, time, getopt, glob
tic = time.time()
sys.dont_write_bytecode = True
import GDAL_OGR_Functions

# Setup directories to import ogr2ogr.py
pwd = os.path.dirname(__file__)                                                 # Store current working directory
sys.path.append(pwd)                                                            # Append current directory to the python path
#os.chdir(pwd)                                                                  # Set the working directory so that we can import ogr2ogr.py
import ogr2ogr                                                                  # Import GDAL OGR2OGR script for use in this function

# Screen print in case invalid parameters are given
use = '''
Usage: %s <Polygon1> <FieldName1> <Polygon2> <FieldName2> [<gridflag>] [<outputNC>] [<VLENflag>]
      <Polygon1>     -> Polygon Shapefile (full path)
      <FieldName1>   -> Polygon identifier field name for Polygon1
      <Polygon2>     -> Polygon Shapefile2 (full path)
      <FieldName2>   -> Polygon identifier field name for Polygon 2
      <gridflag>     -> Optional - Type 'GRID' to indicate Polygon2 was generated from grid2shp.py
      <outputNC>     -> Optional - Full path to output netCDF file (.nc extentsion required)
      <VLENflag>     -> Optional - Type 'VLEN' to indicate output NC file should be variable length array type
'''

def gridtobasin_function(infeatures1, fieldname1, infeatures2, fieldname2, gridflag=0, outputnc=None, raggednc=False):

    # Input directories
    indir = os.path.dirname(infeatures2)                                        # Use directory of output file
    print 'Input Directory: %s' %indir

    # Output Files
    if outputnc is None:
        outncFile = os.path.join(indir, os.path.basename(infeatures1)[:-4]+'.nc')
    else:
        outncFile = outputnc
    projfeatures = os.path.join(indir, os.path.basename(infeatures2).replace(".shp", "_proj.shp"))

    # Check if files exist and delete existing files if necessary
    if os.path.isfile(projfeatures) == True:
        print "Removing existing file: %s" %projfeatures
        remover = [os.remove(infile) for infile in glob.glob(projfeatures.replace('.shp', '.*'))]
    if os.path.isfile(outncFile) == True:
        print "Removing existing file: %s" %outncFile
        os.remove(outncFile)

    # Step 1: Project gridded polygons to coordinate system of input features
    tic1 = time.time()
    print "Step 1: Projecting <Polygon2> to the coordinate system of <Polygon1> starting."
    projfeatures = GDAL_OGR_Functions.project_to_input(infeatures1, projfeatures, infeatures2)
    print "Step 1: Projecting <Polygon2> to the coordinate system of <Polygon1> completed in %s seconds." %(time.time()-tic1)

    # Step 2: Calculate per gridcell feature coverage fractions
    tic1 = time.time()
    print "Step 2: Compute spatial weights starting."
    spatialweights, fieldtype1, fieldtype2 = GDAL_OGR_Functions.spatial_weights(infeatures1, projfeatures, fieldname1, fieldname2, gridflag)
    print "Step 2: Compute spatial weights completed in %s seconds." %(time.time()-tic1)

    # Step 3: Write weights to ragged array netcdf file
    tic1 = time.time()
    print "Step 3: Create correspondence netCDF file starting."
    if raggednc == True:
        GDAL_OGR_Functions.create_ragged_ncfile(spatialweights, outncFile, fieldtype1, fieldtype2, gridflag)
    else:
        GDAL_OGR_Functions.create_ncfile(spatialweights, outncFile, fieldtype1, fieldtype2, gridflag)
    print "Step 3: Create correspondence netCDF file completed in %s seconds." %(time.time()-tic1)

    # Remove intermediate projected polygons
    if os.path.isfile(projfeatures) == True:
        print "Removing intermediate projected shapefile: %s" %projfeatures
        remover = [os.remove(infile) for infile in glob.glob(projfeatures.replace('.shp', '.*'))]

    print "Total time elapsed: %s seconds." %(time.time()-tic)

if __name__ == '__main__':

    print "Starting __main__ function"

    def usage():
      sys.stderr.write(use % sys.argv[0])
      sys.exit(1)
    try:
      (opts, args) = getopt.getopt(sys.argv[1:], 'h')
    except getopt.error:
      usage()

    if len(sys.argv) < 5 or len(sys.argv) > 8:
      usage()
    else:

        # Input arguments
        infeatures1 = sys.argv[1]
        fieldname1 = sys.argv[2]
        infeatures2 = sys.argv[3]
        fieldname2 = sys.argv[4]

        # Default arguments
        gridf = 0
        outnc = None
        vlenflag = False

        # Conditionals for optional input arguments
        if len(sys.argv) > 5:
            if 'GRID' in sys.argv[:]:
                gridf = 1

            if 'VLEN' in sys.argv[:]:
                vlenflag = True

            if os.path.exists(os.path.dirname(sys.argv[5])) == True:
                outnc = sys.argv[5]
            if len(sys.argv) > 6:
                if os.path.exists(os.path.dirname(sys.argv[6])) == True:
                    outnc = sys.argv[6]
            if len(sys.argv) > 7:
                if os.path.exists(os.path.dirname(sys.argv[7])) == True:
                    outnc = sys.argv[7]

        gridtobasin_function(infeatures1, fieldname1, infeatures2, fieldname2, gridflag=gridf, outputnc=outnc, raggednc=vlenflag)