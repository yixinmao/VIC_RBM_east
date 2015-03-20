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
 Name:         grid2shp.py
 Author:       Kevin Sampson
               Associate Scientist
               National Center for Atmospheric Research
               ksampson@ucar.edu

 Created:      10/23/2013
 Modified:     01/28/2015
'''

# Import modules
import time
tic = time.time()
import os, sys, glob, getopt
sys.dont_write_bytecode = True
import GDAL_OGR_Functions

# Screen print in case invalid parameters are given
use = '''
Usage: %s <Raster> <outputSHP> [<Include>]
      <Raster>     -> Raster to be converted to shapefile (full path)
      <outputSHP>  -> Output shapefile (full path)
      <Include>    -> Optional - Type 'INCLUDE' to include cells with NoData values
'''

def grid2shp(inraster, outputSHP, nodata=True):

    # Input directories
    #indir = os.path.dirname(__file__)                                           # Store current working directory
    indir = os.path.dirname(outputSHP)                                          # Use directory of output file
    print 'Input Directory: %s' %indir

    # Output Files
    outputFile = os.path.basename(outputSHP)
    print 'Output Shapefile: %s' %outputFile

    # Check if files exist and delete
    remover = [os.remove(infile) for infile in glob.glob(outputSHP.replace('.shp', '.*')) if os.path.isfile(infile)==True and infile[-4:] in ['.shp', '.shx', '.xml', '.sbx', '.sbn', '.prj', '.dbf']]
    if len(remover) >= 1:
        print "Removed existing file: %s" %outputFile

    # Step 1: Convert gridded dataset to polygon features
    gridpolys = GDAL_OGR_Functions.polygonize_raster(inraster, outputSHP, valfield=nodata)
    print "Created files: %s" %(gridpolys)
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

    print "Length of argv: %s" %len(sys.argv)
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        usage()
    else:
        # Input arguments
        inraster = sys.argv[1]
        outputSHP = sys.argv[2]
        if len(sys.argv) == 4:
            if sys.argv[3] == 'INCLUDE':
                nodata = False
            else:
                nodata = True
        else:
            nodata = True
        grid2shp(inraster, outputSHP, nodata)
