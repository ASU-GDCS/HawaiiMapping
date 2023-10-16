#!/usr/bin/env python
import sys, os
import argparse
import numpy

import rasterio # Use for raster reads

# The following will be used for masking based on features to extract covered pixels
import shapely.geometry as geom
import rasterio.features as feat

from osgeo import ogr # Use this for vector input, since it can handle sql queries
from osgeo import osr # Use this for transforming coordinates from vector to raster CRS

# Keep OGR errors from happening siliently
ogr.UseExceptions()


import tqdm # for reporting progress


## Function to make terminal output occur immediately (output caching can cause problems)
#########################################################################################
def write_and_flush(*args,**kwargs):
    sys.stdout.write(*args,**kwargs)
    sys.stdout.flush()

## Function to process raster specification arguments into fieldname, filename, bandlist
def breakrasterarg(rasterarg):
    ##Split at equals sign and get out field name
    #############################################
    # Check if there is an equal sign and make up a fieldname if not
    if rasterarg.find("=") < 0:
        return "", rasterarg, []
    try:
        fieldname, rightpart = rasterarg.split("=")
    except ValueError as exc:
        raise RuntimeError("Wrong number of '=' found in argstring {}, should be 0 or 1".format(rasterarg))

    ##Check for a : meaning band numbers are specified
    ##################################################
    if rightpart.find(":") < 0:
        # If we get here, then there is no colon and bands are not specified
        return fieldname, rightpart, []

    # If there is a colon, then try to split into proper parts - fails if there is more than 1 :
    try:
        filename, bandspec = rightpart.split(":")
    except ValueError as exc:
        raise RuntimeError("Wrong number of ':' found in argstring {}, should be 0 or 1".format(rasterarg))

    ##For each comma-separated range in the bandspec, try to collect appropriate bands
    ##################################################################################
    bandlist = []
    for rangespec in bandspec.split(","):
        # Ignoring empty strings from double commas
        if len(rangespec) > 0:
            # Split on dash
            rangeends = rangespec.split("-")
            if len(rangeends) == 1:
                # If no dash, then there should just be an integer
                bandlist.append(int(rangeends[0]))
            elif len(rangeends) == 2:
                # If one dash, then there should be two integers
                bandlist.extend(list(range(int(rangeends[0]),int(rangeends[1])+1)))
            else:
                # Something is wrong
                raise RuntimeError("Wrong number of '-' found in band range {}, should be 0 or 1".format(rangespec))
    return fieldname, filename, bandlist


## From the provided list of positive integers figure out how many digits
##  are needed to represent the biggest, i.e. 1-9=1, 10-99=2, 100-999=3
##  Will be used to format field names with band numbers
def max_decimal_digits(intlist):
    return int(numpy.ceil(numpy.log10(numpy.array(intlist)+0.01).max()))


## See if we need to shorten base to meet the field name size limit of
## certain formats (cough, cough, ESRI). If so, try adding sequential
## digits to truncated name until we don't clash with existing fields
def modify_basename(base, maxlen=32, banddigits=4, currentbases=[]):
    # First check original name to see if it is too long
    if (len(base) + banddigits) > maxlen:
        # Compute how long it can be ad check if > 0 
        redlen = maxlen - banddigits
        if redlen < 1:
            raise RuntimeError("Cant reduce fieldname {} to {} characters".format(
                base,
                redlen))
        newbase = base[:redlen]
    else:
        newbase = base

    # Now check that name does not already exist, and return if so
    if newbase not in currentbases:
        return newbase

    # If we get here, then the reduced name already exists, and we have to shorten more
    # so we can add single digits to the end to make it unique
    possible_endings = ["{}".format(i) for i in range(1,10)]

    # Compute a new, new length
    redlen = maxlen - banddigits - 1
    if redlen < 1:
        raise RuntimeError("Cant reduce fieldname {} to {} characters".format(
            base,
            redlen))

    # Go through endings until we find a name that doesn't exist (hopefully)
    for tmpend in possible_endings:
        newbase = "{}{}".format(base[:redlen],tmpend)
        # Now check that name does not already exist, and return if so
        if newbase not in currentbases:
            return newbase

    # If we got here, then we exhausted the endings - just fail. Hopefully user doesn't specify more than 9
    # field names with similar beginnings
    raise RuntimeError("Ran out of suffixes trying to shorten field name {} to {} characters".format(
        base,
        redlen))


## With a given basename and a list of band numbers, create the field names
def create_band_names(base, bands, banddigits=4):
    if len(bands) < 2:
        return [base]
    else:
        # Now add the band numbers to the (maybe shortened) basename
        bnformat = "{{}}{{:0{}d}}".format(banddigits)
        return [bnformat.format(base,b) for b in bands]

def indices_to_window(rows, cols):
    offsetc = min(cols)
    offsetr = min(rows)
    maxc = max(cols)
    maxr = max(rows)
    return rasterio.windows.Window(
               offsetc, offsetr, (maxc - offsetc) + 1, (maxr - offsetr) + 1)

def extent_to_window(ext, rast):
    llr, llc = rast.index(*ext[:2])
    urr, urc = rast.index(*ext[2:])
    return indices_to_window([llr,urr],[llc,urc])

def main():

    ################################
    ## Set up and read the arguments
    ################################
    parser = argparse.ArgumentParser(description="Extract raster data from within vector points or regions")
    parser.add_argument("--nodata","-nodata",type=float,default=None,
                        help="Specifiy no data value in raster (if not identified in metadata)")
    parser.add_argument("--sql","-sql",type=str,default=None,
                        help="Apply sql query to vector layer before collecting stats")
    parser.add_argument("--ignore","-ignore","-i",action='store_true',
                        help="Ignore raster no data value - include pixels marked as invalid")
    parser.add_argument("--alltouched","-alltouched","-a",action='store_true',
                        help="Include all touched pixels in areal vectors (default is just those with center inside)")
    parser.add_argument("--nln","-nln",type=str,default=None,
                        help="Name out output layer - for formats that support this")
    parser.add_argument("--format","-format","-of",type=str,default=None,
                        help="OGR format of output file - default is same as input vector file")
    parser.add_argument("input",help="Input vector file")
    parser.add_argument("output",help="Output vector file")
    parser.add_argument("rasters",nargs='+',
                        help="Entries of the format FIELDNAME=FILENAME[:BANDLIST], where "+\
                             "optional BANDLIST is 1-indexed, and can be a list separated "+\
                             "by commas and can include ranges (i.e, 1:55). All rasters "+\
                             "should bin in same spatial reference. If vector and at least "+\
                             "one raster have crs in metainfo, then coordinate transformation "+\
                             "of vector features will occur")
    args = parser.parse_args()
    # args = parser.parse_args(["SKI_locations_2020_20220119_transoutlines.gpkg","test.gpkg","Depth=../Harish_RF_Biomass_Analysis/data_aoi_extent/output_2m_maps/SKI_blended_depth_0m_to_16m.tif","Northing=../Ethans_work/south_kona_stratsamp/redo_2m/data/aoi_extent/northing.tif","Sand=../Harish_RF_Biomass_Analysis/data_aoi_extent/output_2m_maps/SKI_trunc_sand_prop.tif","Algal=../Harish_RF_Biomass_Analysis/data_aoi_extent/output_2m_maps/SKI_trunc_algal_prop.tif","Live=../Harish_RF_Biomass_Analysis/data_aoi_extent/output_2m_maps/SKI_trunc_live_prop.tif","Planar=../Harish_RF_Biomass_Analysis/data_aoi_extent/output_2m_maps/SKI_trunc_xtrafine_rugosity.tif","Curv=../Harish_RF_Biomass_Analysis/data_aoi_extent/output_2m_maps/SKI_trunc_absprofcurv_0_16.tif"]) 
    ###########################
    ## Try to open input vector
    ###########################
    try:
        indata = ogr.Open( args.input, True )
        if indata is None:
            raise RuntimeError("OGR failed to open input file")
    except Exception as exc:
        write_and_flush("OGR could not not open input file {} - {}.\n".format(args.input, exc))
        sys.exit(1)

    ##############################################
    ## Figure out which format output should be in
    ##############################################
    if args.format is None:
        #Take driver from input
        outdriver = indata.GetDriver()
        args.format = outdriver.GetName()
        write_and_flush("Taking output format from input: {}\n".format(args.format))
    else:
        outdriver = ogr.GetDriverByName(args.format)
        if outdriver is None:
            write_and_flush("Could not open driver of format {}\n".format(args.format))
            sys.exit(99)

    ###################################################################
    ##Try to run SQL statement if given, otherwise just get first layer
    ###################################################################
    if args.sql is not None:
        write_and_flush("Executing SQL statement: {}\n".format(args.sql))
        inlayer = indata.ExecuteSQL(args.sql)
    else:
        inlayer = indata.GetLayer( 0 ) 
    if inlayer is None:
        write_and_flush("Could not get layer/result from {}\n".format(args.input))
        sys.exit(2)

    #########################################
    ## Collect input layer definition and CRS
    #########################################
    indefn = inlayer.GetLayerDefn()
    infieldcount = indefn.GetFieldCount() # How many fields are in the vector input
    incrs = inlayer.GetSpatialRef()

    #######################################################################
    ## Recover raster filenames from arguments and gather info from rasters
    ## Each should be formatted like FIELDNAME=FILENAME[:BAND,BAND-BAND,...]   
    #######################################################################
    
    ##Because the same raster could be specified in multiple entries, rasters and
    ## associated into will be kept in a dictionary indexed by filename
    #############################################################################
    raster_dict = dict()

    # These just hold the rasternames, fieldnames, and bandlists in order
    ordered_rastnames = list()
    ordered_bandlists = list()

    ## Data for each specified field=raster pair will be collected on a pixel-by-
    ## pixel basis, so we'll need to keep track of which raster band goes to each
    ## created field. Here, we start to prepare for this indexing
    #############################################################################
    newfields = [("PixelNum",ogr.OFTInteger),
                 ("FeatPixel",ogr.OFTInteger),
                 ("Rast1Row",ogr.OFTInteger),
                 ("Rast1Col",ogr.OFTInteger),
                 ("X",ogr.OFTReal),
                 ("Y",ogr.OFTReal)]
    metafieldcount = len(newfields)
    current_findex_offset = metafieldcount
    
    used_basenames = list() # This will be used to make sure output field names stay unique

    ##############################################################
    ## Collect information from rasters specified in the arguments
    ##############################################################
    rastcrs = None
    for ri, a in enumerate(args.rasters):
        # Process the raster arguments into fieldname, filename, bandlist
        arg_field, arg_name, arg_bands = breakrasterarg(a)
        # If field name was not specified, just make a generic one
        if arg_field is None:
            arg_field = "image{:02d}".format(ri+1)

        # Add to file list
        ordered_rastnames.append(arg_name)
        # Add band list
        ordered_bandlists.append(arg_bands)

        ##Get raster info and add to dictionary if not there already
        ############################################################
        
        # Open the raster and save metadata
        try:
            tmpref = rasterio.open(arg_name)
        except rasterio.errors.RasterioIOError as exc:
            raise RuntimeError("Could not open raster {} - {}".format(arg_name,str(exc)))

        tmpmeta = tmpref.meta

        # Check CRS and create OSR SpatialReference if it exists
        tmpcrs = None
        if tmpref.crs:
            #We have a crs in the metadata, convert it to osr
            tmpcrs = osr.SpatialReference()
            tmpcrs.ImportFromProj4(tmpref.crs.to_proj4())
            
            # Make this the crs that applies to all rasters if not yet discovered
            if not rastcrs:
                rastcrs = tmpcrs

        # Fill arg_bands if band list was not specified
        ###############################################
        if not arg_bands:
            arg_bands = list(range(1,tmpmeta['count']+1))

        # Create the new band names
        ###########################
        # How many digits do we need for band numbers
        tmpdigits = max_decimal_digits(arg_bands)
        if len(arg_bands) < 2:
            tmpdigits=0

        # might need to truncate filenames
        if args.format == "ESRI Shapefile":
            maxlen = 10
        else:
            maxlen = 32

        # Call make_fieldnames
        red_field_name = modify_basename(arg_field, 
                                         maxlen=maxlen,
                                         banddigits=tmpdigits,
                                         currentbases=used_basenames)

        # Warn if field name is changing
        if red_field_name != arg_field:
            warnings.warn(
                    "Field name {} changed to {} to avoid collision or meet length restriction".format(
                        arg_field,
                        red_field_name))

        # Add the modified field name to our list
        used_basenames.append(red_field_name)

        # Okay collect the fieldnames for this entry and get their indices in output vector
        tmpbandnames = create_band_names(red_field_name, arg_bands, banddigits = tmpdigits)
        tmpindices = [current_findex_offset + i for i in range(len(arg_bands))]

        # Add the fields to the new fields list
        newfields.extend(zip(tmpbandnames, [ogr.OFTReal]*len(arg_bands)))
        
        # Increment the field index offset
        current_findex_offset = len(newfields)

        # Put all this info in the dictionary
        #####################################
        if arg_name in raster_dict:
            raster_dict[arg_name]['Fields'].extend(tmpbandnames)
            raster_dict[arg_name]['Indices'].extend(tmpindices)
            raster_dict[arg_name]['Bands'].extend(arg_bands)
        else:
            raster_dict[arg_name] = { 'Reference':tmpref, 
                                      'CRS':tmpcrs,
                                      'Meta':tmpmeta,
                                      'Bands':arg_bands,
                                      'Fields':tmpbandnames,
                                      'Indices':tmpindices}

    ##############################################################
    ## Create functions to transform from vector CRS to raster CRS
    ##############################################################
    gt = None
    igt = None
    if (incrs is not None) and (rastcrs is not None):
        gt = osr.CoordinateTransformation(incrs,rastcrs)
        igt = osr.CoordinateTransformation(rastcrs,incrs)

    ################################################################################
    ## Try to apply a spatial filter on the input layer based on first raster extent
    ################################################################################
    firstrast = raster_dict[ordered_rastnames[0]]['Reference']
    firstrastmeta = raster_dict[ordered_rastnames[0]]['Meta']
    firstrastenv = list(firstrast.bounds)
    if igt:
        # Get extent of the first raster in input layer CRS coordinates
        transwind = list(igt.TransformPoint(*firstrastenv[:2])[:2])
        transwind.extend(list(igt.TransformPoint(*firstrastenv[2:])[:2]))

        # Apply this as a SpatialFilter to the input layer
        write_and_flush("Clipping input layer to extent of first raster: {}\n".format(transwind))
        inlayer.SetSpatialFilterRect(*transwind)


    ##################################
    ## Create output dataset and layer
    ##################################

    # How many features pass the filter
    infeatcount = inlayer.GetFeatureCount()
    # If none, then exit
    if infeatcount < 1:
        write_and_flush("No unfiltered input features found\n")
        sys.exit(3)

    # Try to open output vector
    outref = outdriver.CreateDataSource(args.output)
    if outref is None:
        write_and_flush("Could not create vector file of format {}\n".format(args.format))
        sys.exit(2)

    # Create a new layer
    if args.nln is None:
        args.nln = os.path.splitext(os.path.basename(args.output))[0]
    outlayer = outref.CreateLayer(args.nln, incrs, ogr.wkbPoint)
    if outlayer is None:
        write_and_flush("Could not create layer {} of type {}\n".format(args.nln,ogr.GeometryTypeToName(inlayer.GetGeomType())))
        sys.exit(2)

    ################################
    ## Create fields of output layer
    ################################

    # Start a batch of changes
    outlayer.StartTransaction()

    # Copy over field definitions from input file
    outfields = []
    write_and_flush("Copying fields from input layer:\n")
    for i in range(indefn.GetFieldCount()):
        fieldDef = indefn.GetFieldDefn(i)
        tmpname = fieldDef.GetNameRef()
        if outlayer.CreateField( fieldDef ) != 0:
            write_and_flush("Can't create field {}\n".format(tmpname))
            sys.exit(5)
        outfields.append(tmpname)

    # Create the new fields
    write_and_flush("Creating new fields for raster data:\n")
    for fname, ftype in newfields:
        fieldDef = ogr.FieldDefn( fname, ftype )
        if outlayer.CreateField( fieldDef ) != 0:
            write_and_flush("Can't create field {} of type {}\n".format(fname, ftype))
            sys.exit(5)
        outfields.append(fname)

    write_and_flush("Created output fields:\n{}\n".format(outfields))
    
    # End the batch of changes
    outlayer.CommitTransaction()

    #################################################
    ## Get raster data from all intersecting features
    #################################################
    
    # Make sure we start at the first feature
    inlayer.ResetReading()

    # A point that we can modify over and over
    tmppoint = ogr.Geometry(ogr.wkbPoint)

    # Start a batch of changes
    outlayer.StartTransaction()
    # Keep track of transactions
    num_trans = 0

    ## Iterate through features, grabbing data from maps 
    #####################################################
    for featnum, infeat in tqdm.tqdm(enumerate(inlayer),
                                    desc="Checking features ",
                                    ncols=80,
                                    total=infeatcount,
                                    mininterval=1):
        # Check for missing features
        if infeat is None:
            continue

        # Read in the feature's geometry 
        ingeom = infeat.GetGeometryRef()
        # Check that it is a geometry
        if ingeom is None:
            infeat.Destroy()
            print("feature {} had no geometry".format(featnum))
            sys.stdout.flush()
            continue
        

        # If we have a crs transformation, apply it
        if gt is not None:
            ingeom.Transform(gt)

        ## Use the first image to get XY points for this feature
        ########################################################
        #Convert the geometry into a shapely shape object
        shapely_geom = geom.shape(eval(ingeom.ExportToJson()))
        
        # Find which indices would be rasterized by this shape
        ######################################################
        # Get the window for this feature
        geom_extent = list(shapely_geom.buffer(firstrast.res[0]).bounds)
        
        # Convert this to a rasterio window and transform
        geom_wind = extent_to_window(geom_extent, firstrast)
        geom_trans = firstrast.window_transform(geom_wind)

        # Build a geometry mask for this subwindow
        geom_mask = feat.geometry_mask(
                [shapely_geom],
                (geom_wind.height, geom_wind.width),
                geom_trans,
                all_touched=args.alltouched,
                invert=True)

        # If no covered pixels, move on to next feature
        if geom_mask.sum() < 1:
            write_and_flush("Feature {} has no pixels inside raster {}\n".format(featnum+1,ordered_rastnames[0]))
            continue

        # Get reative row and index for the valid pixels in this mask
        mask_r, mask_c = numpy.where(geom_mask)

        # What are these with offsets added
        geom_r = mask_r + geom_wind.row_off
        geom_c = mask_c + geom_wind.col_off
        

        # Get the associated X and Y addresses of these indices,
        #  they will be used on all rasters to grab data
        geom_x, geom_y = firstrast.xy(geom_r, geom_c)
        
        # Buld a temporary dataset to store all the raster reads
        tmpdat = numpy.ones((len(newfields),len(geom_x)),dtype=numpy.float32) * numpy.nan
        # Pixel number
        tmpdat[1,:] = range(1,len(geom_r)+1)
        # Row, Col
        tmpdat[2,:] = geom_r
        tmpdat[3,:] = geom_c
        # X, Y
        tmpdat[4,:] = geom_x
        tmpdat[5,:] = geom_y

        #########################################################
        ## Gather data from the rasters at these X, Y coordinates
        #########################################################
        for rindex, rname in enumerate(raster_dict.keys()):
            rdict = raster_dict[rname]

            # Find the address of our geometry x,y in the current map
            feat_r, feat_c = rdict['Reference'].index(geom_x, geom_y)

            # Get a window from these addresses
            feat_window = indices_to_window(feat_r, feat_c)

            # Get the addresses relative to the window offset
            off_r = numpy.array(feat_r) - feat_window.row_off
            off_c = numpy.array(feat_c) - feat_window.col_off

            # Read the current raster within this window and using specified indices
            feat_dat = rdict['Reference'].read(
                    indexes=rdict['Bands'],
                    window=feat_window)

            # Copy the raster values to the tmpdat array
            try:
                tmpdat[rdict['Indices'],:] = feat_dat[:,off_r,off_c].reshape(feat_dat.shape[0],-1)
            except IndexError:
                tmpdat[rdict['Indices'],:] = numpy.nan
            
            # Convert nodata to NaN, unless ignore is specified
            if not args.ignore:
                for bi, b in enumerate(rdict['Bands']):
                    ii = rdict['Indices'][bi]
                    wnotval = numpy.where(tmpdat[ii,:] == rdict['Reference'].nodatavals[b-1])[0]
                    tmpdat[ii,wnotval]  = numpy.nan

        # Clear out any pixels whera all pixel data is nan 
        valpix = numpy.any(numpy.isfinite(tmpdat[metafieldcount:,:]),axis=0)
        tmpdat = tmpdat[:,valpix].T


        ##########################################################
        ## For each valid pixel, create a new output point feature
        ##########################################################
        num_added = 0
        for pix in range(tmpdat.shape[0]):

            # Create a new feature
            outfeat = ogr.Feature(outlayer.GetLayerDefn())

            # Copy over field from original feature
            for j in range(infieldcount):
                outfeat.SetField(j,infeat.GetField(j))

            # Set metadata fields
            outfeat.SetField(infieldcount,int(num_added+1))
            outfeat.SetField(infieldcount+1,int(tmpdat[pix,1]))
            outfeat.SetField(infieldcount+2,int(tmpdat[pix,2]))
            outfeat.SetField(infieldcount+3,int(tmpdat[pix,3]))
            outfeat.SetField(infieldcount+4,float(tmpdat[pix,4]))
            outfeat.SetField(infieldcount+5,float(tmpdat[pix,5]))

            # Set data fields
            for datfield in range(tmpdat.shape[1]-metafieldcount):
                outfeat.SetField(datfield+infieldcount+metafieldcount,float(tmpdat[pix,datfield+metafieldcount]))
                num_trans += 1

            # Update the point geometry
            ###########################
            tmppoint.AddPoint(float(tmpdat[pix,4]),float(tmpdat[pix,5]))

            # Translate back if needed
            if gt is not None:
                tmppoint.Transform(igt)

            # Add the geometry to the feature 
            outfeat.SetGeometry(tmppoint)
            num_added += 1

            # Add new feature to output
            if outlayer.CreateFeature(outfeat) != 0:
                warnings.warn("Can't create feature in output layer")
                sys.exit(6)

            outfeat = None

            # If we reached our buffer size, apply all waiting transactions
            if num_trans > 5000:
                outlayer.CommitTransaction()
                outlayer.StartTransaction()
                num_trans = 0

        # Clean up
        infeat = None
        ingeom = None

    ##One final commit 
    outlayer.CommitTransaction()

    outlayer = None

if __name__ == "__main__":
    main()
