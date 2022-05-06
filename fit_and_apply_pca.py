#!/usr/bin/env python
import os
import sys
import gdal
import numpy as np
import argparse
import pickle
import rasterio as rio
#import gc
#  import psutil
from rasterio.windows import Window

from sklearn.decomposition import IncrementalPCA


#  proc = psutil.Process(os.getpid())
#  def print_mem_usage():
#      print("Memory used by process: {} Mb".format(proc.memory_info()[0]/1024.0/1024))

#  def dump_garbage():
#      """
#      show us what's the garbage about
#      """
    
#      print("\nGLOBAL OBJECTS:")
#      cur = dict(globals())
#      for k, v in cur.items():
#          try:
#              print("{}: {} bytes".format(k,sys.getsizeof(v)))
#          except:
#              pass
        
#      print("\nLOCAL OBJECTS:")
#      cur = dict(locals())
#      for k, v in cur.items():
#          try:
#              print("{}: {} bytes".format(k,sys.getsizeof(v)))
#          except:
#              pass

#      # force collection
#      gc.collect()
#      print("\nGARBAGE OBJECTS:")
#      for x in gc.garbage:
#          s = str(x)
#          if len(s) > 80: s = s[:80]
#          print(s,":",sys.getsizeof(x)," - ",type(x))
#      del gc.garbage[:]


def parse_bands(bandstr):
    bands = []
    cparts = bandstr.split(",")
    for p in cparts:
        sp = p.strip()
        dparts = sp.split("-")
        if len(dparts) > 1:
            for tmpband in range(int(dparts[0]),int(dparts[1])+1):
                bands.append(tmpband)
        else:
            bands.append(int(dparts[0]))
    return bands

def main():
    parser = argparse.ArgumentParser(description="Compute PCA and/or apply transform to input raster")
    parser.add_argument("input_raster", 
            help = "Input raster map")

    parser.add_argument("out_raster", 
            help = "Outfile, computed as residual sum of squares")

    parser.add_argument("--pickled_model", "-p",
            default = '',
            help = "Pickled sklearn model that should be applied only (instead fitting model to input)")

    parser.add_argument("--model_out", "-m",
            default = '',
            help = "Save model to this file name, by default extension is changed on output raster name")
    
    parser.add_argument("--driver", "--format", "-of", "-f",
            default = 'GTiff',
            help = "GDAL format shortname of output map")
    
    parser.add_argument("--bands", "-b",
            default = '',
            help = "List of bands to use - comma separated with ranges, e.g. 1,4,7-10,12")
    
    parser.add_argument("--chunksize", "-c",
            type=float,
            default = 100,
            help = "How much data to read in per chunk (in Mb)")

    parser.add_argument("--co", "-co",
            action='append',
            metavar='KEY=VALUE',
            help = "GDAL creation options for specified format")

    parser.add_argument("--ignore", "-i",
            action='store_true',
            help = "Ignore nodata (consider pixels with nodata values in the transform)")

    parser.add_argument("--build_only", "-o",
            action='store_true',
            help = "Just build the PCA decomposition and save")

    parser.add_argument("--num_comp", "-n",
            type = int,
            default = 0,
            help = "Number of components to apply during transform, default all")

    args = parser.parse_args()


    ##Get information about input raster
    inmap = rio.open(args.input_raster, 'r')

    if args.bands != "":
        bandlist = parse_bands(args.bands)
    else:
        bandlist = [i+1 for i in range(inmap.count)]

    #Number of rows to read in at a time assuming 4-byte
    bytes_per_row = 4 * inmap.width * len(bandlist)
    chunk_size = min(max(1,int(np.floor((args.chunksize*1024*1024) / bytes_per_row))),inmap.height)

    do_fit = True
    if args.pickled_model != '':
        with open(args.pickled_model, "rb") as f:
            model = pickle.load(f)
        do_fit = False
        args.num_comp = min(args.num_comp,model.components_.shape[0],len(bandlist)) if args.num_comp > 0 else min(model.components_.shape[0],len(bandlist))
    else:
        args.num_comp = min(args.num_comp,len(bandlist)) if args.num_comp > 0 else len(bandlist)


    if args.model_out == '':
        args.model_out = os.path.splitext(args.out_raster)[0] + "_pcamod_{}comp".format(args.num_comp) + ".pkl"

    if inmap.nodata is None:
        args.ignore = True

    ##Make a dict from co arguments
    opts={}
    if args.co is not None:
        for v in args.co:
            eparts = v.split("=")
            if len(eparts) < 2:
                raise(ValueError("{} invalid for co argument - must be formatted 'KEY=STR'".format(v)))
            opts[eparts[0].strip()] = eparts[1].strip()


    inmap.close()

    if do_fit: 

        ##Set up the model
        model = IncrementalPCA(batch_size=chunk_size * inmap.width,copy=False)

        ##Set up chunking
        starts = np.arange(0,inmap.height,chunk_size)
        dat = None
        dat = None
        for s_i, s in enumerate(starts):
            tmpwind = Window(0,s,inmap.width,min(chunk_size,inmap.height-s))
            print("Reading in chunk {} of {}: {}".format(s_i + 1, len(starts), tmpwind))
            
            ##Get the data
            inmap = rio.open(args.input_raster, 'r')
            dat = inmap.read(indexes=bandlist,window=tmpwind).reshape(-1,tmpwind.height*tmpwind.width).astype(np.float32)
            inmap.close()

            ##Drop invalid values (and missing data if needed)
            if not args.ignore:
                #Make null into nan
                dat[dat == inmap.nodata] = np.nan
            wval = np.where(np.all(np.isfinite(dat),axis=0))[0]
            print("Dropping {} missing or invalid pixels".format(dat.shape[1]-len(wval)))
            dat = dat[:,wval]
            
            ##Fit this chunk
            if dat.shape[1] > 0:
                print("Adding {} samples to model".format(dat.shape[1]))
                model.partial_fit(dat.T)
                #  print_mem_usage()
                dat = None
            else:
                print("No valid pixels to add")
                dat = None
                continue

            ##Report changes
            print("Variance explained by first PC: {}".format(model.explained_variance_[0]))

        print("Done fitting model. Loadings: {}".format(model.explained_variance_))
    
    ##Save the model to a pickle
    print("Saving model to {}".format(args.model_out))
    with open(args.model_out,'wb') as pout:
        pickle.dump(model,pout)
    

    if not args.build_only: 
        ##Set up output map
        if not args.ignore:
            outmap = rio.open(
                    args.out_raster,
                    mode='w',
                    driver=args.driver,
                    width=inmap.width,
                    height=inmap.height,
                    count=args.num_comp,
                    crs=inmap.crs,
                    transform=inmap.transform,
                    dtype=rio.float32,
                    nodata=-9999,
                    **opts)
        else:
            outmap = rio.open(
                    args.out_raster,
                    mode='w',
                    driver=args.driver,
                    width=inmap.width,
                    height=inmap.height,
                    count=args.num_comp,
                    crs=inmap.crs,
                    transform=inmap.transform,
                    dtype=rio.float32,
                    **opts)


        ##Set up chunking
        starts = np.arange(0,inmap.height,chunk_size)
        for s_i, s in enumerate(starts):
            tmpwind = Window(0,s,inmap.width,min(chunk_size,inmap.height-s))
            print("Reading in chunk {} of {}: {}".format(s_i + 1, len(starts), tmpwind))
            
            ##Get the data
            inmap = rio.open(args.input_raster, 'r')
            dat = inmap.read(indexes=bandlist,window=tmpwind).reshape(-1,tmpwind.height*tmpwind.width).astype(np.float32)
            inmap.close()

            if not args.ignore:
                #Make null into nan
                dat[dat == inmap.nodata] = np.nan
            wval = np.where(np.all(np.isfinite(dat),axis=0))[0]
            print("Clearing {} missing or invalid pixels".format(dat.shape[1]-len(wval)))
            
            if args.ignore:
                out = np.zeros((args.num_comp,tmpwind.height*tmpwind.width))
            else:
                out = np.ones((args.num_comp,tmpwind.height*tmpwind.width)) * -9999

            ##Apply the transform
            print("Transforming data")
            out[:,wval] = model.transform(dat[:,wval].T)[:,:args.num_comp].T
            
            #  print_mem_usage()

            ##Write output
            print("Writing data")
            outmap.write(out.astype(np.float32).reshape(-1,tmpwind.height,tmpwind.width),window=tmpwind)

        print("Done writing pca map")

if __name__ == "__main__":
    main()
