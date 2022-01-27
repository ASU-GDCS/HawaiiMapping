# HawaiiMapping
Code for assisting with the Hawaii mapping efforts

## Scripts
### _sampleraster.py_
```
usage: sampleraster.py [-h] [--nodata NODATA] [--sql SQL] [--ignore]
                       [--alltouched] [--nln NLN] [--format FORMAT]
                       input output rasters [rasters ...]

Extract raster data from within vector points or regions

positional arguments:
  input                 Input vector file
  output                Output vector file
  rasters               Entries of the format FIELDNAME=FILENAME[:BANDLIST],
                        where optional BANDLIST is 1-indexed, and can be a
                        list separated by commas and can include ranges (i.e,
                        1:55). All rasters should bin in same spatial
                        reference. If vector and at least one raster have crs
                        in metainfo, then coordinate transformation of vector
                        features will occur

optional arguments:
  -h, --help            show this help message and exit
  --nodata NODATA, -nodata NODATA
                        Specifiy no data value in raster (if not identified in
                        metadata)
  --sql SQL, -sql SQL   Apply sql query to vector layer before collecting
                        stats
  --ignore, -ignore, -i
                        Ignore raster no data value - include pixels marked as
                        invalid
  --alltouched, -alltouched, -a
                        Include all touched pixels in areal vectors (default
                        is just those with center inside)
  --nln NLN, -nln NLN   Name out output layer - for formats that support this
  --format FORMAT, -format FORMAT, -of FORMAT
                        OGR format of output file - default is same as input
                        vector file

Example:
sampleraster.py sequoia_crowns.shp --format GPKG sequoia_crowns_filled.gpkg refl=seki_refl elev=seki_elev.tif shade=seki_shade mbdat=multiband.tif:1,3,5-6
