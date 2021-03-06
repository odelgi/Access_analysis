#-----------------------------------------------------------------------------------------------------------------------
# Analysis_Chp1 (2020)
# Title: Access analysis data preparation
# del Giorgio, Olivia (https://github.com/odelgi)
# Messager, Mathis (https://github.com/messamat)

#Structure of analysis
# A. SET ENVIRONMENT SETTINGS
# B. SET DIRECTORY STRUCTURE
# C. DEFINE INPUT VARIABLES
# D. DEFINE OUTPUT VARIABLES
# E. PREPARE DATA
#    1. Get and pre-process Landsat imagery
#    2. Project and rasterize polygon of Pelegrinni department (matching grid layout of Landsat)
#    3. Download, project, and clip Hansen et al. forest loss
#    4. Compute deforestation raster for each year
#    5. Rasterize barrier lines
#    6. Assign livelihood-specific weighting to barrier rasters
#    7. For each livelihood, assign pixels to groups whose livelihood-specific buffer sizes do not overlap so that cost
#    distance can be calculated for multiple locations at a time.

#-----------------------------------------------------------------------------------------------------------------------

# Import system modules
import arcpy
from arcpy.sa import *

import os #Module for operating system operations (e.g. os.path.exists)
import re #Module for string pattern matching
import requests #Module to access URLs
import pandas as pd #Module for manipulating tables (non-gis tables mostly) https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html
from usgs import api #Module to grab landsat data
from collections import defaultdict #Module to define dictionaries (data structure) with default values (e.g. list)
import json
import urllib2
import time
import gzip
import tarfile
from functools import wraps
from accesscalc_parallel import *

###---------------------------------------- A. DEFINE FUNCTIONS -----------------------------------------------------###
#Function to retry three times if urllib2.urlopen fails
def retry(ExceptionToCheck, tries=3, delay=3, backoff=2, logger=None):
    """Retry calling the decorated function using an exponential backoff.
    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry
    :param ExceptionToCheck: the exception to check. may be a tuple of
        exceptions to check
    :type ExceptionToCheck: Exception or tuple
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    """
    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck, e:
                    msg = "%s, Retrying in %d seconds..." % (str(e), mdelay)
                    if logger:
                        logger.warning(msg)
                    else:
                        print msg
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry

@retry(urllib2.URLError, tries=3, delay=3, backoff=2)
def urlopen_with_retry(in_url):
    return urllib2.urlopen(in_url)#Retry three times if urllib2.urlopen fails

#Function to download and unzip landsat scenes
def dl_landsattiles(lss, dataset, apikey, outdir, mincctile=True, maxcloudcover=100):
    sublss = defaultdict(dict)
    print('Subsetting scenes to keep those with minimum land cover...')
    if mincctile:
        for tile in lss['data']['results']:
            rowpath = tile['displayId'].split('_')[2]
            if rowpath in sublss:
                if float(tile['cloudCover']) < float(sublss[rowpath]['cloudCover']):
                    sublss[rowpath] = tile
            else:
                sublss[rowpath] = tile
    else:
        sublss = lss

    for tile in sublss.values():
        # Only continue if there is less than set cloud cover on the image
        if float(tile['cloudCover']) < maxcloudcover:
            print(tile['entityId'])
            # Get tile download info
            gettile = api.download(dataset=dataset,
                                   node='EE',
                                   entityids=tile['entityId'],
                                   product='STANDARD',
                                   api_key=apikey)
            tileurl = gettile['data'][0]['url']

            # Get file name
            print('Trying to access {}...'.format(tileurl))
            req = urlopen_with_retry(in_url=tileurl)
            cd = req.headers.get('content-disposition')
            fname = re.findall('(?<=filename=")[A-Za-z0-9_.]+', cd)[0]
            tiledirname = os.path.join(outdir, os.path.splitext(os.path.splitext(fname)[0])[0])
            if not os.path.exists(tiledirname):
                os.mkdir(tiledirname)
            outpath = os.path.join(tiledirname, fname)

            if not os.path.exists(outpath):  # Inspired from https://github.com/yannforget/landsatxplore
                with open(outpath, 'wb') as f:
                    chunk = req.read()
                    f.write(chunk)
                print('Downloaded {}, now unzipping...'.format(outpath))

                tarf = os.path.splitext(outpath)[0]
                with gzip.GzipFile(outpath, 'rb') as input:
                    s = input.read()
                    with open(tarf, 'wb') as output:
                        output.write(s)

                with tarfile.open(tarf) as ftar:
                    ftar.extractall(tiledirname)  # specify which folder to extract to
            else:
                print('{} already exists...'.format(outpath))

###---------------------------------------- A. SET ENVIRONMENT SETTINGS ---------------------------------------------###
#https://pro.arcgis.com/en/pro-app/arcpy/classes/env.htm
arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

analysis_years = ['2000', '2010', '2018'] #Years for analysis

###---------------------------------------- B. SET DIRECTORY STRUCTURE ----------------------------------------------###
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0] #Get project directory based on location of script within rootdir/src/
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

pathsdir =  os.path.join(datadir, 'Paths') #Primary data
data2dir = os.path.join(datadir, 'Secondary_data') #Secondary data

barriersdir = os.path.join(datadir, 'Barriers_Pellegrini')
pathcheckcreate(barriersdir)

#Landsat imagery directory
landsatdir = os.path.join(datadir, 'landsat') #Raw data
pathcheckcreate(landsatdir)

landsatoutdir = os.path.join(resdir, 'landsat') #Processed data
pathcheckcreate(landsatoutdir)

#Global Forest Change directory
forestdir = os.path.join(data2dir, 'Forest_Hansen')
pathcheckcreate(forestdir)
forestoutdir = os.path.join(resdir, 'Base_layers_Pellegrini/Forest_Hansen.gdb')
pathcheckcreate(forestoutdir)

# Source file geodatabase for unweighted barrier features
barrier_indir = os.path.join(barriersdir, "Barriers_clean_coded.gdb")
pathcheckcreate(barrier_indir)

# Output geodatabase to hold basemaps
basemapgdb = os.path.join(resdir,
                          "Base_layers_Pellegrini",
                          "Basemaps_UTM20S.gdb")
pathcheckcreate(basemapgdb)

# Output raster geodatabase for unweighted barrier rasters
barrier_outdir = os.path.join(resdir, "Base_layers_Pellegrini\Boundary_rasters_3030.gdb")
pathcheckcreate(barrier_outdir)

# Output directory for livelihood-specific weighted barrier raster
barrierweight_outdir = os.path.join(resdir,
                                    'Analysis_Chp1_W1',
                                    'W1_3030',
                                    'Barrier_weighting_1_3030')
pathcheckcreate(barrierweight_outdir)

#Output gdb to hold intermediate data product for access calculation
llhoodbuffer_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Buffers_W1.gdb') #Output buffer path for each livelihood
pathcheckcreate(llhoodbuffer_outdir) #Make sure that directories exist to write out buffers

###---------------------------------------- C. DEFINE INPUT VARIABLES -----------------------------------------------###
"""Potential source of data: http://demo-ide.arsat.com.ar/ide-santiago/"""
#UTM20S Coordinate System (CS) oject to use when projecting
csref = arcpy.SpatialReference(32720)
pelleoutline = os.path.join(data2dir,  'Pellegrini_outline.gdb/Pellegri_department_outline_wgs84')
pelleextent_wgs84 = arcpy.Describe(pelleoutline).Extent #command to get layer's metadata

#Table containing data on livelihood-specific use area radius and barrier weighting
weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')

###---------------------------------------- D. DEFINE OUTPUT VARIABLES ----------------------------------------------###
pelleUTM20S = os.path.join(basemapgdb, "Pellegri_department_outline_UTM20S")
pelleras = os.path.join(basemapgdb,"Pellegri_department_UTM20S") #Reference raster for snapping throughout analysis
pellepoints = os.path.join(llhoodbuffer_outdir, 'pellefishpoints')
fishpoints = os.path.join(resdir, 'Analysis_Chp1_W1', 'Buffers_W1.gdb', 'fishnetpoints')

#Create dictionary that associates each year (key) with the path of the barrier layer for that year (value)
barrier_cleancodedlist = {i: os.path.join(barrier_indir, 'Paths_{}_coded_clean_UTM20S'.format(i))
                          for i in analysis_years}

###---------------------------------------- E. PREPARE DATA ---------------------------------------------------------###
#-----------------------------------------------------------------------------------------------------------------------
### 1. Get and pre-process Landsat imagery
# Initialize a new API instance and get an access key
# The user credentials that will be used to authenticate access to the data
with open("configs.json") as json_data_file:  # https://martin-thoma.com/configuration-files-in-python/
    authdat = json.load(json_data_file)
#Get temporary API key
usgs_api_key = api.login(
        str(authdat["username"]),
        str(authdat["password"]),
        save=False,
        catalogId='EE')['data']

#Get list of scenes from Landsat 7 for 2000
lss_LC7_2000 = api.search(dataset='LANDSAT_ETM_C1',
                          node='EE',
                          ll={"longitude":pelleextent_wgs84.XMin, "latitude":pelleextent_wgs84.YMin},
                          ur={"longitude":pelleextent_wgs84.XMax, "latitude":pelleextent_wgs84.YMax},
                          start_date='2000-01-01', end_date='2000-12-31',
                          api_key = usgs_api_key)

#Download and unzip landsat 7 scenes with least cloud cover for 2000
dl_landsattiles(lss=lss_LC7_2000, dataset='LANDSAT_ETM_C1', apikey=usgs_api_key,
                outdir=landsatdir, mincctile=True, maxcloudcover=100)

#Get list of scenes from Landsat 7 for 2010
lss_LC7_2010 = api.search(dataset='LANDSAT_ETM_C1',
                          node='EE',
                          ll={"longitude":pelleextent_wgs84.XMin, "latitude":pelleextent_wgs84.YMin},
                          ur={"longitude":pelleextent_wgs84.XMax, "latitude":pelleextent_wgs84.YMax},
                          start_date='2010-01-01', end_date='2010-12-31',
                          api_key = usgs_api_key)

#Download and unzip landsat 7 scenes with least cloud cover for 2010
dl_landsattiles(lss=lss_LC7_2010, dataset='LANDSAT_ETM_C1', apikey=usgs_api_key,
                outdir=landsatdir, mincctile=True, maxcloudcover=100)

#Get list of scenes from Landsat 8
lss_LC8 = api.search(dataset='LANDSAT_8_C1',
                     node='EE',
                     ll={"longitude":pelleextent_wgs84.XMin, "latitude":pelleextent_wgs84.YMin},
                     ur={"longitude":pelleextent_wgs84.XMax, "latitude":pelleextent_wgs84.YMax},
                     start_date='2018-01-01', end_date='2018-12-31',
                     api_key = usgs_api_key)

#Download and unzip landsat 8 scenes with least cloud cover for 2018
dl_landsattiles(lss=lss_LC8, dataset='LANDSAT_8_C1', apikey=usgs_api_key,
                outdir=landsatdir, mincctile=True, maxcloudcover=100)

#Create composites
for scenedir in os.listdir(landsatdir):
    scene_composite = os.path.join(landsatoutdir, '{}_composite.tif'.format(scenedir))
    if not arcpy.Exists(scene_composite):
        print('Create {}...'.format(scene_composite))
        scenebandslist = getfilelist(dir=os.path.join(landsatdir, scenedir),
                                     repattern='(LC08|LE07)_L1[GTP]{2}.*B[0-9]+[.]TIF$')
        arcpy.CompositeBands_management(scenebandslist, scene_composite)

#Mosaick and project scenes
landsat_raw = getfilelist(dir=landsatoutdir, repattern='(LC08|LE07)_L1[GTP]{2}.*_composite[.]tif$') #For explanation of repattern, check https://docs.python.org/3/howto/regex.html
landsat_years = set([os.path.split(tilename)[1][17:21] for tilename in landsat_raw])
#os.path.split splits path into a two-item list with path root and file name as the first and second objects in the list
#[1] grabs the second item in that list (file name)
#[17:21] subset the file name string to grab from the 18th to the 22nd character
#for tilename in landsat raw loops through the objects in landsat_raw and in each loop assign that object to tilename variable
#set takes unique values in list

landsatproj = {}
for year in landsat_years:
    #if year in analysis_years:
    inraslist = [tilename for tilename in landsat_raw if os.path.split(tilename)[1][17:21] == year] #Grab only images from 'year'
    outras = os.path.join(landsatoutdir, 'landsatmosaic_{}.tif'.format(year))

    if not os.path.exists(os.path.join(landsatoutdir, outras)):
        print('Moisaicking landsat images from {}...'.format(year))
        arcpy.MosaicToNewRaster_management(input_rasters=inraslist,
                                           output_location= os.path.split(outras)[0],
                                           raster_dataset_name_with_extension=os.path.split(outras)[1],
                                           pixel_type='16_BIT_UNSIGNED',
                                           number_of_bands=arcpy.Describe(inraslist[0]).bandCount,
                                           mosaic_method='MAXIMUM')

    #Re-project to UTM 20S
    landsatproj[year] = '{}_UTM20S.tif'.format(os.path.splitext(outras)[0]) #os.path.splitext() method in Python is used to split the path name into a pair root and ext. Here, ext stands for extension
    if not os.path.exists(landsatproj[year]):
        print('Projecting landsat images from {}...'.format(year))
        arcpy.ProjectRaster_management(in_raster=outras,
                                       out_raster= landsatproj[year],
                                       out_coor_system = csref,
                                       resampling_type='NEAREST')

#-----------------------------------------------------------------------------------------------------------------------
### 2. Project and rasterize polygon of Pelegrinni department (matching grid layout of Landsat)
if not arcpy.Exists(pelleUTM20S):
    print('Projecting study area polygon...')
    arcpy.Project_management(in_dataset=pelleoutline,
                             out_dataset=pelleUTM20S,
                             out_coor_system=csref)

pelleextent_utm = arcpy.Describe(pelleUTM20S).Extent #command to get layer's metadata
refraster = landsatproj['2020']
arcpy.env.snapRaster = refraster
if not arcpy.Exists(pelleras):
    print('Converting study area polygon to raster...')
    arcpy.FeatureToRaster_conversion(in_features=pelleUTM20S, field="Code",
                                     out_raster=pelleras,
                                     cell_size=refraster)
arcpy.env.snapRaster = pelleras #Reference raster for snapping throughout analysis #####################################

#-----------------------------------------------------------------------------------------------------------------------
### 3. Download, project, and clip Hansen et al. forest loss ###
GFCoutclip = {}
for product in ['datamask', 'gain', 'lossyear', 'treecover2000']:

    #Download
    GFCURL = "https://storage.googleapis.com/earthenginepartners-hansen/" \
             "GFC-2018-v1.6/Hansen_GFC-2018-v1.6_{}_20S_070W.tif".format(product)

    GFCout = os.path.join(forestdir, 'Hansen_GFC_v16_{}.tif'.format(product))

    if not os.path.exists(GFCout):
        print("downloading {}".format(GFCURL))
        with open(GFCout, "wb") as local_file:
            local_file.write(requests.get(GFCURL, stream=True).content)
    else:
        print('{} already exists...'.format(GFCout))

    GFCoutproj = os.path.join(forestoutdir, 'Hansen_GFC_v16_{}proj'.format(product))

    #Project to UTM 20S
    if not arcpy.Exists(GFCoutproj):
        print("projecting {}".format(GFCout))
        arcpy.ProjectRaster_management(in_raster=GFCout,
                                       out_raster= GFCoutproj,
                                       out_coor_system = csref,
                                       resampling_type='NEAREST',
                                       cell_size=arcpy.Describe(refraster).meanCellHeight)

    #Clip to study area
    GFCoutclip[product] = os.path.join(forestoutdir, 'Hansen_GFC_v16_{}clip'.format(product))
    if not arcpy.Exists(GFCoutclip[product]):
        print("clipping {}".format(GFCoutclip[product]))
        ExtractByMask(in_raster=GFCoutproj,
                      in_mask_data=pelleras).save(GFCoutclip[product])

#-----------------------------------------------------------------------------------------------------------------------
### 4. Compute deforestation raster for each year

# Compute yearly forest cover
forestyearly = {}
for year in analysis_years:
    #Assign output path to foresyearly dictionary with year as the key
    forestyearly[year] = os.path.join(forestoutdir, 'Hansen_GFC_v16_treecover{}'.format(year))

    if not arcpy.Exists(forestyearly[year]): #First check that output does not exist
        print('Produce {}...'.format(forestyearly[year]))

        forestloss_sum = Con((Raster(GFCoutclip['lossyear']) < int(year[2:])) & #If loss year < year (0-18)
                             (Raster(GFCoutclip['lossyear']) > 0),  # and loss year > 0
                             0, #Assign 0 to output (no forest)
                             Con(Raster(GFCoutclip['treecover2000']) >= 1, 1, 0)) #Otherwise, if treecover2000 > 1%, assign 1 (forest), else assign 0 (no forest)
        forestloss_sum.save(forestyearly[year])

#-----------------------------------------------------------------------------------------------------------------------
### 5. Rasterize barrier lines
#For each year that we decided to analyze
barrier_outras = {} #Create empty dictionary to populate with output raster layer names in the loop associated to each year
for year in analysis_years:
    #Create name of output raster (value) for the given year (key)
    barrier_outras[year] = os.path.join(barrier_outdir,
                                        "{}ras".format(os.path.split(barrier_cleancodedlist[year])[1]))

    if not arcpy.Exists(barrier_outras[year]):
        print('Rasterize {}...'.format(barrier_cleancodedlist[year]))
        arcpy.PolylineToRaster_conversion(in_features=barrier_cleancodedlist[year],
                                          value_field='Code',
                                          out_rasterdataset=barrier_outras[year],
                                          priority_field='Code',
                                          cellsize=refraster)

#-----------------------------------------------------------------------------------------------------------------------
### 6. Assign livelihood-specific weighting to barrier rasters
weightingmeta = pd.read_excel(weighting_table, sheetname='Metadata') #Read weighting excel table, only Metadata sheet
weightingcodes = weightingmeta.loc[weightingmeta['Variable name'].str.contains("BW"), #Subset metadata sheet to only keep rows for which 'Variable name' column string contains 'BW'
                                   ['Variable name', 'Code']] #For those rows, only keep Variable name and Code columns
weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

livelihoods = weightingpd['Livelihood'].tolist()
livelihoods.remove('Combined_livelihood')

barrierweight_outras = {} #Dictionary that will be populated with the path to the weighted barrier data for each livelihood-year combination
#Iterate through the different livelihood names (as contained in the excel sheet)
for llhood in livelihoods:
    print(llhood)

    #Create the gdb required for weighted barriers for each livelihood
    outllhood_gdb = os.path.join(resdir,
                                 'Analysis_Chp1_W1',
                                 'W1_3030',
                                 'Barrier_weighting_1_3030',
                                 '{}_bw1_3030.gdb'.format(llhood))
    pathcheckcreate(outllhood_gdb)

    #Iterate through years of analysis
    for year in analysis_years:
        #Create path
        barrierweight_outras[llhood+year] = os.path.join(outllhood_gdb, '{0}_bw1_{1}'.format(llhood, year))

        if not arcpy.Exists(barrierweight_outras[llhood+year]):
            print('Reclassifying, output raster:{}...'.format(barrierweight_outras[llhood+year]))

            #Subselect weighting code table to only keep values for the livelihood of interest
            subtable = weightingpd.loc[weightingpd['Livelihood']==llhood]
            #Create a reclassify template (list) with, for each barrier type, a list whose first element is the barrier code and the second element is the barrier weight
            reclastemplate = [[int(row['Code']), int(10*subtable[row['Variable name']])] for index, row in weightingcodes.iterrows()]
            reclastemplate.append(['NODATA',0])
            outreclas = ExtractByMask(
                in_raster= arcpy.sa.Reclassify(in_raster=barrier_outras[year],
                                               reclass_field="Value",
                                               remap=arcpy.sa.RemapValue(reclastemplate),
                                               missing_values='NODATA'),
                in_mask_data=pelleras)
            outreclas.save(barrierweight_outras[llhood+year])

#-----------------------------------------------------------------------------------------------------------------------
### 7. Assign livelihood-specific weighting to barrier rasters

#Create gdbs for cost raster
costtab_outgdb = {}
for llhood in livelihoods:
    #Output cost distance
    for year in analysis_years:
        costtab_outgdb[llhood+year] = os.path.join(resdir, 'Analysis_Chp1_W1',
                                                   'W1_3030', 'Cost_distance_W1_3030',
                                                   'Cost_distance_{}_w1'.format(llhood),
                                                   'CD_{0}_{1}_w1.gdb'.format(llhood, year))
        pathcheckcreate(costtab_outgdb[llhood+year])
del llhood

#For each livelihood, assign pixels to groups so that can run non-overlapping cost distance within livelihood-specific buffer sizes
bufferad = {}
for llhood in livelihoods:
    print(llhood)
    bufferad[llhood] = float(weightingpd.loc[weightingpd['Livelihood']==llhood, 'Buffer_max_rad']) #Get livelihood_specific buffer radius from table
    fishout = os.path.join(resdir, 'Analysis_Chp1_W1', 'Buffers_W1.gdb', 'fishnet{}'.format(llhood))
    arcpy.env.outputCoordinateSystem = csref

    if not arcpy.Exists(fishout):
        print('Create fishnet for {}...'.format(llhood))
        arcpy.CreateFishnet_management(out_feature_class = fishout,
                                       origin_coord = '{0} {1}'.format(pelleextent_utm.XMin, pelleextent_utm.YMin),
                                       y_axis_coord = '{0} {1}'.format(pelleextent_utm.XMin, pelleextent_utm.YMax),
                                       cell_width = (bufferad[llhood]-bufferad[llhood]%arcpy.Describe(refraster).meanCellWidth)*2 + arcpy.Describe(refraster).meanCellWidth*5,
                                       cell_height = (bufferad[llhood]-bufferad[llhood]%arcpy.Describe(refraster).meanCellHeight)*2 + arcpy.Describe(refraster).meanCellHeight*5,
                                       labels = 'NO_LABELS',
                                       template = pelleextent_utm,
                                       geometry_type = 'POLYGON')

    fishras = os.path.join(resdir, 'Analysis_Chp1_W1',
                           'Buffers_W1.gdb',
                           'fishnetras{}'.format(llhood))

    #Change raster tile size to make sure that cell IDs are ordered within a fishnet cell
    arcpy.env.tileSize = "{0} {0}".format(int(((bufferad[llhood] - bufferad[llhood] % arcpy.Describe(refraster).meanCellWidth) * 2 +
                                               arcpy.Describe(refraster).meanCellWidth * 5) / arcpy.Describe(refraster).meanCellWidth))

    if not arcpy.Exists(fishras):
        arcpy.env.extent = fishout
        print('Converting fishnet to raster for {}...'.format(llhood))
        arcpy.PolygonToRaster_conversion (in_features = fishout,
                                          value_field = arcpy.Describe(fishout).OIDFieldName,
                                          out_rasterdataset = fishras,
                                          cell_assignment = 'CELL_CENTER',
                                          cellsize = refraster)

    #Convert study area raster to point (one point for each 30x30 m cell)
    if not arcpy.Exists(fishpoints):
        print('Converting fishnet raster back to points...')
        arcpy.RasterToPoint_conversion(fishras, fishpoints, "Value")

        #Add new field to point layer whose value will be 0 if point falls within study area boundaries and <null> otherwise
        print('Overlap points with Pellegrinni outline')
        ExtractMultiValuesToPoints(in_point_features=fishpoints,
                                   in_rasters= [[pelleras, "pelleoverlap"]],
                                   bilinear_interpolate_values="NONE")
        arcpy.AlterField_management(fishpoints,
                                    field='grid_code',
                                    new_field_name='gc_{}'.format(llhood),
                                    new_field_alias='gc_{}'.format(llhood))
    else:
        if not 'gc_{}'.format(llhood) in [f.name for f in arcpy.ListFields(fishpoints)]:
            print('Extract fishnet cell number to points for {}...'.format(llhood))
            ExtractMultiValuesToPoints(in_point_features=fishpoints,
                                       in_rasters= [[fishras, 'grid_code']],
                                       bilinear_interpolate_values="NONE")
            arcpy.AlterField_management(fishpoints,
                                        field='grid_code',
                                        new_field_name='gc_{}'.format(llhood),
                                        new_field_alias='gc_{}'.format(llhood))
    #Assign group to each point
    print('Compute point coordinates...')
    if not ('POINT_X' in [f.name for f in arcpy.ListFields(fishpoints)]):
        arcpy.AddGeometryAttributes_management(Input_Features=fishpoints,
                                               Geometry_Properties='POINT_X_Y_Z_M')

    if not 'group{}'.format(llhood) in [f.name for f in arcpy.ListFields(fishpoints)]: #Check whether livelihood-specific grouping field already exists in fish points' attribute table
        print('Assign non-overlapping group number to points for {}...'.format(llhood))
        arcpy.AddField_management(in_table=fishpoints,
                                  field_name='group{}'.format(llhood),
                                  field_type='LONG')
        groupdict = defaultdict(lambda: 1) #Create a dictionary for which key corresponds to the fishnet cell number and the default value is 1
        with arcpy.da.UpdateCursor(in_table=fishpoints,
                                   field_names=['gc_{}'.format(llhood), 'group{}'.format(llhood), 'POINT_X', 'POINT_Y'],
                                   sql_clause = (None, 'ORDER BY POINT_X, POINT_Y')) as cursor: #Create a cursor to query attribute table of fishpoint
            x=0
            for row in cursor: #Iterate through points first column-wise, then row-wise
                if x % 100000 == 0: #Homemade progress bar: print row number every 100,000th row
                    print(x)
                row[1] = groupdict[row[0]] #Assign to 'group{}'.format(llhood) field the dictionary value based on which fishnet cell the point is in (expressed by 'gc_llhood' in fishpoints)
                groupdict[row[0]] +=1  #Add one to group value in dictionary for that fishnet cell (so that the next point that fall within that cell gets the same group number + 1)
                cursor.updateRow(row) #Write value to table
                x+=1

    arcpy.env.extent = pelleextent_utm

# Subset points to only keep those that fall within the department of Pellegrini
if not arcpy.Exists(pellepoints):
    arcpy.MakeFeatureLayer_management(in_features=fishpoints, out_layer='pellepoints_lyr',
                                      where_clause='pelleoverlap = 0')
    arcpy.CopyFeatures_management(in_features='pellepoints_lyr', out_feature_class=pellepoints)

#-----------------------------------------------------------------------------------------------------------------------
#Pre-process areas for assessing access in specific communities
commugdb = os.path.join(datadir, 'Community_boundaries' , 'Developed_areas.gdb')
commuptsgdb = os.path.join(datadir, 'Community_boundaries', 'Community_locations_byyear.gdb')
commupolys = getfilelist(commugdb, repattern='Developed_areas_.*', gdbf=True, nongdbf=False)
commu_censopts = getfilelist(commuptsgdb, repattern='Ubicacion_censo_UTM20S_.*', gdbf=True, nongdbf=False)

chunkdir_commu = os.path.join(resdir, 'Analysis_Chp1_W1', 'inputdata_communities')
chunkdir_commudat = os.path.join(chunkdir_commu, 'datapts.gdb')
pathcheckcreate(chunkdir_commudat)
commumerge = os.path.join(chunkdir_commudat, 'developed_areas_merge')
commu_censoptsmerge = os.path.join(chunkdir_commudat, 'Ubicacion_censo_UTM20S_merge')
commudiss = os.path.join(chunkdir_commudat, 'developed_areas_dissolve')
commupts = os.path.join(chunkdir_commudat, 'pellepoints_developedareas')

if not arcpy.Exists(commupts):
    arcpy.Merge_management(commupolys, commumerge) #Merge all developed areas across years (i.e. one layer with all separate polygons for all years)
    arcpy.Merge_management(commu_censopts, commu_censoptsmerge) #Same for census community point locations

    #Subset developed area polygons to onlly keep those that intersect a census community loation
    arcpy.MakeFeatureLayer_management(commumerge, 'commumerge_lyr')
    arcpy.SelectLayerByLocation_management('commumerge_lyr',
                                           overlap_type='INTERSECT',
                                           select_features=commu_censoptsmerge)

    #Dissolve developed area polygons to have a single layer without edges (removing overlapping poly)
    arcpy.Dissolve_management(in_features='commumerge_lyr', out_feature_class=commudiss)

    #Subset pixel/fishnet points to keep only those intersecting community developed areas
    arcpy.MakeFeatureLayer_management(pellepoints, 'pelleptlyr')
    arcpy.SelectLayerByLocation_management('pelleptlyr',
                                           overlap_type='INTERSECT',
                                           select_features=commudiss)
    arcpy.CopyFeatures_management('pelleptlyr', commupts)

#-----------------------------------------------------------------------------------------------------------------------
hhpts = os.path.join(datadir, 'Household_points/hh_points.gdb/hh_points')
chunkdir_hhdat = os.path.join(resdir, 'Analysis_Chp1_W1', 'inputdata_hh', 'datapts.gdb')
hhpts_proj = os.path.join(chunkdir_hhdat, 'hh_pointsproj')
hhras = os.path.join(chunkdir_hhdat, 'hh_pointsras')
hhraspts = os.path.join(chunkdir_hhdat, 'hh_pointraspts')
pathcheckcreate(chunkdir_hhdat)
pellehhpts = os.path.join(chunkdir_hhdat, 'pellepoints_households')

if not arcpy.Exists(hhpts):
    arcpy.env.snapRaster = pelleras
    arcpy.Project_management(hhpts, out_dataset=hhpts_proj, out_coor_system=pelleras)
    arcpy.PointToRaster_conversion(in_features=hhpts_proj, value_field='participant_code',
                                   out_rasterdataset=hhras, cellsize=pelleras)
    arcpy.RasterToPoint_conversion(hhras, hhraspts)

    arcpy.MakeFeatureLayer_management(pellepoints, 'pelleptlyr')
    arcpy.SelectLayerByLocation_management('pelleptlyr',
                                           overlap_type='INTERSECT',
                                           select_features=hhraspts)
    arcpy.CopyFeatures_management('pelleptlyr', pellehhpts)



