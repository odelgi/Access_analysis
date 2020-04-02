#------------------------------------------------------------------------------
# Analysis_Chp1 (2020)
# Title: Access analysis
# del Giorgio, Olivia
# Messager, Mathis (https://github.com/messamat)

# Model description: Produces a raster with "access index" values for each cell

"""Usage: 
where

"""
# -----------------------------------------------------------------------------
# Import system modules
import arcpy
from arcpy.sa import *
import numpy as np
import os #Module for operating system operations (e.g. os.path.exists)
import re #Module for string pattern matching
import requests #Module to access URLs
import pandas as pd #Module for manipulating tables (non-gis tables mostly) https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html
from usgs import api #Module to grab landsat data
from collections import defaultdict
from collections import Counter
import time
import math
from internal_functions import *

### Set environment settings ###
#https://pro.arcgis.com/en/pro-app/arcpy/classes/env.htm
arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

# File extension for the output rasters
# output_ext = ".img"

analysis_years = ['2000', '2010', '2018'] #Years for analysis

#### Directory structure ###
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

data2dir = os.path.join(datadir, 'Secondary_data') #Secondary data
pathsdir =  os.path.join(datadir, 'Paths') #Primary data

barriersdir = os.path.join(datadir, 'Barriers_Pellegrini')
pathcheckcreate(barriersdir)

#Landsat imagery directory
landsatdir = os.path.join(datadir, 'landsat')
pathcheckcreate(landsatdir)

landsatoutdir = os.path.join(resdir, 'landsat')
pathcheckcreate(landsatoutdir)

#Global Forest Change directory
forestdir = os.path.join(data2dir, 'Forest_Hansen')
pathcheckcreate(forestdir)
forestoutdir = os.path.join(resdir, 'Base_layers_Pellegrini/Forest_Hansen.gdb')
pathcheckcreate(forestoutdir)

# Source file geodatabase for unweighted barrier features
barrier_indir = os.path.join(barriersdir, "Barriers_clean_coded.gdb")
pathcheckcreate(barrier_indir)

# Output raster geodatabase for unweighted barrier rasters
barrier_outdir = os.path.join(resdir, "Base_layers_Pellegrini\Boundary_rasters_3030.gdb")
pathcheckcreate(barrier_outdir)

# Output directory for livelihood-specific weighted barrier raster
barrierweight_outdir = os.path.join(resdir,
                                    'Analysis_Chp1_W1',
                                    'W1_3030',
                                    'Barrier_weighting_1_3030')
pathcheckcreate(barrierweight_outdir)

#Output buffer and points gdbs
llhoodbuffer_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Buffers_W1.gdb') #Output buffer path for each livelihood
pathcheckcreate(llhoodbuffer_outdir) #Make sure that directories exist to write out buffers

#Output director for HPC directories containing buffers and points to process through cost distance
inputdir_HPC = os.path.join(resdir, 'Analysis_Chp1_W1', 'inputdata_HPC')
pathcheckcreate(inputdir_HPC)

#### Input variables #### 
"""Potential source of data: http://demo-ide.arsat.com.ar/ide-santiago/"""
#UTM20S Coordinate System (CS) oject to use when projecting
csref = arcpy.SpatialReference(32720)
pelleoutline = os.path.join(data2dir,
                            'Pellegrini_outline.gdb/Pellegri_department_outline_wgs84')
pelleextent_wgs84 = arcpy.Describe(pelleoutline).Extent #command to get layer's metadata
weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')

### Output variables ####
basemapgdb = os.path.join(resdir,
                          "Base_layers_Pellegrini",
                          "Basemaps_UTM20S.gdb")
pathcheckcreate(basemapgdb)

pelleUTM20S = os.path.join(basemapgdb, "Pellegri_department_outline_UTM20S")
pelleras = os.path.join(basemapgdb,"Pellegri_department_UTM20S") #Reference raster for snapping throughout analysis
pellepoints = os.path.join(llhoodbuffer_outdir, 'pellefishpoints')
fishpoints = os.path.join(resdir, 'Analysis_Chp1_W1', 'Buffers_W1.gdb', 'fishnetpoints')

#Create dictionary that associates each year (key) with the path of the barrier layer for that year (value)
barrier_cleancodedlist = {i: os.path.join(barrier_indir, 'Paths_{}_coded_clean_UTM20S'.format(i))
                          for i in analysis_years}

###################### PREPARE DATA ###########################################
### Get Landsat imagery ###
# Initialize a new API instance and get an access key
#### TO CONTINUE -- NEED API KEY ######
"""api.search(dataset='LANDSAT_8', node='EE', 
                ll={"longitude":pelleextent_wgs84.XMin, "latitude":pelleextent_wgs84.YMin},
                ur={"longitude":pelleextent_wgs84.XMax, "latitude":pelleextent_wgs84.YMax},
                start_date='2000-01-01', end_date='2018-12-31')
"""
landsat_raw = getfilelist(dir=landsatdir, repattern='LC08_L1[GTP]{2}.*B1[.]TIF$') #For explanation of repatter, check https://docs.python.org/3/howto/regex.html
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
                                           number_of_bands=1,
                                           mosaic_method='MAXIMUM')

    #Re-project to UTM 20S
    landsatproj[year] = '{}_UTM20S.tif'.format(os.path.splitext(outras)[0]) #os.path.splitext() method in Python is used to split the path name into a pair root and ext. Here, ext stands for extension
    if not os.path.exists(landsatproj[year]):
        print('Projecting landsat images from {}...'.format(year))
        arcpy.ProjectRaster_management(in_raster=outras,
                                       out_raster= landsatproj[year],
                                       out_coor_system = csref,
                                       resampling_type='NEAREST')


### Project and rasterize polygon of Pelegrinni department (matching grid layout of Landsat) ### 
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
arcpy.env.snapRaster = pelleras #Reference raster for snapping throughout analysis

#------------------------------------------------------------------------------
### Download, project, and clip Hansen et al. forest loss ###
GFCoutclip = {}
for product in ['datamask', 'gain', 'lossyear', 'treecover2000']:

    #Download
    GFCURL = "https://storage.googleapis.com/earthenginepartners-hansen/" \
             "GFC-2018-v1.6/Hansen_GFC-2018-v1.6_{}_20S_070W.tif".format(product)

    GFCout = os.path.join(forestdir, 'Hansen_GFC_v16_{}.tif'.format(product))

    if not os.path.exists(GFCout):
        print "downloading {}".format(GFCURL)
        with open(GFCout, "wb") as local_file:
            local_file.write(requests.get(GFCURL, stream=True).content)
    else:
        print('{} already exists...'.format(GFCout))

    GFCoutproj = os.path.join(forestoutdir, 'Hansen_GFC_v16_{}proj'.format(product))

    #Project to UTM 20S
    if not arcpy.Exists(GFCoutproj):
        print "projecting {}".format(GFCout)
        arcpy.ProjectRaster_management(in_raster=GFCout,
                                       out_raster= GFCoutproj,
                                       out_coor_system = csref,
                                       resampling_type='NEAREST',
                                       cell_size=arcpy.Describe(refraster).meanCellHeight)

    #Clip to study area
    GFCoutclip[product] = os.path.join(forestoutdir, 'Hansen_GFC_v16_{}clip'.format(product))
    if not arcpy.Exists(GFCoutclip[product]):
        print "clipping {}".format(GFCoutclip[product])
        ExtractByMask(in_raster=GFCoutproj,
                      in_mask_data=pelleras).save(GFCoutclip[product])

#------------------------------------------------------------------------------
# Creating resource weighting rasters for each livelihood based on deforestation.

# Compute yealry forest cover
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

#------------------------------------------------------------------------------
### Rasterize barrier lines ###
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

#------------------------------------------------------------------------------
#Assign weighting to barrier rasters
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

#------------------------------------------------------------------------------
#Prepare data for generation of cost distance rasters

#Create gdbs for cost raster
costtab_outgdb = {}
access_outgdb = {}
for llhood in livelihoods:
    #Output cost distance
    for year in analysis_years:
        costtab_outgdb[llhood+year] = os.path.join(resdir, 'Analysis_Chp1_W1',
                                                   'W1_3030', 'Cost_distance_W1_3030',
                                                   'Cost_distance_{}_w1'.format(llhood),
                                                   'CD_{0}_{1}_w1.gdb'.format(llhood, year))
        pathcheckcreate(costtab_outgdb[llhood+year])

    # Output access gdb
    access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                         'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
    pathcheckcreate(access_outgdb[llhood])
del llhood

#LOOP: for each livelihood, assign pixels to groups so that can run non-overlapping cost distance within livelihood-specific buffer
bufferad = {}
for llhood in livelihoods:
    #llhood = u'Bovine_livestock'
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
    arcpy.env.tileSize = "{0} {0}".format(int(((bufferad[llhood]-bufferad[llhood]%arcpy.Describe(refraster).meanCellWidth)*2 + arcpy.Describe(refraster).meanCellWidth*5)/arcpy.Describe(refraster).meanCellWidth))

    if not arcpy.Exists(fishras):
        arcpy.env.extent = fishout
        print('Converting fishnet to raster for {}...'.format(llhood))
        arcpy.PolygonToRaster_conversion (in_features = fishout,
                                          value_field = arcpy.Describe(fishout).OIDFieldName,
                                          out_rasterdataset = fishras,
                                          cell_assignment = 'CELL_CENTER',
                                          cellsize = refraster)

    #Convert study area raster to point (one point for each cell)
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

########################################################################################################################
########################################################################################################################

"""
Each job on Béluga should have a duration of at least one hour (five minutes for test jobs) 
and a user cannot have more than 1000 jobs, running and queued, at any given moment. 
The maximum duration for a job on Béluga is 7 days (168 hours).
A user can access 40 cores at a time, with each core 2 x Intel Gold 6148 Skylake @ 2.4 GHz, in total 1 x SSD 480G

To do:
- In calcgroup function:
    - Join zonal statistics results into table
    - Benchmark intersection of buffers with pellegrini department to speed up analysis + others

- Test new set up in parallel
- Zip it all up
"""
# LOOP: for each livelihood, for each group, create separate folder-points-buffers to run cost-distance in parallel on HPC

### ------ Add index to points, get list of groups, and run analysis on 10 groups to check speed ----- ####
fishgroups = defaultdict(set)
testdir = os.path.join(inputdir_HPC, 'testdir')
pathcheckcreate(testdir)
testgdb = os.path.join(testdir, 'test.gdb')
pathcheckcreate(testgdb)
grp_process_time = defaultdict(float)

for llhood in livelihoods:
    print('Assessing access calculationg run time for {}...'.format(llhood))
    if 'group{}'.format(llhood) not in [i.name for i in arcpy.Describe(pellepoints).indexes]:
        print('Adding index to pellepoints...')
        arcpy.AddIndex_management(pellepoints, fields='group{}'.format(llhood),
                                  index_name='group{}'.format(llhood))  # Add index to speed up querying

    if llhood not in fishgroups:
        print('Getting groups...')
        fishgroups[llhood] = {row[0] for row in arcpy.da.SearchCursor(pellepoints, 'group{}'.format(llhood))}

        # Output points for 10 groups for each livelihood
    for group in list(fishgroups[llhood])[0:5]:
        print(group)

        outpoints = os.path.join(testdir, 'testpoints_{0}_{1}_{2}.shp').format(llhood, int(bufferad[llhood]), group)
        if not arcpy.Exists(outpoints):
            # Subset points based on group (so that their buffers don't overlap) and only keep points that overlap study area
            arcpy.MakeFeatureLayer_management(in_features=pellepoints, out_layer='subpoints{}'.format(group),
                                              where_clause='{0} = {1}'.format('group{}'.format(llhood), group))
            arcpy.CopyFeatures_management('subpoints{}'.format(group), outpoints)

        # Test time that each group takes for each livelihood

        inbw = {yr: barrierweight_outras['{0}{1}'.format(llhood, yr)] for yr in analysis_years}

        tic = time.time()
        # Get subpoints
        accesscalc3(inllhood=llhood,
                    ingroup=group,
                    inpoints=outpoints,
                    inbuffer_radius=bufferad[llhood],
                    inlimits=pelleras,
                    inyears=analysis_years,
                    inbarrierweight_outras=inbw,
                    inforestyearly=forestyearly,
                    costtab_outgdb=testgdb)
        toc = time.time()
        print(toc - tic)
        grp_process_time[llhood] = grp_process_time[llhood] + (toc - tic) / 5

### ------ Compute number of chunks to divide each livelihood in to process each chunk with equal time ------###
numcores = 40  # Number of chunks to divide processing into
maxdays = 1 #Max number of days that processes can be run at a time

#Filter groups to process based on those already processed
processeddict = defaultdict(list)
groupstoprocess = defaultdict(list)
for llhood in livelihoods:
    print(llhood)
    #Add group for every year
    for year in analysis_years:
        arcpy.env.workspace = costtab_outgdb[llhood + year]
        for tab in arcpy.ListTables():
            processeddict[llhood].append(
                int(re.search('(?<=w1[_])[0-9]*'.format(llhood, year), tab).group()))
    #Only keep groups that have not been processed for all years in analysis_years
    groupstoprocess[llhood] = [g for g in fishgroups[llhood] if not \
        Counter(processeddict[llhood])[g] == len(analysis_years)]

#Assess amount of time reguired and number of chunks
totaltime = sum([grp_process_time[llhood] * len(groupstoprocess[llhood]) for llhood in livelihoods])
print('Total processing times among {0} cores: {1} days...'.format(
    numcores, totaltime/float(3600*24*numcores))) #Total time if process is divided into numcores chunks at same speed
numchunks = math.ceil(totaltime / float(3600 * 24 * maxdays))
print('Total number of chunks for each to be processed within {0} days among {1} cores: {2}...'.format(
    maxdays, numcores, numchunks))

### ------ Assign groups to chunks ------###
llhood_chunks = {}
formatdir = os.path.join(inputdir_HPC, 'Beluga_input{}'.format(time.strftime("%Y%m%d")))
formatdir_data = os.path.join(formatdir, 'data')
formatdir_results = os.path.join(formatdir, 'results')
pathcheckcreate(formatdir_data, verbose=True)
pathcheckcreate(formatdir_results, verbose=True)

for llhood in livelihoods:
    print(llhood)
    llhood_chunks[llhood] = math.ceil(numchunks * grp_process_time[llhood] * len(groupstoprocess[llhood]) / totaltime)
    print('    Number of daily-core chunks to divide {0} groups into: {1}...'.format(
        llhood, llhood_chunks[llhood]))

    groupchunklist = groupindexing(grouplist=list(groupstoprocess[llhood]), chunknum=llhood_chunks[llhood])

    #Output points and ancillary data to chunk-specific gdb
    for chunk in range(0, len(groupchunklist)):
        print(chunk)
        outchunkgdb = os.path.join(formatdir_data, '{0}{1}_{2}.gdb'.format(llhood, int(bufferad[llhood]), chunk))
        pathcheckcreate(outchunkgdb, verbose=True)
        outchunkpoints= os.path.join(outchunkgdb, 'subpoints{0}{1}_{2}'.format(llhood, int(bufferad[llhood]), chunk))

        print('Copying points...')
        arcpy.CopyFeatures_management(
            arcpy.MakeFeatureLayer_management(in_features=pellepoints, out_layer='pointslyr',
                                              where_clause='group{0} IN {1}'.format(llhood,tuple(groupchunklist[chunk]))),
            outchunkpoints)

        print('Copying ancillary data...')
        for yr in analysis_years:
            #Copy barrier raster
            arcpy.CopyRaster_management(barrierweight_outras[llhood + yr],
                                        os.path.join(outchunkgdb, os.path.split(barrierweight_outras[llhood + yr])[1]))
            #Copy forest cover
            arcpy.CopyRaster_management(forestyearly[yr],
                                        os.path.join(outchunkgdb, os.path.split(forestyearly[yr])[1]))