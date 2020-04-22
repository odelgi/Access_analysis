# Import system modules
import arcpy
from arcpy.sa import *
import numpy as np
import os #Module for operating system operations (e.g. os.path.exists)
import re #Module for string pattern matching
from collections import defaultdict
from collections import Counter
import time
import pandas as pd
import math
from accesscalc_parallel import *

### Set environment settings ###
#https://pro.arcgis.com/en/pro-app/arcpy/classes/env.htm
arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

# File extension for the output rasters
# output_ext = ".img"

analysis_years = ['2000', '2010', '2018'] #Years for analysis
numyears = len(analysis_years)

#### Directory structure ###
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

inputdir_HPC = os.path.join(resdir, 'Analysis_Chp1_W1', 'inputdata_HPC')
pathcheckcreate(inputdir_HPC)

outstats = os.path.join(resdir, 'Analysis_Chp1_W1/W1_3030/Cost_distance_W1_3030')

#Repreate paths
weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')
livelihoods = weightingpd['Livelihood'].tolist()
livelihoods.remove('Combined_livelihood')
llhoodbuffer_outdir = os.path.join(resdir, 'Analysis_Chp1_W1',
                                   'Buffers_W1.gdb')  # Output buffer path for each livelihood
pellepoints = os.path.join(llhoodbuffer_outdir, 'pellefishpoints')

#Recreate paths
barrierweight_outdir = os.path.join(resdir,'Analysis_Chp1_W1','W1_3030','Barrier_weighting_1_3030')
forestoutdir = os.path.join(resdir, 'Base_layers_Pellegrini/Forest_Hansen.gdb')
barrierweight_outras = {}; forestyearly = {}; bufferad = {}; costtab_outgdb = {}
for llhood in livelihoods:
    print(llhood)
    bufferad[llhood] = float(weightingpd.loc[weightingpd['Livelihood'] == llhood, 'Buffer_max_rad'])
    outllhood_gdb = os.path.join(resdir,'Analysis_Chp1_W1','W1_3030','Barrier_weighting_1_3030','{}_bw1_3030.gdb'.format(llhood))
    for year in analysis_years:
        barrierweight_outras[llhood+year] = os.path.join(outllhood_gdb, '{0}_bw1_{1}'.format(llhood, year))
        forestyearly[year] = os.path.join(forestoutdir, 'Hansen_GFC_v16_treecover{}'.format(year))
        costtab_outgdb[llhood+year] = os.path.join(resdir, 'Analysis_Chp1_W1',
                                                   'W1_3030', 'Cost_distance_W1_3030',
                                                   'Cost_distance_{}_w1'.format(llhood),
                                                   'CD_{0}_{1}_w1.gdb'.format(llhood, year))

# LOOP: for each livelihood, for each group, create separate folder-points-buffers to run cost-distance in parallel on HPC

### ------ Add index to points, get list of groups, and run analysis on 20 groups to check speed ----- ####
fishgroups = defaultdict(set)
testdir = os.path.join(inputdir_HPC, 'testdir')
pathcheckcreate(testdir)
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
    for group in list(fishgroups[llhood])[0:9]:
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
        accesscalc(inllhood=llhood,
                    ingroup=group,
                    inpoints=outpoints,
                    inbuffer_radius=bufferad[llhood],
                    inyears=analysis_years,
                    inbarrierweight_outras=inbw,
                    inforestyearly=forestyearly,
                    costtab_outdir=testdir)
        toc = time.time()
        print(toc - tic)
        grp_process_time[llhood] = grp_process_time[llhood] + (toc - tic) / 5

### ------ Compute number of chunks to divide each livelihood in to process each chunk with equal time ------###
numcores = 4  # Number of chunks to divide processing into
maxdays = 4 #Max number of days that processes can be run at a time

#Get all processed tables
tablist = getfilelist(dir=outstats, repattern=".*[.]dbf$", gdbf = False, nongdbf = True)
tablist.extend(getfilelist(dir=outstats, gdbf = True, nongdbf = False))

tables_pd = pd.concat([pd.Series(tablist),
                       pd.Series(tablist).apply(lambda x: os.path.splitext(os.path.split(x)[1])[0]).
                      str.split('_', expand=True)],
                      axis=1)
tables_pd.columns = ['path', 'dataset', 'llhood1', 'llhood2', 'year', 'weighting', 'group']
tables_pd['llhood'] = tables_pd['llhood1'] + '_' + tables_pd['llhood2']
tables_pd = tables_pd.drop(labels=['llhood1', 'llhood2'], axis=1)

processed_pd = tables_pd.groupby(['llhood', 'group']).filter(lambda x: x['year'].nunique() == 3).\
    drop_duplicates(subset=['llhood', 'group'])

#Filter groups to process based on those already processed
groupstoprocess = defaultdict(list)
for llhood in livelihoods:
    subprocessed_set = set(processed_pd.loc[processed_pd['llhood'] == llhood]['group'])
    #Only keep groups that have not been processed for all years in analysis_years
    groupstoprocess[llhood] = {g for g in fishgroups[llhood] if g not in subprocessed_set}

#Assess amount of time reguired and number of chunks
totaltime = sum([grp_process_time[llhood] * len(groupstoprocess[llhood]) for llhood in livelihoods])
print('Total processing times among {0} cores: {1} days...'.format(
    numcores, totaltime/float(3600*24*numcores))) #Total time if process is divided into numcores chunks at same speed
numchunks = math.ceil(totaltime / float(3600 * 24 * maxdays))
print('Total number of chunks for each to be processed within {0} days among {1} cores: {2}...'.format(
    maxdays, numcores, numchunks))

### ------ Assign groups to chunks ------###
llhood_chunks = {}
formatdir = os.path.join(inputdir_HPC, 'Black_input{}'.format(time.strftime("%Y%m%d")))
formatdir_data = os.path.join(formatdir, 'data')
formatdir_results = os.path.join(formatdir, 'results')
pathcheckcreate(formatdir_data, verbose=True)
pathcheckcreate(formatdir_results, verbose=True)

for llhood in livelihoods:
    print(llhood)
    llhood_chunks[llhood] = math.ceil(numchunks * grp_process_time[llhood] * len(groupstoprocess[llhood]) / totaltime)
    print('    Number of chunks to divide {0} groups into: {1}...'.format(
        llhood, llhood_chunks[llhood]))

    groupchunklist = groupindexing(grouplist=list(groupstoprocess[llhood]), chunknum=llhood_chunks[llhood])

    #Output points and ancillary data to chunk-specific gdb
    for chunk in range(0, len(groupchunklist)):
        print(chunk)
        outchunkgdb = os.path.join(formatdir_data, '{0}{1}_{2}.gdb'.format(llhood, int(bufferad[llhood]), chunk))
        if (arcpy.Exists(outchunkgdb)):
            arcpy.Delete_management(outchunkgdb)
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
            if llhood == 'Charcoal_production':
                arcpy.CopyRaster_management(forestyearly[yr],
                                            os.path.join(outchunkgdb, os.path.split(forestyearly[yr])[1]))