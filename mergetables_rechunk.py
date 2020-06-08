# Import system modules
import arcpy
from arcpy.sa import *
import numpy as np
import os #Module for operating system operations (e.g. os.path.exists)
from shutil import copyfile
import re #Module for string pattern matching
from collections import defaultdict
from collections import Counter
import time
import pandas as pd
import math
from accesscalc_parallel import *

#Define function
def groupindexing(grouplist, chunknum):
    # Set of groups
    # Chunk size to divide groups into
    x = np.array(grouplist)
    bin_step = max(grouplist) / chunknum
    bin_edges = np.arange(min(i for i in grouplist if i is not None),
                          max(i for i in grouplist if i is not None) + bin_step,
                          bin_step)
    bin_number = bin_edges.size - 1
    cond = np.zeros((x.size, bin_number), dtype=bool)
    for i in range(bin_number):
        cond[:, i] = np.logical_and(bin_edges[i] <= x,
                                    x < bin_edges[i + 1])
    return [list(x[cond[:, i]]) for i in range(bin_number)]

def tabmerge_dict(tablist):
    outdict = {}
    for tab in tablist:
        print(tab)
        for row in arcpy.da.SearchCursor(tab, ['Value', 'MEAN']):
            outdict[row[0]] = row[1]
    return(outdict)

### Set environment settings ###
#https://pro.arcgis.com/en/pro-app/arcpy/classes/env.htm
arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

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
basemapgdb = os.path.join(resdir, "Base_layers_Pellegrini", "Basemaps_UTM20S.gdb")
pelleras = os.path.join(basemapgdb, "Pellegri_department_UTM20S")

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
### ------- Get all processed tables -----------------------------------------------------------------------------------
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

### ------ Create a raster of access for each llhood and year by aggregating all tables (yielding an access value for each pixel-point) ---------
access_outgdb = {}
refraster = pelleras
fishgroups = defaultdict(set)

# Iterate over each livelihood
for llhood in tables_pd['llhood'].unique():
    access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                         'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
    pathcheckcreate(access_outgdb[llhood])
    # Iterate over each year
    for year in tables_pd['year'].unique():
        # Path of output access raster
        access_outras = os.path.join(access_outgdb[llhood], 'accessras_W1_{0}{1}'.format(llhood, year))

        # Perform analysis only if output raster doesn't exist
        if not arcpy.Exists(access_outras):
            print("Processing {}...".format(access_outras))

            # Aggregate values across all pixels-points for that livelihood-year
            print('Aggregating zonal statistics tables...')
            merged_dict = tabmerge_dict(tables_pd.loc[(tables_pd['llhood'] == llhood) &
                                                      (tables_pd['year'] == year), 'path'])

            if len(merged_dict) > 0:
                # Join all statistics tables of access to pellepoints (a point for each 30x30 m pixel in Pellegrini department)
                print('Joining tables to points...')
                accessfield = 'access{0}{1}'.format(llhood, year)
                if not accessfield in [f.name for f in arcpy.ListFields(pellepoints)]:
                    print('Create {} field'.format(accessfield))
                    arcpy.AddField_management(in_table=pellepoints, field_name=accessfield, field_type='FLOAT')

                with arcpy.da.UpdateCursor(pellepoints, ['pointid', accessfield, 'group{}'.format(llhood)]) as cursor:
                    x = 0
                    for row in cursor:
                        if x % 100000 == 0:
                            print(x)
                        try:
                            row[1] = merged_dict[row[0]]
                            fishgroups[llhood].append(row[2])
                        except:
                            print('pointid {} was not found in dictionary'.format(row[0]))
                        cursor.updateRow(row)
                        x += 1

                # Convert points back to raster
                if len(merged_dict) == x:
                    print('Converting points to raster...')
                    arcpy.PointToRaster_conversion(in_features=pellepoints,
                                                   value_field=accessfield,
                                                   cellsize=refraster,
                                                   out_rasterdataset=access_outras)
            else:
                print('No zonal statistics available for that livelihood for that year...')

        else:
            print('{} already exists...'.format(access_outras))

### ------ Run analysis on 10 groups to check speed ----- ####
testdir = os.path.join(inputdir_HPC, 'testdir')
pathcheckcreate(testdir)
grp_process_time = defaultdict(float)
for llhood in livelihoods:
    print('Assessing access calculation run time for {}...'.format(llhood))

    if 'group{}'.format(llhood) not in [i.name for i in arcpy.Describe(pellepoints).indexes]:
        print('Adding index to pellepoints...')
        arcpy.AddIndex_management(pellepoints, fields='group{}'.format(llhood),
                                  index_name='group{}'.format(llhood))  # Add index to speed up querying

    # Output points for 10 groups for each livelihood
    for group in list(fishgroups[llhood])[0:10]:
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
        grp_process_time[llhood] = grp_process_time[llhood] + (toc - tic) / 10.0

### ------ Compute number of chunks to divide each livelihood in to process each chunk with equal time ------###
numcores = 14  # Number of chunks to divide processing into
maxdays = 1 #Max number of days that processes can be run at a time

#Filter groups to process based on those already processed
groupstoprocess = defaultdict(list)
for llhood in livelihoods:
    subprocessed_set = set(processed_pd.loc[processed_pd['llhood'] == llhood]['group'].astype('int'))
    #Only keep groups that have not been processed for all years in analysis_years
    groupstoprocess[llhood] = {g for g in fishgroups[llhood] if g not in subprocessed_set}

#Assess amount of time reguired and number of chunks
totaltime = sum([grp_process_time[llhood] * len(groupstoprocess[llhood]) for llhood in livelihoods])
print('Total processing times among {0} cores: {1} days...'.format(
    numcores, totaltime/float(3600.0*24*numcores))) #Total time if process is divided into numcores chunks at same speed
numchunks = math.ceil(totaltime / float(3600.0 * 24 * maxdays))
print('Total number of chunks for each to be processed within {0} days among {1} cores: {2}...'.format(
    maxdays, numcores, numchunks))

### ------ Assign groups to chunks ------###
llhood_chunks = {}
formatdir = os.path.join(inputdir_HPC, 'White_input{}'.format(time.strftime("%Y%m%d")))
formatdir_data = os.path.join(formatdir, 'data')
formatdir_results = os.path.join(formatdir, 'results')
formatdir_src = os.path.join(formatdir, 'src')
pathcheckcreate(formatdir_data, verbose=True)
pathcheckcreate(formatdir_results, verbose=True)
pathcheckcreate(formatdir_src, verbose=True)

#Copy processing file to directory
in_processingscript = os.path.join(rootdir, 'src', 'Chap1_Analysis1', 'accesscalc_parallel.py')
out_processingscript = os.path.join(formatdir_src, 'accesscalc_parallel.py')
copyfile(in_processingscript, out_processingscript)

for llhood in livelihoods:
    if len(groupstoprocess[llhood]) > 0:
        print(llhood)
        llhood_chunks[llhood] = math.ceil(numchunks * grp_process_time[llhood] * len(groupstoprocess[llhood]) / totaltime)
        print('    Number of chunks to divide {0} groups into: {1}...'.format(
            llhood, llhood_chunks[llhood]))

        #groupchunklist = groupindexing(grouplist=list(groupstoprocess[llhood]), chunknum=llhood_chunks[llhood])
        interval = int(math.ceil(len(groupstoprocess[llhood])/ float(numchunks)))
        groupchunklist = [list(groupstoprocess[llhood])[i:(i + interval)] for i
                   in range(0, len(groupstoprocess[llhood]), interval)]

        #Output points and ancillary data to chunk-specific gdb
        for chunk in range(0, len(groupchunklist)):
            print(chunk)
            outchunkgdb = os.path.join(formatdir_data, '{0}{1}_{2}.gdb'.format(llhood, int(bufferad[llhood]), chunk))
            if not (arcpy.Exists(outchunkgdb)):
                pathcheckcreate(outchunkgdb, verbose=True)
                outchunkpoints= os.path.join(outchunkgdb, 'subpoints{0}{1}_{2}'.format(llhood, int(bufferad[llhood]), chunk))

                print('Copying points...')
                if len(groupchunklist[chunk])>0:
                    print(len(groupchunklist[chunk]))
                    arcpy.CopyFeatures_management(
                        arcpy.MakeFeatureLayer_management(
                            in_features=pellepoints, out_layer='pointslyr',
                            where_clause='group{0} IN {1}'.format(llhood,tuple(groupchunklist[chunk]))), #[i for i in groupchunklist[chunk] if i is not None]
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
            else:
                print('{} already exists...'.format(outchunkgdb))
    else:
        print('All groups for {} have already been processed...'.format(llhood))