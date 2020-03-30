import arcpy
from arcpy.sa import *
import time
import os
import re
import pandas as pd
import multiprocessing
from functools import partial
import traceback
import cProfile

from internal_functions import *

arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")


#accesscalc_grouped(infishgroup=90000, inpoints=pellepoints, inllhood=llhood, inbuffer_radius=bufferad,
#                                             inyears=analysis_years, inbarrierweight_outras=barrierweight_outras,
#                                             inforestyearly=forestyearly, incosttab_outgdb=costtab_outgdb)

if __name__ == '__main__':
    #### Directory structure ###
    rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
    datadir = os.path.join(rootdir, 'data')
    resdir = os.path.join(rootdir, 'results')

    # Functions

    # Parameter
    analysis_years = ['2000', '2010', '2018']
    refraster = os.path.join(resdir,
                             "Base_layers_Pellegrini",
                             "Basemaps_UTM20S.gdb",
                             "Pellegri_department_UTM20S")  # Reference raster for snapping throughout analysis

    # Define paths for input variables
    fishpoints = os.path.join(resdir, 'Analysis_Chp1_W1', 'Buffers_W1.gdb', 'fishnetpoints')
    weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
    llhoodbuffer_outdir = os.path.join(resdir, 'Analysis_Chp1_W1',
                                       'Buffers_W1.gdb')  # Output buffer path for each livelihood
    pellepoints = os.path.join(llhoodbuffer_outdir, 'pellefishpoints')

    weightingmeta = pd.read_excel(weighting_table,
                                  sheetname='Metadata')  # Read weighting excel table, only Metadata sheet
    weightingcodes = weightingmeta.loc[weightingmeta['Variable name'].str.contains(
        "BW"),  # Subset metadata sheet to only keep rows for which 'Variable name' column string contains 'BW'
                                       ['Variable name',
                                        'Code']]  # For those rows, only keep Variable name and Code columns
    weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

    livelihoods = weightingpd['Livelihood'].tolist()
    livelihoods.remove('Combined_livelihood')

    forestoutdir = os.path.join(resdir, 'Base_layers_Pellegrini/Forest_Hansen.gdb')

    barrierweight_outras = {}
    forestyearly = {}
    costtab_outgdb = {}
    access_outgdb = {}
    access_outras = {}
    for llhood in livelihoods:
        outllhood_gdb = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Barrier_weighting_1_3030',
                                     '{}_bw1_3030.gdb'.format(llhood))
        pathcheckcreate(outllhood_gdb)

        access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                             'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
        pathcheckcreate(access_outgdb[llhood])

        for year in analysis_years:
            barrierweight_outras[llhood + year] = os.path.join(outllhood_gdb, '{0}_bw1_{1}'.format(llhood, year))
            forestyearly[year] = os.path.join(forestoutdir, 'Hansen_GFC_v16_treecover{}'.format(year))
            costtab_outgdb[llhood + year] = os.path.join(resdir, 'Analysis_Chp1_W1',
                                                         'W1_3030', 'Cost_distance_W1_3030',
                                                         'Cost_distance_{}_w1'.format(llhood),
                                                         'CD_{0}_{1}_w1.gdb'.format(llhood, year))
            pathcheckcreate(costtab_outgdb[llhood + year])

            access_outras[llhood + year] = os.path.join(access_outgdb[llhood],
                                                        'accessras_W1_{0}{1}'.format(llhood, year))
    del llhood;
    del year

    #
    #for llhood in livelihoods:
    llhood = 'Bovine_livestock'
    if not any([arcpy.Exists(layer) for layer in [access_outras[llhood + y] for y in analysis_years]]):
        print('Processing {}...'.format(llhood))
        bufferad = float(weightingpd.loc[weightingpd[
                                             'Livelihood'] == llhood, 'Buffer_max_rad'])  # Get livelihood_specific buffer radius from table
        print('Adding index to pellepoints...')
        if 'group{}'.format(llhood) not in [i.name for i in arcpy.Describe(pellepoints).indexes]:
            arcpy.AddIndex_management(pellepoints, fields='group{}'.format(llhood),
                                      index_name='group{}'.format(llhood))  # Add index to speed up querying

        print('Getting groups...')
        fishgroups = {row[0] for row in arcpy.da.SearchCursor(pellepoints, 'group{}'.format(llhood))}

        ############# PARALLEL START ###################################################################################
        print('Launch parallel processing')
        tic = time.time()
        p = multiprocessing.Pool(int(multiprocessing.cpu_count()))
        accesscalc_grouped_partial = partial(accesscalc_grouped,
                                             inpoints=pellepoints, inllhood=llhood, inbuffer_radius=bufferad,
                                             inyears=analysis_years, inbarrierweight_outras=barrierweight_outras,
                                             inforestyearly=forestyearly, incosttab_outgdb=costtab_outgdb)
        p.map(accesscalc_grouped_partial, fishgroups)
        p.close()
        print(time.time() - tic)
        ############# PARALLEL END #####################################################################################

        # Process: Merge tables and join to points (within loop)
        for year in analysis_years:
            print(year)
            accessdict = {}
            for dirpath, dirnames, filenames in arcpy.da.Walk(costtab_outgdb[llhood + year],
                                                              datatype="Table"):  # Retrieve the names of all tables in livelihood-specific access geodatabase
                for tab in [os.path.join(dirpath, f) for f in filenames]:
                    print(tab)
                    for row in arcpy.da.SearchCursor(tab, ['pointid', 'MEAN']):
                        accessdict[row[0]] = row[1]

            outfield = 'access{0}{1}'.format(llhood, year)
            if not outfield in [f.name for f in arcpy.ListFields(
                    pellepoints)]:  # Make sure that livelihood specific access field doesn't already exist for that year in pellepoints
                print('Create {} field'.format(outfield))
                arcpy.AddField_management(in_table=pellepoints, field_name=outfield, field_type='FLOAT')

                print('Write values out to point dataset...')
                with arcpy.da.UpdateCursor(pellepoints, ['pointid', outfield]) as cursor:
                    x = 0
                    for row in cursor:
                        if x % 100000 == 0:
                            print(x)
                        try:
                            row[1] = accessdict[row[0]]
                        except:
                            print('pointid {} was not found in dictionary'.format(row[0]))
                        cursor.updateRow(row)
                        x += 1

            # Convert points back to raster
            arcpy.PointToRaster_conversion(in_features=pellepoints, value_field=outfield, cellsize=refraster,
                                           out_rasterdataset=access_outras)
