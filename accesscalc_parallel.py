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

arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

def accesscalc_grouped(infishgroup, inpoints, inllhood, inbuffer_radius, inyears, inbarrierweight_outras,
                       inforestyearly, incosttab_outgdb):
    print('Processing group # {}...'.format(infishgroup))

    #tmpdir = os.path.join(os.path.dirname(incosttab_outgdb[inllhood+'2000']), 'tmp_{}'.format(str(infishgroup)))

    #try:
        #os.mkdir(tmpdir)
    arcpy.env.scratchWorkspace = r"in_memory"  # Not having a separate scratch workspace can lead to bad locking issues

    tic = time.time()

    # Buffer subsetted points based on livelihood-specific buffer distance
    arcpy.Buffer_analysis(in_features='subpoints{}'.format(group),
                          out_feature_class=os.path.join(groupdir, 'buffers.shp'),
                          buffer_distance_or_field=bufferad,
                          dissolve_option='NONE',
                          method='PLANAR')

    # Iterate through years of analysis
    for year in inyears:
        print(year)
        # Compute cost distance and access
        arcpy.env.mask = 'in_memory/subbuffers{}'.format(infishgroup)
        if inllhood == 'Charcoal_production':
            # forest resource weighting*(1/(1+cost))
            accessras = Int(100 * Raster(inforestyearly[year]) * \
                            (1 / (1 + CostDistance(in_source_data='subpoints{}'.format(infishgroup),
                                                   in_cost_raster=inbarrierweight_outras[inllhood + year]))))
        else:
            # (1/(1+cost))
            accessras = Int(100 * (1 / (1 + CostDistance(in_source_data='subpoints{}'.format(infishgroup),
                                                         in_cost_raster=inbarrierweight_outras[inllhood + year]))))

        # Zonal statistics based on buffer (using pointid, the unique ID of each point for that livelihood)
        # Compute mean access within livelihood-specific buffer and writes it out to table
        ZonalStatisticsAsTable(in_zone_data='in_memory/subbuffers{}'.format(infishgroup),
                               zone_field='pointid',
                               in_value_raster=accessras,
                               out_table=os.path.join(incosttab_outgdb[inllhood + year],
                                                      'CD_{0}_{1}_w1_{2}'.format(inllhood, year, infishgroup)),
                               ignore_nodata='DATA',
                               statistics_type='MEAN')
    toc = time.time()
    print('Took {} s...'.format(toc - tic))

    #except Exception:
    #    traceback.print_exc()
        #arcpy.Delete_management(tmpdir)

    #if arcpy.Exists(tmpdir):
    #   arcpy.Delete_management(tmpdir)

#accesscalc_grouped(infishgroup=90000, inpoints=pellepoints, inllhood=llhood, inbuffer_radius=bufferad,
#                                             inyears=analysis_years, inbarrierweight_outras=barrierweight_outras,
#                                             inforestyearly=forestyearly, incosttab_outgdb=costtab_outgdb)

if __name__ == '__main__':
    #### Directory structure ###
    rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
    datadir = os.path.join(rootdir, 'data')
    resdir = os.path.join(rootdir, 'results')

    # Functions

    def getfilelist(dir, repattern):
        """Function to iteratively go through all subdirectories inside 'dir' path
        and retrieve path for each file that matches "repattern"""
        return [os.path.join(dirpath, file)
                for (dirpath, dirnames, filenames) in os.walk(dir)
                for file in filenames if re.search(repattern, file)]

    def pathcheckcreate(path):
        """"Function that takes a path as input and:
          1. Checks which directories and .gdb exist in the path
          2. Creates the ones that don't exist"""

        dirtocreate = []
        # Loop upstream through path to check which directories exist, adding those that don't exist to dirtocreate list
        while not os.path.exists(os.path.join(path)):
            dirtocreate.append(os.path.split(path)[1])
            path = os.path.split(path)[0]

        dirtocreate.reverse()

        # After reversing list, iterate through directories to create starting with the most upstream one
        for dir in dirtocreate:
            # If gdb doesn't exist yet, use arcpy method to create it and then stop the loop to prevent from trying to create anything inside it
            if os.path.splitext(dir)[1] == '.gdb':
                print('Create {}...'.format(dir))
                arcpy.CreateFileGDB_management(out_folder_path=path,
                                               out_name=dir)
                break

            # Otherwise, if it is a directory name (no extension), make a new directory
            elif os.path.splitext(dir)[1] == '':
                print('Create {}...'.format(dir))
                path = os.path.join(path, dir)
                os.mkdir(path)

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
