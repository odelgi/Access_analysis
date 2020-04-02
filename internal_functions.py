import os
import sys
import traceback
import re
import arcpy
from arcpy.sa import *
import time
import numpy as np
from concurrent.futures import TimeoutError

# Functions
def getfilelist(dir, repattern):
    """Function to iteratively go through all subdirectories inside 'dir' path
    and retrieve path for each file that matches "repattern"""
    return [os.path.join(dirpath, file)
            for (dirpath, dirnames, filenames) in os.walk(dir)
            for file in filenames if re.search(repattern, file)]

def pathcheckcreate(path, verbose=True):
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
            if verbose:
                print('Create {}...'.format(dir))
            arcpy.CreateFileGDB_management(out_folder_path=path,
                                           out_name=dir)
            break

        # Otherwise, if it is a directory name (no extension), make a new directory
        elif os.path.splitext(dir)[1] == '':
            if verbose:
                print('Create {}...'.format(dir))
            path = os.path.join(path, dir)
            os.mkdir(path)

def groupindexing(grouplist, chunknum):
    # Set of groups
    # Chunk size to divide groups into
    x = np.array(grouplist)
    bin_step = max(grouplist) / chunknum
    bin_edges = np.arange(min(grouplist),
                          max(grouplist) + bin_step,
                          bin_step)
    bin_number = bin_edges.size - 1
    cond = np.zeros((x.size, bin_number), dtype=bool)
    for i in range(bin_number):
        cond[:, i] = np.logical_and(bin_edges[i] <= x,
                                    x < bin_edges[i + 1])
    return [list(x[cond[:, i]]) for i in range(bin_number)]

def accesscalc(inllhood, ingroup, inpoints, inbuffer_radius, inyears, inbarrierweight_outras,
               inforestyearly, costtab_outgdb):
    # Buffer subsetted points based on livelihood-specific buffer distance
    #print('Bufferring...')
    try:
        arcpy.Buffer_analysis(in_features=inpoints,
                              out_feature_class=r'in_memory/subbuffers{}'.format(ingroup),
                              buffer_distance_or_field=inbuffer_radius,
                              dissolve_option='NONE',
                              method='PLANAR')

        templateras = inforestyearly[inyears[0]]
        arcpy.env.SnapRaster = templateras
        buffras_memory = r'in_memory/subbufferas{}'.format(ingroup)
        #print('Conversion to raster...')
        arcpy.FeatureToRaster_conversion(in_features=r'in_memory/subbuffers{}'.format(ingroup),
                                         field='pointid',
                                         out_raster=buffras_memory,
                                         cell_size=templateras)

        # Iterate through years of analysis
        #print('Iterating through years...')
        for year in inyears:
            #print(year)
            # Compute cost distance and access
            arcpy.env.mask = buffras_memory
            if inllhood == 'Charcoal_production':
                # forest resource weighting*(1/(1+cost))
                accessras = Int(100 * Raster(inforestyearly[year]) * \
                                (1 / (1 + CostDistance(in_source_data=inpoints,
                                                       in_cost_raster=inbarrierweight_outras[year]))))
            else:
                # (1/(1+cost))
                accessras = Int(100 * (1 / (1 + CostDistance(in_source_data=inpoints,
                                                             in_cost_raster=inbarrierweight_outras[year]))))

            # Zonal statistics based on buffer (using pointid, the unique ID of each point for that livelihood)
            # Compute mean access within livelihood-specific buffer and writes it out to table
            #accessras.save(os.path.join(costtab_outgdb, 'test_{0}_{1}_w1_{2}'.format(inllhood, year, ingroup))
            ZonalStatisticsAsTable(in_zone_data=buffras_memory,
                                   zone_field='Value',
                                   in_value_raster=accessras,
                                   out_table=os.path.join(costtab_outgdb,
                                                          'CD_{0}_{1}_w1_{2}'.format(inllhood, year, ingroup)),
                                   ignore_nodata='DATA',
                                   statistics_type='MEAN')

        #Clean up
        try:
            del accessras
        except:
            pass
        arcpy.ClearEnvironment("mask")
        arcpy.Delete_management("in_memory")

    except:
        #Clean up memory, environment, and temporary files
        try:
            del accessras
        except:
            pass

        arcpy.ClearEnvironment("mask")
        arcpy.Delete_management("in_memory")

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Print Python error messages
        print("PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1]))
        print("ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n")

def accesscalc_chunkedir(inchunkgdb, inyears, outgdb):
    arcpy.env.workspace = inchunkgdb
    arcpy.env.scratchWorkspace = 'in_memory'
    rootname = os.path.splitext(os.path.split(inchunkgdb)[1])[0]
    llhood = re.search('[aA-zZ]*[_][aA-zZ]*', rootname).group()
    bufferrad = re.search('[0-9]*(?=_[0-9]*$)', rootname).group()
    inpoints = 'subpoints{}'.format(rootname)
    barrierweight_outras = {year: '{0}_bw1_{1}'.format(llhood, year) for year in inyears}
    forestyearly = {year: 'Hansen_GFC_v16_treecover{}'.format(year) for year in inyears}
    print(rootname)

    if 'group{}'.format(llhood) not in [i.name for i in arcpy.Describe(inpoints).indexes]:
        arcpy.AddIndex_management(inpoints, fields='group{}'.format(llhood),
                                  index_name='group{}'.format(llhood))  # Add index to speed up querying

    # Iterate through groups
    grouplist = {row[0] for row in arcpy.da.SearchCursor(inpoints, 'group{}'.format(llhood))}
    for group in grouplist:
        print(group)
        #tic = time.time()
        arcpy.MakeFeatureLayer_management(in_features=inpoints, out_layer='pointslyr_{}'.format(group),
                                          where_clause='{0} = {1}'.format('group{}'.format(llhood), group))
        accesscalc(inllhood=llhood,
                   ingroup=group,
                   inpoints='pointslyr_{}'.format(group),
                   inbuffer_radius=bufferrad,
                   inyears=inyears,
                   inbarrierweight_outras=barrierweight_outras,
                   inforestyearly=forestyearly,
                   costtab_outgdb=outgdb)

        arcpy.Delete_management('pointslyr_{}'.format(group))
        #print(time.time() - tic)
    arcpy.ClearEnvironment('workspace')

def task_done(future):
    try:
        future.result()  # blocks until results are ready
    except TimeoutError as error:
        print("Function took longer than %d seconds" % error.args[1])
    except Exception as error:
        print("Function raised %s" % error)
        print(error.traceback)  # traceback of the function

if __name__ == '__main__':
    ####### TROUBLESHOOTING AND PROFILING CODE #############
    import cProfile

    arcpy.env.overwriteOutput = 'True'
    arcpy.CheckOutExtension("Spatial")

    # Define directory structure
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0],
                             'results\Analysis_Chp1_W1\inputdata_HPC/Beluga_input20200331')  # To update for final run
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    arcpy.env.workspace = datadir
    ingdbs = arcpy.ListWorkspaces('*', workspace_type='FileGDB')
    outstats = os.path.join(resdir, 'outstats.gdb')
    pathcheckcreate(outstats)
    analysis_years = ['2000', '2010', '2018']

    # for gb in ingdbs:
    #     print(gb)
    #     arcpy.env.workspace = gb
    #     for i in arcpy.ListRasters('CostDis_subp*'):
    #         arcpy.Delete_management(i)

    accesscalc_chunkedir2(ingdbs[3], inyears = analysis_years, outgdb = outstats)


    # for gb in ingdbs[2:]:
    #     accesscalc_chunkedir(gb,
    #                          inyears = analysis_years,
    #                          outgdb = outstats)


    inyears = ['2000', '2010', '2018']
    outgdb=outstats
    arcpy.env.workspace = inchunkgdb
    rootname = os.path.splitext(os.path.split(inchunkgdb)[1])[0]
    llhood = re.search('[aA-zZ]*[_][aA-zZ]*', rootname).group()
    bufferrad = re.search('[0-9]*(?=_[0-9]$)', rootname).group()
    inpoints = 'subpoints{}'.format(rootname)
    barrierweight_outras = {year: '{0}_bw1_{1}'.format(llhood, year) for year in inyears}
    forestyearly = {year: 'Hansen_GFC_v16_treecover{}'.format(year) for year in inyears}

    if 'group{}'.format(llhood) not in [i.name for i in arcpy.Describe(inpoints).indexes]:
        arcpy.AddIndex_management(inpoints, fields='group{}'.format(llhood),
                                  index_name='group{}'.format(llhood))  # Add index to speed up querying

    grouplist = {row[0] for row in arcpy.da.SearchCursor(inpoints, 'group{}'.format(llhood))}
    group = grouplist[0]
    arcpy.MakeFeatureLayer_management(in_features=inpoints, out_layer='pointslyr_{}'.format(group),
                                      where_clause='{0} = {1}'.format('group{}'.format(llhood), group))

    cProfile.run("""
    accesscalc2(inllhood=llhood,
                ingroup=group,
                inpoints='pointslyr_{}'.format(group),
                inbuffer_radius=bufferrad,
                inyears=inyears,
                inbarrierweight_outras=barrierweight_outras,
                inforestyearly=forestyearly,
                costtab_outgdb=outgdb)""")

    cProfile.run("""
    accesscalc2_jit(inllhood=llhood,
                ingroup=group,
                inpoints='pointslyr_{}'.format(group),
                inbuffer_radius=bufferrad,
                inyears=inyears,
                inbarrierweight_outras=barrierweight_outras,
                inforestyearly=forestyearly,
                costtab_outgdb=outgdb)""")

    #Extra stuff
    # def access_mapalgebra(inllhood, inforestyearly, inbarrierweight_outras, year):
    #     if inllhood == 'Charcoal_production':
    #         # forest resource weighting*(1/(1+cost))
    #         return (Int(100 * Raster(inforestyearly[year]) * \
    #                     (1 / (1 + CostDistance(in_source_data=inpoints,
    #                                            in_cost_raster=inbarrierweight_outras[year])))))
    #     else:
    #         # (1/(1+cost))
    #         return(Int(100 * (1 / (1 + CostDistance(in_source_data=inpoints,
    #                                                 in_cost_raster=inbarrierweight_outras[year])))))
