import os
import sys
import gc
import traceback
import re
import arcpy
from arcpy.sa import *
import time
import numpy as np
from concurrent.futures import TimeoutError
from multiprocessing import cpu_count
import psutil
from functools import partial
from pebble import ProcessPool

arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

#To do:
#   - Add logging: from multiprocessing_logging import install_mp_handler
#   - Convert to Snakemake?

#To launch:
#   If python not in environment variables, in terminal:
#   C;
#   cd \Python27\ArcGISx6410.7
#   python D:\Users\odelgi\MSwork\Black_input20200402\src\accesscalc_parallel.py #for black
#   python C:\Users\odelgi\White_input20200518\src\accesscalc_parallel.py

# Functions
#Get all files in a ArcGIS workspace (file or personal GDB)
def getwkspfiles(dir, repattern=None):
    arcpy.env.workspace = dir
    filenames_list = (arcpy.ListDatasets() or []) + (arcpy.ListTables() or [])  # Either LisDatsets or ListTables may return None so need to create empty list alternative
    if not repattern == None:
        filenames_list = [filen for filen in filenames_list if re.search(repattern, filen)]
    return ([os.path.join(dir, filen) for filen in filenames_list])
    arcpy.ClearEnvironment('workspace')

def getfilelist(dir, repattern=None, gdbf=True, nongdbf=True):
    """Function to iteratively go through all subdirectories inside 'dir' path
    and retrieve path for each file that matches "repattern"
    gdbf and nongdbf allows the user to choose whether to consider ArcGIS workspaces (GDBs) or not or exclusively"""

    try:
        if arcpy.Describe(dir).dataType == 'Workspace':
            if gdbf == True:
                print('{} is ArcGIS workspace...'.format(dir))
                filenames_list = getwkspfiles(dir, repattern)
            else:
                raise ValueError(
                    "A gdb workspace was given for dir but gdbf=False... either change dir or set gdbf to True")
        else:
            filenames_list = []

            if gdbf == True:
                for (dirpath, dirnames, filenames) in os.walk(dir):
                    for in_dir in dirnames:
                        fpath = os.path.join(dirpath, in_dir)
                        if arcpy.Describe(fpath).dataType == 'Workspace':
                            print('{} is ArcGIS workspace...'.format(fpath))
                            filenames_list.extend(getwkspfiles(dir=fpath, repattern=repattern))

            if nongdbf == True:
                for (dirpath, dirnames, filenames) in os.walk(dir):
                    for file in filenames:
                        if repattern is None:
                            filenames_list.append(os.path.join(dirpath, file))
                        else:
                            if re.search(repattern, file):
                                filenames_list.append(os.path.join(dirpath, file))
        return (filenames_list)

    # Return geoprocessing specific errors
    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages(2))
    # Return any other type of error
    except:
        # By default any other errors will be caught here
        e = sys.exc_info()[1]
        print(e.args[0])

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

def accesscalc(inllhood, ingroup, inpoints, inbuffer_radius, inyears, inbarrierweight_outras,
               inforestyearly, costtab_outdir):
    # Buffer subsetted points based on livelihood-specific buffer distance
    #print('Bufferring...')
    try:
        arcpy.Buffer_analysis(in_features=inpoints,
                              out_feature_class=r'in_memory/subbuffers{}'.format(ingroup),
                              buffer_distance_or_field=inbuffer_radius,
                              dissolve_option='NONE',
                              method='PLANAR')

        templateras = inbarrierweight_outras[inyears[0]]
        arcpy.env.SnapRaster = templateras
        buffras_memory = r'in_memory/subbufferas{}'.format(ingroup)
        #print('Conversion to raster...')
        arcpy.FeatureToRaster_conversion(in_features=r'in_memory/subbuffers{}'.format(ingroup),
                                         field='pointid',
                                         out_raster=buffras_memory,
                                         cell_size=templateras)

        # Iterate through years of analysis
        #print('Iterating through years...')
        arcpy.env.mask = buffras_memory
        tempaccessdict = {}
        for year in inyears:
            outtab_year = os.path.join(costtab_outdir, 'CD_{0}_{1}_w1_{2}.dbf'.format(inllhood, year, ingroup))
            dictkey = '{0}{1}'.format(ingroup, year)
            if not arcpy.Exists(outtab_year):
                #print(year)
                # Compute cost distance and access
                if inllhood == 'Charcoal_production':
                    # Round(100*forest resource weighting*(1/(1+cost))
                    tempaccessdict[dictkey] = \
                        Int(100 * Raster(inforestyearly[year]) * \
                            (1 / (1 + CostDistance(in_source_data=inpoints,
                                                   in_cost_raster=Raster(inbarrierweight_outras[year])+0.0001))) +
                            0.5)
                else:
                    # Round(100*(1/(1+cost))
                    tempaccessdict[dictkey] = \
                        Int(100 * \
                            (1 / (1 + CostDistance(in_source_data=inpoints,
                                                   in_cost_raster=Raster(inbarrierweight_outras[year])+0.0001))) +
                            0.5)

                # Zonal statistics based on buffer (using pointid, the unique ID of each point for that livelihood)
                # Compute mean access within livelihood-specific buffer and writes it out to table
                #accessras.save(os.path.join(costtab_outdir, 'test_{0}_{1}_w1_{2}'.format(inllhood, year, ingroup))
                ZonalStatisticsAsTable(in_zone_data=buffras_memory,
                                       zone_field='Value',
                                       in_value_raster=tempaccessdict[dictkey],
                                       out_table=outtab_year,
                                       ignore_nodata='DATA',
                                       statistics_type='MEAN')

            #Clean up
            try:
                arcpy.Delete_management(tempaccessdict[dictkey])
                del tempaccessdict[dictkey]
            except:
                pass

        arcpy.ClearEnvironment("mask")
        arcpy.Delete_management(buffras_memory)
        arcpy.Delete_management("in_memory")
        del tempaccessdict

    except:
        #Clean up memory, environment, and temporary files
        arcpy.ClearEnvironment("mask")
        arcpy.Delete_management("in_memory")

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Print Python error messages
        print("PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1]))
        print("ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n")

        try:
            arcpy.Delete_management(tempaccessdict[dictkey])
            del tempaccessdict[dictkey]
        except:
            pass

def accesscalc_chunkedir(inchunkgdb, inyears, outdir):
    arcpy.env.workspace = inchunkgdb
    arcpy.env.scratchWorkspace = 'in_memory'

    rootname = os.path.splitext(os.path.split(inchunkgdb)[1])[0]
    llhood = re.search('[aA-zZ]*[_][aA-zZ]*', rootname).group()
    bufferrad = re.search('[0-9]*(?=_[0-9]*$)', rootname).group()
    inpoints = 'subpoints{}'.format(rootname)
    barrierweight_outras = {year: '{0}_bw1_{1}'.format(llhood, year) for year in inyears}
    forestyearly = {year: 'Hansen_GFC_v16_treecover{}'.format(year) for year in inyears}
    print(rootname)

    outdir_chunk = os.path.join(outdir, rootname)
    if not os.path.exists(outdir_chunk):
        os.mkdir(outdir_chunk)

    if 'group{}'.format(llhood) not in [i.name for i in arcpy.Describe(inpoints).indexes]:
        arcpy.AddIndex_management(inpoints, fields='group{}'.format(llhood),
                                  index_name='group{}'.format(llhood))  # Add index to speed up querying

    # Iterate through groups
    grouplist = {row[0] for row in arcpy.da.SearchCursor(inpoints, 'group{}'.format(llhood))}
    for group in grouplist:
        #If all output tables don't already exist
        if not all(arcpy.Exists(os.path.join(outdir_chunk, 'CD_{0}_{1}_w1_{2}.dbf'.format(llhood, year, group)))
                   for year in inyears):
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
                       costtab_outdir=outdir_chunk)

            arcpy.Delete_management('pointslyr_{}'.format(group))
            #print(time.time() - tic)
        else:
            print('Group {} was already analyzed...'.format(group))

        gc.collect()

    del grouplist
    arcpy.Delete_management("in_memory")
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
    #Define directory structure
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0])#To update for final run
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    arcpy.env.workspace = datadir
    ingdbs = arcpy.ListWorkspaces('*', workspace_type='FileGDB')

    #---- Run in parallel for all gdb -----#
    analysis_years = ['2000', '2010', '2018'] #Define years of analysis

    #Get number of cpus to run unto - should work on HPC and personal computer
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK')) if \
        os.environ.get('SLURM_CPUS_PER_TASK') is not None else \
        psutil.cpu_count(logical=False)

    maxruntime = 7*24*3600  # Define maximum running time of worker processes

    print('Launch parallel processing on {0} cores with {1}s timeout...'.format(ncpus, maxruntime))
    with ProcessPool(max_workers=ncpus) as pool:  # Create pool of worker processes (with N # of physical cores)
        #Assign chunked gdbs to worker processes, keeping analysis years and outstats arguments constant with 'partial'
        future = pool.map(partial(accesscalc_chunkedir, inyears=analysis_years, outdir=resdir),
                          ingdbs,
                          #chunksize=int(0.5 + len(ingdbs) / float(ncpus)),  # Divide all gdbs among chunks upfront so that timeout doesn't lead to new worker process
                          timeout=maxruntime) #Raise timeout error after maxruntime
        task_done(future) #Return timeout error
