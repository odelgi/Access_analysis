#-----------------------------------------------------------------------------------------------------------------------
# Analysis_Chp1 (2020)
# Title: Get access and deforestation statistics for communities
# del Giorgio, Olivia (https://github.com/odelgi)
# Messager, Mathis (https://github.com/messamat)

# Compute access and deforestation statistics for each community/developed area in Department of Pellegrini
# Compute statistics for every combination of analysis year and community location across years.
#For example, for a given community it will compute access and deforestation for 200, 2010, and 2018 within the boundaries
#of that community from 2000. This would allow us for example to see whether the community changed a lot because access changed
#within its former boundaries from 2000 to 2010.

#-----------------------------------------------------------------------------------------------------------------------

import arcpy
from arcpy.sa import *
import pandas as pd
from accesscalc_parallel import *

### Set environment settings ###
# https://pro.arcgis.com/en/pro-app/arcpy/classes/env.htm
arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

#### Directory structure ###
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

defo_outgdb = os.path.join(resdir, 'Analysis_Chp1_W1', 'Deforestation_index')

commugdb_clip = os.path.join(datadir, 'Community_boundaries' , 'Developed_areas', 'Developed_areas.gdb')
commustats_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Communities_access')
commudefo_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Communities_deforestation')

#Define livelihoods to analyze
weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

livelihoods = weightingpd['Livelihood'].tolist()
livelihoods.remove('Combined_livelihood')

#Define years to analyze
analysis_years = ['2000', '2010', '2018']

#Get list of community boundaries
commupoly_list = [os.path.join(commugdb_clip,  'Area_comunidades_{}_Clip'.format(y)) for y in analysis_years]

### ---------------------------------- GET ACCESS STATISTICS -------------------------------------------------------####
#Get paths to livelihood-year specific access rasters
access_outgdb = {}
access_outras = {}

for llhood in livelihoods:
    print(llhood)
    commullhoodstats_outgdb = os.path.join(commustats_outdir, 'commuaccess_{}.gdb'.format(llhood))
    pathcheckcreate(commullhoodstats_outgdb)

    for year in analysis_years:
        # Path of access raster
        access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                             'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
        access_outras[llhood + year] = os.path.join(access_outgdb[llhood],
                                                    'accessras_W1_{0}{1}'.format(llhood, year))

        #For each yearly community/developed area polygon set
        for commupoly in commupoly_list:
            if arcpy.Exists(access_outras[llhood+year]):
                arcpy.env.snapRaster = access_outras[llhood + year]
                #Name of output table
                outstats = os.path.join(commullhoodstats_outgdb,
                                        '{0}_access_{1}{2}'.format(os.path.split(commupoly)[1], llhood, year))
                print('Processing {}...'.format(outstats))
                #Compute zonal stats
                ZonalStatisticsAsTable(in_zone_data=commupoly,
                                       zone_field=arcpy.Describe(commupoly).OIDFieldName,
                                       in_value_raster=access_outras[llhood+year],
                                       out_table=outstats,
                                       ignore_nodata="NODATA")
            else:
                print('{} does not exist...'.format(access_outras[llhood + year]))

### ---------------------------------- GET DEFORESTATION STATISTICS ------------------------------------------------####
#Get paths to livelihood-year specific access rasters
for bufferad in list(weightingpd['Buffer_max_rad'])[:-1]:
    defo_outgdb = os.path.join(resdir, 'Analysis_Chp1_W1', 'Deforestation_index', 'Deforestation_index{}k.gdb'.format(bufferad/1000))
    commudefostats_outgdb = os.path.join(commudefo_outdir, 'commudeforestation_{}k.gdb'.format(bufferad/1000))
    pathcheckcreate(commudefostats_outgdb)

    for year in analysis_years:
        deforas = os.path.join(defo_outgdb, 'deforestation_index{0}k_{1}'.format(bufferad/1000, year))
        #Path of input community/developed area polygons (not merged) for that year
        for commupoly in commupoly_list:
            if arcpy.Exists(deforas):
                arcpy.env.snapRaster = deforas
                #Name of output table
                outstats = os.path.join(commudefostats_outgdb,
                                        '{0}_deforestation_{1}k_{2}'.format(os.path.split(commupoly)[1], bufferad/1000, year))
                print('Processing {}...'.format(outstats))
                #Compute zonal stats
                ZonalStatisticsAsTable(in_zone_data=commupoly,
                                       zone_field=arcpy.Describe(commupoly).OIDFieldName,
                                       in_value_raster=deforas,
                                       out_table=outstats,
                                       ignore_nodata="NODATA")
            else:
                print('{} does not exist...'.format(deforas))