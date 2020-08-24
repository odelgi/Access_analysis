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

commugdb_byyear = os.path.join(datadir, 'Community_boundaries' , 'Developed_areas.gdb')
commustats_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Communities_access')
commudefo_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Communities_deforestation')

hhpts = os.path.join(datadir, 'Household_points/hh_points.gdb/hh_points')
hhstats_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Households_access')
hhdefo_outdir = os.path.join(resdir, 'Analysis_Chp1_W1', 'Households_deforestation')

#Define livelihoods to analyze
weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

livelihoods = weightingpd['Livelihood'].tolist()
livelihoods.remove('Combined_livelihood')

#Define years to analyze
analysis_years = ['2000', '2010', '2018']

#Get list of community boundaries
commupoly_list = [os.path.join(commugdb_byyear,  'Developed_areas_{}a'.format(y)) for y in analysis_years]

### -------------------- GET ACCESS STATISTICS FOR COMMUNITIES AND HOUSEHOLDS ------------------------------------------####
#Get paths to livelihood-year specific access rasters
access_outgdb = {}
access_outras = {}
overwritestats = True

for llhood in livelihoods:
    print(llhood)
    commullhoodstats_outgdb = os.path.join(commustats_outdir, 'commuaccess_{}.gdb'.format(llhood))
    pathcheckcreate(commullhoodstats_outgdb)

    hhllhoodstats_outgdb = os.path.join(hhstats_outdir, 'hhaccess_{}.gdb'.format(llhood))
    pathcheckcreate(hhllhoodstats_outgdb)

    for year in analysis_years:
        # Path of access raster
        access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                             'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
        access_outras[llhood + year] = os.path.join(access_outgdb[llhood],
                                                    'accessras_W1_{0}{1}'.format(llhood, year))

        #If access raster exists
        if arcpy.Exists(access_outras[llhood + year]):
            arcpy.env.snapRaster = access_outras[llhood + year]

            # For each yearly community/developed area polygon set and household points
            for commupoly in commupoly_list:
                #Name of output table
                outstats_commu = os.path.join(commullhoodstats_outgdb,
                                        '{0}_access_{1}{2}'.format(os.path.split(commupoly)[1], llhood, year))
                print('Processing {}...'.format(outstats_commu))
                #Compute zonal stats
                ZonalStatisticsAsTable(in_zone_data=commupoly,
                                       zone_field=arcpy.Describe(commupoly).OIDFieldName,
                                       in_value_raster=access_outras[llhood+year],
                                       out_table=outstats_commu,
                                       ignore_nodata="DATA")

            outstats_hh = os.path.join(hhllhoodstats_outgdb,
                                        'households_access_{0}{1}'.format(llhood, year))
            ZonalStatisticsAsTable(in_zone_data=hhpts,
                                   zone_field=arcpy.Describe(hhpts).OIDFieldName,
                                   in_value_raster=access_outras[llhood + year],
                                   out_table=outstats_hh,
                                   ignore_nodata="DATA")

        else:
            print('{} does not exist...'.format(access_outras[llhood + year]))

### ---------------------------------- GET DEFORESTATION STATISTICS ------------------------------------------------####
#Get paths to livelihood-year specific access rasters
for bufferad in list(weightingpd['Buffer_max_rad'])[:-1]:
    defo_outgdb = os.path.join(resdir, 'Analysis_Chp1_W1', 'Deforestation_index', 'Deforestation_index{}k.gdb'.format(bufferad/1000))
    commudefostats_outgdb = os.path.join(commudefo_outdir, 'commudeforestation_{}k.gdb'.format(bufferad/1000))
    pathcheckcreate(commudefostats_outgdb)

    hhdefostats_outgdb = os.path.join(hhdefo_outdir, 'hhdeforestation_{}k.gdb'.format(bufferad/1000))
    pathcheckcreate(hhdefostats_outgdb)


    for year in analysis_years:
        deforas = os.path.join(defo_outgdb, 'deforestation_index{0}k_{1}'.format(bufferad/1000, year))

        if arcpy.Exists(deforas):
            arcpy.env.snapRaster = deforas

            #Path of input community/developed area polygons (not merged) for that year
            for commupoly in commupoly_list:
                    #Name of output table
                    outstats_commu = os.path.join(commudefostats_outgdb,
                                            '{0}_deforestation_{1}k_{2}'.format(os.path.split(commupoly)[1], bufferad/1000, year))
                    print('Processing {}...'.format(outstats_commu))
                    #Compute zonal stats
                    ZonalStatisticsAsTable(in_zone_data=commupoly,
                                           zone_field=arcpy.Describe(commupoly).OIDFieldName,
                                           in_value_raster=deforas,
                                           out_table=outstats_commu,
                                           ignore_nodata="DATA")

            #Process stats for households
            outstats_hh = os.path.join(hhdefostats_outgdb,
                                        'households_deforestation_{0}k_{1}'.format(bufferad/1000, year))
            ZonalStatisticsAsTable(in_zone_data=hhpts,
                                   zone_field=arcpy.Describe(hhpts).OIDFieldName,
                                   in_value_raster=deforas,
                                   out_table=outstats_hh,
                                   ignore_nodata="DATA")

        else:
            print('{} does not exist...'.format(deforas))