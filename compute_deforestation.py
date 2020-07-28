#-----------------------------------------------------------------------------------------------------------------------
# Analysis_Chp1 (2020)
# Title: Compute deforestation
# del Giorgio, Olivia (https://github.com/odelgi)
# Messager, Mathis (https://github.com/messamat)

# Compute percent deforestation within livelihood-specific buffer around every location in Department of Pellegrini
# by computing a running mean of binary forest cover raster using focal statistics.
#-----------------------------------------------------------------------------------------------------------------------

from accesscalc_parallel import *
import pandas as pd

#-----------------------------------------------------------------------------------------------------------------------
### Set environment settings
# https://pro.arcgis.com/en/pro-app/arcpy/classes/env.htm
arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

#-----------------------------------------------------------------------------------------------------------------------
#### Set directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

basemapgdb = os.path.join(resdir, "Base_layers_Pellegrini", "Basemaps_UTM20S.gdb")
pelleras = os.path.join(basemapgdb, "Pellegri_department_UTM20S")
forestoutdir = os.path.join(resdir, 'Base_layers_Pellegrini/Forest_Hansen.gdb')

weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

#Define years to analyze
analysis_years = ['2000', '2010', '2018']

#-----------------------------------------------------------------------------------------------------------------------
#Compute deforestation index for each year
arcpy.env.snapRaster = pelleras

#For each livelihood-specific buffer size
for bufferad in list(weightingpd['Buffer_max_rad'])[:-1]:
    defo_outgdb = os.path.join(resdir, 'Analysis_Chp1_W1', 'Deforestation_index', 'Deforestation_index{}k.gdb'.format(bufferad/1000))
    pathcheckcreate(defo_outgdb)

    #For each analysis year
    for year in analysis_years:
        outdefo = os.path.join(defo_outgdb, 'deforestation_index{0}k_{1}'.format(bufferad/1000, year))

        if not arcpy.Exists(outdefo):
            print('Processing {}...'.format(outdefo))
            forestyearly = os.path.join(forestoutdir, 'Hansen_GFC_v16_treecover{}'.format(year))

            #Compute percent forested within radius around each pixel in area
            Con(~IsNull(pelleras),
                FocalStatistics(in_raster=Raster(forestyearly),
                                neighborhood=NbrCircle(5000, "MAP"),
                                statistics_type="MEAN",
                                ignore_nodata="DATA")
                ).save(outdefo)
