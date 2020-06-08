import arcpy
import os
import pandas as pd
from collections import defaultdict
from accesscalc_parallel import *

if __name__ == 'main':
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0])  # To update for final run
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    llhoodbuffer_outdir = os.path.join(resdir, 'Analysis_Chp1_W1','Buffers_W1.gdb')  # Output buffer path for each livelihood
    basemapgdb = os.path.join(resdir, "Base_layers_Pellegrini", "Basemaps_UTM20S.gdb")
    pelleras = os.path.join(basemapgdb, "Pellegri_department_UTM20S")
    pellepoints = os.path.join(llhoodbuffer_outdir, 'pellefishpoints')
    outstats = os.path.join(resdir, 'Analysis_Chp1_W1/W1_3030/Cost_distance_W1_3030')
    outmerge = os.path.join(resdir, 'outstats_merge.gdb')
    pathcheckcreate(outmerge)

    refraster = pelleras

    weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
    weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

    livelihoods = weightingpd['Livelihood'].tolist()
    livelihoods.remove('Combined_livelihood')

    #Get panda df of tables
    print('Getting all zonal statistics tables...')
    tablist = getfilelist(dir=outstats, repattern=".*[.]dbf$", gdbf=False, nongdbf=True)
    tablist.extend(getfilelist(dir=outstats, gdbf=True, nongdbf=False))

    tables_pd = pd.concat([pd.Series(tablist),
                          pd.Series(tablist).apply(lambda x: os.path.splitext(os.path.split(x)[1])[0]).
                          str.split('_', expand=True)],
                          axis=1)
    tables_pd.columns = ['path', 'dataset', 'llhood1', 'llhood2', 'year', 'weighting', 'group']
    tables_pd['llhood'] = tables_pd['llhood1'] + '_' + tables_pd['llhood2']
    tables_pd = tables_pd.drop(labels=['llhood1', 'llhood2'], axis=1)

    tables_pd.groupby('year')['path'].nunique()
    processed_pd = tables_pd.groupby(['llhood', 'group']).filter(lambda x: x['year'].nunique() == 3). \
        drop_duplicates(subset=['llhood', 'group'])

    #Create a raster of access for each llhood and year by aggregating all tables (yielding an access value for each pixel-point)
    access_outgdb = {}

    #Iterate over each livelihood
    for llhood in tables_pd['llhood'].unique():
        access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                             'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
        pathcheckcreate(access_outgdb[llhood])
        #Ierate over each year
        for year in tables_pd['year'].unique():
            #Path of output access raster
            access_outras = os.path.join(access_outgdb[llhood], 'accessras_W1_{0}{1}'.format(llhood, year))

            #Perform analysis only if output raster doesn't exist
            if not arcpy.Exists(access_outras):
                print("Processing {}...".format(access_outras))

                #Aggregate values across all pixels-points for that livelihood-year
                print('Aggregating zonal statistics tables...')
                merged_dict = tabmerge_dict(tables_pd.loc[(tables_pd['llhood']==llhood) &
                                                          (tables_pd['year']==year), 'path'])

                #Join all statistics tables of access to pellepoints (a point for each 30x30 m pixel in Pellegrini department)
                print('Joining tables to points...')
                accessfield = 'access{0}{1}'.format(llhood, year)
                if not accessfield in [f.name for f in arcpy.ListFields(pellepoints)]:
                    print('Create {} field'.format(accessfield))
                    arcpy.AddField_management(in_table=pellepoints, field_name=accessfield, field_type='FLOAT')

                    with arcpy.da.UpdateCursor(pellepoints, ['pointid', accessfield]) as cursor:
                        x = 0
                        for row in cursor:
                            if x % 100000 == 0:
                                print(x)
                            try:
                                row[1] = merged_dict[row[0]]
                            except:
                                print('pointid {} was not found in dictionary'.format(row[0]))
                            cursor.updateRow(row)
                            x += 1

                # Convert points back to raster
                print('Converting points to raster...')
                arcpy.PointToRaster_conversion(in_features=pellepoints, value_field=accessfield, cellsize=refraster,
                                               out_rasterdataset=access_outras)

            else:
                print('{} already exists...'.format(access_outras))