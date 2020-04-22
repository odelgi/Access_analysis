import arcpy
import os
import pandas as pd
from collections import defaultdict
from accesscalc_parallel import *

def tabmerge_dict(tablist):
    outdict = {}
    for tab in tablist:
        print(tab)
        for row in arcpy.da.SearchCursor(tab, ['Value', 'MEAN']):
            outdict[row[0]] = row[1]
    return(outdict)

if __name__ == 'main':
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0])  # To update for final run
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    llhoodbuffer_outdir = os.path.join(resdir, 'Analysis_Chp1_W1',
                                       'Buffers_W1.gdb')  # Output buffer path for each livelihood
    pellepoints = os.path.join(llhoodbuffer_outdir, 'pellefishpoints')
    outstats = os.path.join(resdir, 'Analysis_Chp1_W1/W1_3030/Cost_distance_W1_3030')
    outmerge = os.path.join(resdir, 'outstats_merge.gdb')
    pathcheckcreate(outmerge)

    weighting_table = os.path.join(datadir, 'Weighting_scheme.xlsx')
    weightingpd = pd.read_excel(weighting_table, sheetname='Weighting_1')

    livelihoods = weightingpd['Livelihood'].tolist()
    livelihoods.remove('Combined_livelihood')

    #Get panda df of tables
    # arcpy.env.workspace = resdir
    # tablist = arcpy.ListTables(table_type='dBASE')
    tablist = getfilelist(dir=outstats, repattern= ".*[.]dbf$")

    tables_pd = pd.concat([pd.Series(tablist),
                          pd.Series(tablist).apply(lambda x: os.path.splitext(os.path.split(x)[1])[0]).
                          str.split('_', expand=True)],
                          axis=1)
    tables_pd.columns = ['path', 'dataset', 'llhood1', 'llhood2', 'year', 'weighting', 'group']

    tables_pd.groupby('year')['path'].nunique()


    tablist2 = getfilelist(dir="D:/Users/odelgi/MSwork/Black_input20200419/results", repattern= ".*[.]dbf$")

    merged_dict = {}
    for llhood in tables_pd['llhood1'].unique():
        for year in tables_pd['year'].unique():
            print(year)
            merged_dict[llhood+year] = tabmerge_dict(tables_pd.loc[(tables_pd['llhood1']==llhood) &
                                                                   (tables_pd['year']==year), 'path'])


    access_outgdb = {}
    for llhood in livelihoods:
        # Output access gdb
        access_outgdb[llhood] = os.path.join(resdir, 'Analysis_Chp1_W1', 'W1_3030', 'Access_W1_3030',
                                             'Access_W1_{0}'.format(llhood), 'Access_W1_{0}.gdb'.format(llhood))
        pathcheckcreate(access_outgdb[llhood])

        for year in tables_pd['year'].unique():
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
                            row[1] = merged_dict[llhood+year][row[0]]
                        except:
                            print('pointid {} was not found in dictionary'.format(row[0]))
                        cursor.updateRow(row)
                        x += 1

        # Convert points back to raster
        arcpy.PointToRaster_conversion(in_features=pellepoints, value_field=accessfield, cellsize=refraster,
                                       out_rasterdataset=access_outras)

#############
for k, v in processeddict.iteritems():
    print(k)
    print(len(v))

for k, v in groupstoprocess.iteritems():
    print(k)
    print(len(v))