import arcpy
import os
import pandas as pd
from collections import defaultdict
from internal_functions import *

formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0])  # To update for final run
resdir = os.path.join(formatdir, 'results')
outstats = os.path.join(resdir, 'outstats.gdb')

tablist = getfilelist(outstats, '*')
tables_pd = pd.concat(pd.Series(tablist),
                      pd.Series(tablist).str.split('_', expand=True),
                      axis=1)
tables_pd.columns = ['path', 'llhood', 'year', 'weighting', 'group']



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