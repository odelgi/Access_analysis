from accesscalc_parallel import *

if __name__ == '__main__':
    ####### TROUBLESHOOTING AND PROFILING CODE #############
    import cProfile

    arcpy.env.overwriteOutput = 'True'
    arcpy.CheckOutExtension("Spatial")

    # Define directory structure
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0])  # To update for final run
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
    inchunkgdb = ingdbs[23]
    inyears = analysis_years
    outdir = resdir
    accesscalc_chunkedir(inchunkgdb, inyears, outdir)


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
