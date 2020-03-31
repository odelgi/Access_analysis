import multiprocessing
from functools import partial
from internal_functions import *

arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

if __name__ == '__main__':
    #Define directory structure
    formatdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    arcpy.env.workspace = datadir
    ingdbs = arcpy.ListWorkspaces('*', workspace_type='FileGDB')
    pathcheckcreate(os.path.join(resdir, 'outstats.gdb'))

    #Run in parallel for all gdb
    analysis_years = ['2000', '2010', '2018']
    print('Launch parallel processing')
    tic = time.time()
    p = multiprocessing.Pool(int(multiprocessing.cpu_count()))
    accesscalc_chunkedir_partial = partial(accesscalc_chunkedir,
                                           inyears=analysis_years,
                                           outgdb =pathcheckcreate(os.path.join(resdir, 'outstats.gdb')))
    p.map(accesscalc_chunkedir_partial, ingdbs)
    p.close()
    print(time.time() - tic)
