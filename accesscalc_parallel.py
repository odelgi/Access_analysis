import multiprocessing
import psutil
from functools import partial
from internal_functions import *

arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

#To do:
#   - Add timeout
#   - Add logging: from multiprocessing_logging import install_mp_handler
#   - Convert to Snakemake?

if __name__ == '__main__':
    #Define directory structure
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0],
                             'results\Analysis_Chp1_W1\inputdata_HPC/Beluga_input20200331')#To update for final run
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    arcpy.env.workspace = datadir
    ingdbs = arcpy.ListWorkspaces('*', workspace_type='FileGDB')
    outstats = os.path.join(resdir, 'outstats.gdb')
    pathcheckcreate(outstats)

    #Run in parallel for all gdb
    analysis_years = ['2000', '2010', '2018']
    print('Launch parallel processing')
    tic = time.time()
    p = multiprocessing.Pool(psutil.cpu_count(logical=False)) #Create pool of worker processes (with N # of physical cores)
    accesscalc_chunkedir_partial = partial(accesscalc_chunkedir,
                                           inyears=analysis_years,
                                           outgdb =outstats)
    p.map(accesscalc_chunkedir_partial, ingdbs)
    p.close()
    print(time.time() - tic)
