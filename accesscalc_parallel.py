from multiprocessing import cpu_count
import psutil
from functools import partial
from internal_functions import *
from pebble import ProcessPool

arcpy.env.overwriteOutput = 'True'
arcpy.CheckOutExtension("Spatial")

#To do:
#   - Add logging: from multiprocessing_logging import install_mp_handler
#   - Convert to Snakemake?

if __name__ == '__main__':
    #Define directory structure
    formatdir = os.path.join(os.path.dirname(os.path.abspath(__file__)).split('\\src')[0])#To update for final run
    datadir = os.path.join(formatdir, 'data')
    resdir = os.path.join(formatdir, 'results')
    arcpy.env.workspace = datadir
    ingdbs = arcpy.ListWorkspaces('*', workspace_type='FileGDB')
    outstats = os.path.join(resdir, 'outstats.gdb')
    pathcheckcreate(outstats)

    #---- Run in parallel for all gdb -----#
    analysis_years = ['2000', '2010', '2018'] #Define years of analysis

    #Get number of cpus to run unto - should work on HPC and personal computer
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK')) if \
        os.environ.get('SLURM_CPUS_PER_TASK') is not None else \
        psutil.cpu_count(logical=False)

    maxruntime = 24*3600  # Define maximum running time of worker processes

    print('Launch parallel processing on {0} cores with {1}s timeout...'.format(ncpus, maxruntime))
    with ProcessPool(max_workers=ncpus) as pool:  # Create pool of worker processes (with N # of physical cores)
        #Assign chunked gdbs to worker processes, keeping analysis years and outstats arguments constant with 'partial'
        future = pool.map(partial(accesscalc_chunkedir, inyears=analysis_years, outgdb=outstats),
                          ingdbs,
                          chunksize=int(0.5 + len(ingdbs) / float(ncpus)),  # Divide all gdbs among chunks upfront so that timeout doesn't lead to new worker process
                          timeout=maxruntime) #Raise timeout error after maxruntime
        task_done(future) #Return timeout error
