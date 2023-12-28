## A livelihood-specific access index for the analysis of the social impacts of commodity frontier expansion in the Argentine Gran Chaco
We developed a novel measure of access to rural livelihood opportunities to analyze 
the social impacts of commodity frontier expansion. This approach is demonstrated here
in Python and implemented for the Department of Pellegrini in the Argentine Gran Chaco.

This repository accompanies the following peer-reviewed publication 
[del Giorgio O., Messager M. L., le Polain de Waroux Y (2021) Fenced off: Measuring growing restrictions on resource access for smallholders in the Argentine Chaco. Applied Geography](https://www.sciencedirect.com/science/article/pii/S0143622821001466).

## Project description
We propose a novel approach to evaluate changes in access to land for smallholders stemming from gradual changes in land control along commodity frontiers. We apply this approach in the Argentine Gran Chaco, a region that has experienced amongst the highest global rates of deforestation for agriculture in recent years. Our findings suggest that access to natural resources for smallholders has been reduced far beyond what would be expected if only looking at deforestation, and that the degree to which access has decreased differs between livelihood activities.    
[Method diagram]

## Getting Started
#### Software prerequisites (what we ran this on)
* 64-bit Windows 10 
* Python 2.7 as part of ESRI ArcGIS for Desktop 10.8 
[Background Geoprocessing (64 Bit)](https://desktop.arcgis.com/en/arcmap/10.3/analyze/executing-tools/64bit-background.htm)
* Python IDE (e.g. [PyCharm](https://www.jetbrains.com/pycharm/promo/?gclid=CjwKCAjwmf_4BRABEiwAGhDfSaGbHthcudKiCTLaWZWj7-MdcAC_4mIWlJ8vWtQTbBAFptMMC_8hvBoC6v4QAvD_BwE))

#### Hardware recommendation
* This analysis was run at 30 m resolution for the Department of Pellegrini, which spans 7,330 km<sup>2</sup>.
* All of the analysis was run on a workstation with the following specs:
    * Intel® Core™ i9-9940X X-series Processor (14 physical cores, 19.25M Cache, up to 4.50 GHz)
    * 64 GB RAM
    * Samsung SSD 970 ECO Plus 500GB
* We estimate that, without significant issues, the analysis would take between 45 and 60 days to run 24/7 on all 14 cores.

#### USGS API access
To be done. Write a config.json file.

#### Input data
To obtain all input data and products, download the accompanying [data repository](url).

The scripts take care of downloading reference Landsat imagery and Forest Loss data 

#### Project structure
The scripts contained in this Github repository will only function if contextualized within the specific directory
structure provided in the data repository (`data`,`results`, and `src`subdirectories form the basis of the project folder structure, 
following [Wilson et al. 2017](https://doi.org/10.1371/journal.pcbi.1005510) and others). The scripts will automatically
detect the absolute path of the project directory and conduct analysis using all relative path otherwise.

The data repository already includes a snapshot of the scripts used to run the 
analysis at the time of the publication of the accompanying peer-reviewed article. If there have been updates to the source code since publication; these will be detailed in this README. 
In this case, download this repository into the `src` folder within the
_[Updates since peer-reviewed article publication](#Updates since peer-reviewed article publication)_ section of this 
README file.

## Workflow
The core of this analysis is based on a "cost distance" calculation. The cost distance tells us the difficulty 
for a smallholder to access the resources necessary to perform a given livelihood for each 30x30 m pixel in the 
department of Pellegrini. Given the degree of hindrance to movement that difference types of barriers represent for a 
given livelihood, a cost distance calculation computes the least-accumulative cost pasth from each pixel to all pixels
within a livelihood-specific distance. For instance, etc. -- to edit, it's kind of gibberish right now.

This methodology implies that the cost distance of neighboring pixels cannot be concurrently computed
(in one single raster processing step), as the area of use of neighboring pixels overlap. Therefore, calculating cost
distance for every pixel is very computing intensive and requires us to divide the area in non-overlapping groups of 
locations so that the buffers around the locations within a given group don't overlap. Further, because the size of 
the area of use differs between livelihoods, different groupings are performed for different livelihoods. Livelihoods
requiring a larger area of use require dividing the area into a greater number of groups, each containing a smaller 
number of locations.

Further, to speed up processing, a parallel processing workflow was implemented (dividing the task into many chunks to 
be analyzed on different processor cores at the same time). This required further separating groups
of locations into different chunks so that each chunk could be processed by a single core.

The analysis requires that the scripts be run in the following order 
(more detail on workflow within each script is available therein):
1. [main_analysis_1.py](https://github.com/odelgi/Access_analysis/blob/master/main_analysis_1.py) to dowload and 
pre-process data.
2. [mergetables_rechunk.py](https://github.com/odelgi/Access_analysis/blob/master/mergetables_rechunk.py) to prepare data
for parallel processing. 
3. [accesscalc_parallel.py](https://github.com/odelgi/Access_analysis/blob/master/accesscalc_parallel.py) to run cost 
distance calculation and compute access statistics in parallel. This script can only be run from command line or
 imported for functions.
4. [mergetables_rechunk.py](https://github.com/odelgi/Access_analysis/blob/master/mergetables_rechunk.py) again to
 compile cost distance statistics across all runs, livelihoods, and years, generate continuous access rasters, determine
 the locations for which access still needs to be calculated, and re-prepare data (chunking) for re-running 
 [accesscalc_parallel.py](https://github.com/odelgi/Access_analysis/blob/master/accesscalc_parallel.py).
 
    Because this analysis is so computing intensive, it will probably take a few cycles of running [mergetables_rechunk.py](https://github.com/odelgi/Access_analysis/blob/master/mergetables_rechunk.py)
    and  [accesscalc_parallel.py](https://github.com/odelgi/Access_analysis/blob/master/accesscalc_parallel.py) before obtaining
    access statistics for the entire study area. 
    Once this is done. Run the next steps.
5. [compute_deforestation.py](https://github.com/odelgi/Access_analysis/blob/master/compute_deforestation.py) to compute
a conventional measure of access loss based only on a binary map of deforestation.
6. [get_communitystats.py](https://github.com/odelgi/Access_analysis/blob/master/get_communitystats.py) to compute yearly 
access and deforestation indices averaged over the geographic span of community polygons.

## Possible speed-ups
To be completed. 

## Updates since peer-reviewed article publication
None

## Code developers
* [Mathis Loïc Messager](https://github.com/messamat)
* [Olivia del Giorgio](https://github.com/odelgi)

## License
To be completed:
A short snippet describing the license (MIT, Apache etc)
MIT © [Yourname]()

## Acknowledgments
* Tim Elrick and the McGill Geographic Information Center
