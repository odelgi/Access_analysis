## A livelihood-specific access index for the analysis of the social impacts of commodity frontier expansion in the Argentine Gran Chaco
We developed a novel measure of access to rural livelihood opportunities to analyze 
the social impacts of commodity frontier expansion. This approach is demonstrated here
in Python and implemented for the Department of Pellegrini in the Argentine Gran Chaco.

This repository accompanies the following peer-reviewed publication 
[del Giorgio O., Messager M. L., le Polain de Waroux Y (year) title. _journal_](link) and
the accompanying [data repository](url).

## Project description
[Short abstract-like paragraph]    
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
The scripts contained in this Github repository will only function if contextualized within a specific directory
structure (_dat_, _results_, and _src_ directories form the basis of the project structure, 
following [Wilson et al. 2017](https://doi.org/10.1371/journal.pcbi.1005510) and others), provided in the data repository. 
The data repository already includes a snapshot of the scripts used to run the 
analysis at the time of the publication of the accompanying peer-reviewed article.

If there have been updates to the source code since publication; these will be detailed in this README. 
In this case, download this repository into the /src folder within the data repository.

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



## Possible speed-ups

## Updates since peer-reviewed article publication
None

## Code developers
* [Mathis Loïc Messager](https://github.com/messamat)
* [Olivia del Giorgio](https://github.com/odelgi)

## License
A short snippet describing the license (MIT, Apache etc)

MIT © [Yourname]()

## Acknowledgments
* Tim Elrick and the McGill Geographic Information Center