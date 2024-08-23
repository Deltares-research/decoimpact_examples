# D-Eco Impact examples
This github contains Python based examples for runs, pre-, postprocessing for [D-Eco Impact](https://www.deltares.nl/en/software-and-data/products/d-eco-impact). D-Eco Impact is an open source spatial temporal ecological postprocessing tool and the source code can be found [here](https://github.com/Deltares/D-EcoImpact).
D-Eco Impact is mainly developed to determine the effect of abiotic and biotic conditions on ecological suitability.  

The examples provided on this GitHub have been developed in projects. We encourage you to add your own, so that other users can benefit from these.


## Content

### Pre-processing
D-Eco Impact makes use of the [UGRID NetCDF format](https://ugrid-conventions.github.io/ugrid-conventions/) for input data . Data that is not provided in this format will need to be pre-processed for use in D-Eco Impact. The "preprocessing" folder contains examples for these pre-processing steps.

#### Model output

##### Hydrodynamic models
* [Delft3D 4 results](https://www.deltares.nl/en/software-and-data/products/delft3d-4-suite)
* [SCHISM](https://ccrm.vims.edu/schismweb/)
* [SHYFEM](https://www.cmcc.it/models/shyfem-shallow-water-hydrodynamic-finite-element-model)

##### Hydrological models
* [wflow](https://deltares.github.io/Wflow.jl)

### Examples
Examples for the operation and use of D-Eco Impact to quantify ecological effect. The "examples" folder contains examples for the use of D-Eco Impact.


### Post-processing
D-Eco Impact exports results in the [UGRID NetCDF format](https://ugrid-conventions.github.io/ugrid-conventions/). For presentation or further analysis of these results post-processing is needed. The "postprocessing" folder contains examples for these postprocessing steps.

#### Export results to ESRI shapefile
Usefull for visual spatial presentation in GIS software. Please note that the time or depth varying aspect of the UGRID file can not be exported to this format. 

### Tests
Subsetted data input examples from projects used to demonstrate D-Eco Impact. The "tests" folder contains example input files for testing the use of D-Eco Impact.

### Documentation
Documentation of tools and applications of D-Eco Impact that decribe D-Eco Impact, are helpfull for cloud computing, visualistion and data inspection. The "documentation" folder contains documentation usefull for D-Eco Impact.

## Acknowledgements
These examples make substantial use of the following packages, which have proven to be beneficial to pre- and postprocessing D-EcoImpact data:
* [xarray](https://docs.xarray.dev/)
* [netcdf4](https://unidata.github.io/netcdf4-python/)
* [xugrid](https://deltares.github.io/xugrid/)
* [DFMtools](https://deltares.github.io/dfm_tools/)
* [HydroMT](https://github.com/Deltares/hydromt)

