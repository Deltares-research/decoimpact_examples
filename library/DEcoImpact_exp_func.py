"""Implement DEcoImpact model class"""
# Implement model class following model API

import os
from os.path import join, dirname, basename, isfile, isdir
from typing import Any, Union, Optional
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import logging
import rioxarray  # required for rio accessor

#TODO: temporary implementation outside of HYDROMT-core
import xugrid as xu
import netCDF4 as nc
import pathlib

from hydromt.models.model_mesh import MeshModel

__all__ = ["DEIModel"]

logger = logging.getLogger(__name__)

class DEIOutput(MeshModel):
    """This is the D-Eco Impact output class"""

    _NAME = "decoimpact"
    _CONF = "decoimpact.ini"
#    _DATADIR = DATADIR
    _GEOMS = {}
    _MESHES = {}
    _FOLDERS = [
        "geoms",
    ]

    def __init__(
        self,
        root=None,
        mode="w",
        config_fn=None,
        data_libs=None,
        logger=logger,
    ):
        super().__init__(
            root=root,
            mode=mode,
            config_fn=config_fn,
            data_libs=data_libs,
            logger=logger,
        )
        #TODO : Hydromt-core set correct mode
        self._read = True
        self._mesh2d = None
        self._input_mesh2d_path = None
        self._categories = {}

    def set_paths_to_DEIresults(
            self,
            input_list_path : Union[list, str] = None,           
    ):
        """Method to retrieve the paths DEI results of importance for
        post-processing.
        Parameters
        ----------
        input list_path:
            Only a netcdf path reference(string, ends with .nc), a list of
            netcdfs (list, contains [.nc,.nc,.nc,..]) or a text file
            containing the netcdf paths is allowed as input.
        """
        #get location of netcdf files to convert to shape
        if(isinstance(input_list_path, list)):
            if(input_list_path[0].endswith(".nc")):
                list_netcdf_paths = input_list_path
            else:
                raise ValueError("Only a netcdf(.nc), a list of netcdfs (.nc,.nc,.nc,..)"+\
                                "or a text file containing the netcdf paths is allowed as input. "+\
                                "Now supplied is : " + str(input_list_path))
        elif(isinstance(input_list_path, str)):
            if(input_list_path.endswith(".nc")):
                list_netcdf_paths = input_list_path.split(",")
            elif(input_list_path.endswith(".txt")):
                with open(input_list_path,'r') as f:
                    list_netcdf_paths = [line.rstrip('\n') for line in f]
            else:
                raise ValueError("Only a netcdf(.nc), a list of netcdfs (.nc,.nc,.nc,..)"+\
                            "or a text file containing the netcdf paths is allowed as input. "+\
                            "Now supplied is : " + str(input_list_path))
        else:
            raise ValueError("Only a netcdf(.nc), a list of netcdfs (.nc,.nc,.nc,..)"+\
                            "or a text file containing the netcdf paths is allowed as input. "+\
                            "Now supplied is : " + str(input_list_path))
                            
        #check if path points to valid XUGRID mesh
        for path_found in list_netcdf_paths:
            try:
                xu.open_dataset(path_found)
            except Exception as e:
                raise ValueError("path does not point to a valid Ugrid netCDF file: "+\
                                 str(path_found) +" . Error message while loading :"+\
                                 str(e))

        self._input_mesh2d_path = list_netcdf_paths
        return()

    # COMPONENTS
    def set_as_base_mesh(
        self,
        mesh : Union[xr.DataArray, xu.UgridDataset],
        crs : Union["pyproj.CRS", str] = None,
    ):
        """Method to set current mesh as basis for model input.
        Parameters
        ----------
        mesh:
            A XUGRID UgridDataset
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.      
        """
        
        #TODO Hydromt-core functionality
        #clean up the mesh so all connectivity is present
        #add connectivity for faces
        if(crs == None):
            self.logger.warning("No CRS set")
            pass
        else:
            mesh.ugrid.set_crs(crs=crs, allow_override=False)

        #self.set_mesh(mesh)
        self._mesh2d = mesh
        
        return(mesh)

    def translate_UgridNetCDFs_to_Shape(
            self,
            list_variables, 
            output_path,
    ):
        if(self._input_mesh2d_path == None):
            raise ValueError("Make sure to load input path locations of D-Eco Impact "+\
                             "results to visualize, by running self.set_paths_to_DEIresults()")
        compiled_shapefile_name = "complete_file.shp"
        store_shapefile_names = []
        store_shapecol_names = []
        for variable in list_variables:
            if(len(variable) > 10):
                shapecol_name = variable[:9]
                print("Shortened variable name '" + variable + "' to '"+ shapecol_name)
                
            else:
                shapecol_name = variable
            store_shapecol_names.append(shapecol_name)
        
        for ugrid_ncdf in self._input_mesh2d_path:
            print("Start with "+ugrid_ncdf)
            mesh = xu.open_dataset(
                        ugrid_ncdf,
                    )
            crs_mesh = next(iter(mesh.ugrid.crs.values()))
            polygons_ugrid = mesh.ugrid.grid.to_shapely(
                            dim = mesh.ugrid.grid.face_dimension
                        )
            gdf = gpd.GeoDataFrame(index = range(len(polygons_ugrid)), crs = crs_mesh, geometry=polygons_ugrid)
            for nr, variable in enumerate(list_variables):
                check_for_time = [True if ("time" in dimname) else False for dimname in mesh[variable].dims]
                if(any(check_for_time) and (check_for_time[0]  == True)):
                    print("Time dimension supplied in data variable "+variable+", only first timestep is used for exporting.")
                    gdf[store_shapecol_names[nr]] = mesh[variable][0,:].values
                elif(any(check_for_time) and (check_for_time[0]  == False)):
                    raise ValueError("Time dimension supplied in data variable "+variable+", but is not first dimension.")
                elif(len(mesh[variable].dims) > 1):
                    raise ValueError("More than one dimension supplied in data variable "+variable+", needs to be limited to one.")
                else:
                    gdf[store_shapecol_names[nr]] = mesh[variable][:].values
                
            shapefile_name = variable_name2 = os.path.splitext(os.path.basename(ugrid_ncdf))[0]
            path_shapefile_name = os.path.join(output_path, shapefile_name)
            store_shapefile_names.append(path_shapefile_name)
            gdf.to_file(path_shapefile_name)
            print("Exported "+path_shapefile_name)
        
        gdf_comp = pd.concat([
            gpd.read_file(shp)
            for shp in store_shapefile_names
        ]).pipe(gpd.GeoDataFrame)
        path_complete_shapefile_name = os.path.join(output_path, compiled_shapefile_name)
        gdf_comp.to_file(path_complete_shapefile_name)
        print("Exported "+path_complete_shapefile_name)

    def write_mesh(self, 
        fn_temp: str = "mesh/mesh_temp.nc",
        fn_ecoimpact: str = "mesh/ecoimpact_input.nc"
    ) -> None:
        """Write model grid data to netcdf file at <root>/<fn>

        key-word arguments are passed to :py:meth:`xarray.Dataset.ugrid.to_netcdf`

        Parameters
        ----------
        fn : str, optional
            filename relative to model root, by default 'grid/grid.nc'
        """
        if self._mesh2d is None:
            self.logger.debug("No mesh data found, skip writing.")
            return
        self._assert_write_mode
        # filename
        _fn_temp = join(self.root, fn_temp)
        _fn_ecoimpact = join(self.root, fn_ecoimpact)
        if not isdir(dirname(_fn_temp)):
            os.makedirs(dirname(_fn_temp))
        self.logger.debug(f"Writing file {fn_temp}")

        #write to location
        self._mesh2d.ugrid.to_netcdf(_fn_temp)

        #correct to useable UGRID consistent file
        ds_in = _fn_temp
        ds_out = _fn_ecoimpact
        self.correct_ugrid_file_for_ecoimpact(
            ds_in = ds_in, 
            ds_out = ds_out)
        self.logger.debug(f"Writing file {fn_ecoimpact}")

        #Write translation table categorical output