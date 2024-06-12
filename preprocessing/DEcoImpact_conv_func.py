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

class DEIModel(MeshModel):
    """This is the habitat model class"""

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
        self._categories = {}



    # COMPONENTS
    def setup_ugrid(
        self,
        region: dict,
        bounds: Optional[tuple] = None,
        name: str = None, 
        crs: Union["pyproj.CRS", str] = None,
    ):
        """
        Creates a 2D regular grid from an existing grid or a 2D unstructured mesh from a existing mesh.

        Adds/Updates model layers:
        * **grid** grid mask: add grid mask to grid object

        Parameters
        ----------
        region : dict
            Dictionary describing region of interest, e.g.:
            * {'ugrid': 'path/to/grid_file'}

            For habitat_grid model region must be of kind [grid].
        bounds : array-like of floats, optional
            (xmin, ymin, xmax, ymax) bounding box to clip an extract of the grid. Should be in the same crs as grid.
            By default None (no clipping).
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.    
        """
        kind = next(iter(region))  # first key of region
        if kind == "ugrid":
            mesh = self.read_ugrid(
                fn = region[kind],
                bounds = bounds,
                as_geom = True,
                crs = crs,
            )
        
        elif kind == "grid":
            mesh = self.read_grid(
                fn = region[kind],
                bounds = bounds,
                name = name,
                as_geom = True,
                crs = crs,
            )
            
        else:
            raise ValueError(
                "Only 'urid' or 'grid' region kind is supported for decoimpact_ugrid model"
            )

        # # If bounds, clip the grid to bounds
        # if bounds:
        #     mesh = grid.raster.clip_bbox(bbox=bounds, align=None, buffer=0)
        #     # Add grid mask to grid object
        #     self.set_mesh(mesh)
        #     # Update region and add to geoms
        #     geom = gpd.GeoDataFrame(
        #         geometry=[box(*grid.raster.bounds)], crs=grid.raster.crs
        #     )
        #     self.set_geoms(geom, name="region")

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

    def correct_mesh_connections(
        self,
        fn : Union[str, xr.DataArray, xr.Dataset],
        mesh : Union[xr.DataArray, xu.UgridDataset],
    ):
        """Method to correct mesh connections where necessary for face, nodes and edges.
        Parameters
        ----------
        mesh:
            A XUGRID UgridDataset
        """
        
        #TODO Hydromt-core functionality
        #clean up the mesh so all connectivity is present
        #add connectivity for faces

        if(not hasattr(mesh,"_coord_names")):       
            return(mesh)
        
        face_coord_present = False
        node_coord_present = False
        edge_coord_present = False
        

        face_coord_present = all(x in mesh._coord_names for x in ["mesh2d_face_x","mesh2d_face_y"])
        node_coord_present = all(x in mesh._coord_names for x in ["mesh2d_node_x","mesh2d_node_y"])
        edge_coord_present =  all(x in mesh._coord_names for x in ["mesh2d_edge_x","mesh2d_edge_y"])     
        
        faces_assessed = False
        nodes_assessed = False
        edges_assessed = False
        
        if(any(x not in mesh._coord_names for x in ["mesh2d_face_x","mesh2d_face_y","mesh2d_nFaces"])\
           and faces_assessed == False):     
            if(node_coord_present == True and face_coord_present == False):
                #try adding the edges
                try:
                    mesh =  mesh.ugrid.assign_edge_coords()
                    self.logger.debug(f"Assigning edge coords succesful {fn}")
                    edges_assessed = True
                except:
                    self.logger.debug(f"Assigning edge coords failed {fn}")
            #try adding the faces
            try:
                mesh = mesh.ugrid.assign_face_coords()
                self.logger.debug(f"Assigning face coords succesful {fn}")
                faces_assessed = True

                #rename nmesh2d_face
                if(x in list(mesh.coords.keys()) for x in ["nmesh2d_face"]):  # for sfincs
                    mesh = mesh.rename({
                            "nmesh2d_face":"mesh2d_nFaces",
                    })
            except:
                self.logger.debug(f"Assigning face coords failed {fn}")
        #if(any(x not in mesh._coord_names for x in ["mesh2d_node_x","mesh2d_node_y","mesh2d_nNodes"])\ # because of SFINX case
        if(any(x not in list(mesh.keys()) for x in ["mesh2d_node_x","mesh2d_node_y"])\
           and nodes_assessed == False):  
            #try adding the nodes
            try:
                mesh = mesh.ugrid.assign_node_coords()
                self.logger.debug(f"Assigning node coords succesful {fn}")
                nodes_assessed = True
            except:
                self.logger.debug(f"Assigning node coords failed {fn}")
        
        if(any(x not in mesh._coord_names for x in ["mesh2d_edge_x","mesh2d_edge_y","mesh2d_nEdges"])\
           and edges_assessed == False): 
            #try adding the edges
            try:
                mesh =  mesh.ugrid.assign_edge_coords()
                self.logger.debug(f"Assigning edge coords succesful {fn}")
                nodes_assessed = True
            except:
                self.logger.debug(f"Assigning edge coords failed {fn}")
        
        return(mesh)
    
    def correct_to_current_ugrid_format(
        self,
        fn : Union[str, xr.DataArray, xr.Dataset], 
    ):
        """Method to update the UGRID format when necessary.
        Parameters
        ----------
        mesh:
            A XUGRID UgridDataset
        """
        
        #TODO Hydromt-core functionality
        #clean up the mesh so all connectivity is present
        #add connectivity for faces
        fn_temp = "mesh/prep_mesh.nc"
        nc_temp = join(self.root, fn_temp)

        with nc.Dataset(fn) as src, nc.Dataset(nc_temp, "w") as dst:
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                if(name == "nNetElem"):
                     dst.createDimension(
                         "mesh2d_nFaces", (len(dimension) if not dimension.isunlimited() else None))
                # elif(name == "nFLowElem"):
                #     dst.createDimension(
                #         "mesh2d_nFaces", (len(dimension) if not dimension.isunlimited() else None))
                elif(name == "nNetElemMaxNode"):
                    dst.createDimension(
                        "mesh2d_nMax_face_nodes", (len(dimension) if not dimension.isunlimited() else None))
                # elif(name == "nFlowLinkPts"):
                #     dst.createDimension(
                #         "Two", (len(dimension) if not dimension.isunlimited() else None))
                elif(name == "nNetLinkPts"):
                    dst.createDimension(
                        "Two", (len(dimension) if not dimension.isunlimited() else None))
                elif(name == "nNetNode"):
                    dst.createDimension(
                        "mesh2d_nNodes", (len(dimension) if not dimension.isunlimited() else None))
                # elif(name == "nFlowLink"):
                #     dst.createDimension(
                #         "mesh2d_nEdges", (len(dimension) if not dimension.isunlimited() else None))
                elif(name == "nNetLink"):
                   dst.createDimension(
                       "mesh2d_nEdges", (len(dimension) if not dimension.isunlimited() else None))
                #elif(name == "nFlowElem" ): # already created with nNetElem
                #    pass
                #elif(name == "nNetLink" or name == "nNetLinkPts"): # already created with nFlowLink and nFlowLinkPts
                #    pass
                else:
                    dst.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))
                
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                
                #Based on old NetCDF UGRID standard translation
                if(name not in ["BndLink",\
                                #"FlowElemContour_x", "FlowElemContour_y",
                                "FlowElemDomain","FlowElemGlobalNr",\
                                "FlowLink_latu","FlowLink_lonu","FlowLinkDomain","FLowLinkType",\
                                "NetLink_xu","NetLink_yu","NetLinkContour_x","NetLinkContour_y",\
                                "NetLinkType","NetNode_z",\
                                #"NetNode_lat","NetNode_lon",
                                #"NetLink",\
                                "FLowElem_zcc","FlowLinkType"]):
                    if(len(variable.dimensions) > 0):
                        
                        new_dimensions = variable.dimensions
                        if("nNetElem" in variable.dimensions):
                            new_dimensions = [tuple(s if s != "nNetElem" else "mesh2d_nFaces" for s in new_dimensions)][0]
                        #if("nFlowElem" in variable.dimensions):
                        #    new_dimensions = [tuple(s if s != "nFlowElem" else "mesh2d_nFaces" for s in new_dimensions)][0]
                        #if("nFlowLink" in variable.dimensions):
                        #    new_dimensions = [tuple(s if s != "nFlowLink" else "mesh2d_nEdges" for s in new_dimensions)][0]
                        if("nNetElemMaxNode" in variable.dimensions):
                            new_dimensions = [tuple(s if s != "nNetElemMaxNode" else "mesh2d_nMax_face_nodes" for s in new_dimensions)][0]
                        if("nNetLinkPts" in variable.dimensions):
                            new_dimensions = [tuple(s if s != "nNetLinkPts" else "Two" for s in new_dimensions)][0]
                        #if("nFlowLinkPts" in variable.dimensions):
                        #    new_dimensions = [tuple(s if s != "nFlowLinkPts" else "Two" for s in new_dimensions)][0]
                        if("nNetNode" in variable.dimensions):
                            new_dimensions = [tuple(s if s != "nNetNode" else "mesh2d_nNodes" for s in new_dimensions)][0]
                        if("nNetLink" in variable.dimensions):
                            new_dimensions = [tuple(s if s != "nNetLink" else "mesh2d_nEdges" for s in new_dimensions)][0]
                    else:
                        new_dimensions = variable.dimensions

                    #rename variables
                    new_names_variables = {
                        "NetElemNode"          : "mesh2d_face_nodes", #old ugrid standard
                        "FlowLink"             : "mesh2d_edge_faces",
                        "NetLink"              : "mesh2d_edge_nodes",
                        "NetNode_x"            : "mesh2d_node_x",
                        "NetNode_y"            : "mesh2d_node_y",
                        "FlowElem_xcc"         : "mesh2d_face_x",
                        "FlowElem_ycc"         : "mesh2d_face_y",
                        "FlowLink_xu"          : "mesh2d_edge_x",
                        "FlowLink_yu"          : "mesh2d_edge_y",
                        "Mesh2D"               : "mesh2d",
                        #"nmesh2d_face"         : "mesh2d_nFaces",
                        #"nmesh2d_node"         : "nmesh2d_nNodes",
                        "FlowElem_bac"         : "mesh2d_flowelem_ba",
                        "FlowElem_bl"          : "mesh2d_flowelem_bl",
                        "mesh2d_edge_nodes_xu" : "mesh2d_edge_x",
                        "mesh2d_edge_nodes_yu" : "mesh2d_edge_y",
                        "FlowElemContour_x"    : "mesh2d_face_x_bnd",
                        "FlowElemContour_y"    : "mesh2d_face_y_bnd",
                    }
                    #create the variable
                    if name in list(new_names_variables.keys()):
                        x = dst.createVariable(new_names_variables[name], variable.datatype, new_dimensions)                        
                        new_name = new_names_variables[name]
                    else:
                        x = dst.createVariable(name, variable.datatype, new_dimensions)
                        new_name = name
                    
                    #copy variable attributes
                    if name == "Mesh2D":
                        mesh_dict = src[name].__dict__
                        mesh_dict = {k: mesh_dict[k] for k in set(list(mesh_dict.keys())) - set(["coordinates"])}
                        mesh_dict["cf_role"] = "mesh_topology"
                        mesh_dict["face_dimensions"] = "mesh2d_nFaces"
                        mesh_dict["node_dimensions"] = "mesh2d_nNodes"
                        mesh_dict["edge_dimensions"] = "mesh2d_nEdges"
                        mesh_dict["face_coordinates"] = "mesh2d_face_x mesh2d_face_y"
                        mesh_dict["node_coordinates"] = "mesh2d_node_x mesh2d_node_y"
                        mesh_dict["edge_coordinates"] = "mesh2d_edge_x mesh2d_edge_y"
                        mesh_dict["edge_node_connectivity"] = "mesh2d_edge_nodes"
                        mesh_dict["face_node_connectivity"] = "mesh2d_face_nodes"
                        mesh_dict["edge_face_connectivity"] = "mesh2d_edge_faces"
                        mesh_dict["max_face_nodes_dimensions"] = "mesh2d_nMax_face_nodes"
                        mesh_dict["name"] = "mesh2d"
                        mesh_dict["long_name"] = "Topology data of 2D mesh"

                        #   dst[name].edge_coordinates = "mesh2d_face_x mesh2d_face_y"       
                    
                    else:
                        # copy variable attributes all at once via dictionary
                        mesh_dict = src[name].__dict__
                        if("coordinates" in list(mesh_dict.keys())):
                            #replace coordinate names with new names
                            for word, initial in new_names_variables.items():
                                mesh_dict["coordinates"] = mesh_dict["coordinates"].replace(word,initial)
                        if("bounds" in list(mesh_dict.keys())):
                            #replace bounds names with new names
                            for word, initial in new_names_variables.items():
                                mesh_dict["bounds"] = mesh_dict["bounds"].replace(word,initial)
                                #mesh_dict = {k: mesh_dict[k] for k in set(list(mesh_dict.keys())) - set(["bounds"])}
                    
                    if name in list(new_names_variables.keys()):
                        if new_names_variables[name] == "mesh2d_face_nodes":
                            mesh_dict = src[name].__dict__
                            mesh_dict["cf_role"] = "face_node_connectivity"
                            mesh_dict["long_name"] = "Vertex nodes of mesh faces (counterclockwise)";
                            #if(src[name][:].min() < 0 and "_FillValue" not in list(mesh_dict.keys())):
                            if("_FillValue" not in list(mesh_dict.keys())):
                                dtype_str = src[name][:].dtype.str[1:]
                                mesh_dict["_FillValue"] = nc.default_fillvals[dtype_str]

                        elif new_names_variables[name] == "mesh2d_edge_faces":
                            mesh_dict = src[name].__dict__
                            mesh_dict["cf_role"] = "edge_face_connectivity"
                            mesh_dict["long_name"] = "Neighboring faces of mesh edges";
                            #mesh_dict["coordinates"] = "mesh2d_edge_y mesh2d_edge_x"
                            #if(src[name][:].min() < 0 and "_FillValue" not in list(mesh_dict.keys())):
                            if("_FillValue" not in list(mesh_dict.keys())):
                                dtype_str = src[name][:].dtype.str[1:]
                                mesh_dict["_FillValue"] = nc.default_fillvals[dtype_str]

                        elif new_names_variables[name] == "mesh2d_edge_nodes":
                            mesh_dict = src[name].__dict__
                            mesh_dict["cf_role"] = "edge_node_connectivity"
                            mesh_dict["long_name"] = "Start and end nodes of mesh edges";
                            #mesh_dict["coordinates"] = "mesh2d_edge_y mesh2d_edge_x"
                            if("_FillValue" not in list(mesh_dict.keys())):
                                dtype_str = src[name][:].dtype.str[1:]
                                mesh_dict["_FillValue"] = nc.default_fillvals[dtype_str]

                    #set attributes
                    dst[new_name].setncatts(mesh_dict)
                    
                    #add the data
                    # if(name == "NetElemNode"):
                    #     dst[new_name][:] = src[name][:]
                    #     dst[new_name][:][dst[new_name][:] < 0] = 0
                    #else:
                    dst[new_name][:] = src[name][:]                           

        return(nc_temp)    


    def correct_ugrid_file_for_ecoimpact(
        self,
        ds_in : Union[str, nc._netCDF4.Dataset],
        ds_out : Union[str, nc._netCDF4.Dataset],
    ):
        """Method to correct XUGRID created NetCDF files to the UGRID standard required
        by D-EcoImpact.
        Parameters
        ----------
        ds_in:
            Path to read a NetCDF UGRID file inconsistent with D-EcoImpact 
        ds_out:
            Path to write a NetCDF UGRID file corrected for D-EcoImpact
        """
        #TODO: temporary implementation outside of HYDROMT-core

        #add metadata
        with nc.Dataset(ds_in) as src, nc.Dataset(ds_out, "w") as dst:
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
                
                #midpoint between to edge nodes. Search node x position based on 
                # mesh2d_edge_nodes [0,1] and retrieve x coord from mesh2d_node_x (mesh2d_nEdges=variable)
                #dst.createDimension( "mesh2d_edge_x")
                #midpoint between to edge nodes. Search node y position based on 
                # mesh2d_edge_nodes [0,1] and retrieve y coord from mesh2d_node_x (mesh2d_nEdges=variable)
                #dst.createDimension( "mesh2d_edge_y")
                #midpoint of face. Search node x position based on 
                # mesh2d_face_nodes [0,1,2,3] and retrieve x coord from mesh2d_node_x.
                # Average all.(mesh2d_nFaces=variable)
                #dst.createDimension( "mesh2d_face_x")
                #midpoint of face. Search node y position based on 
                # mesh2d_face_nodes [0,1,2,3] and retrieve y coord from mesh2d_node_y.
                # Average all.(mesh2d_nFaces=variable)
                #dst.createDimension( "mesh2d_face_y")
                #coords of face nodes. Search node X position based on 
                # mesh2d_face_nodes [0,1,2,3] and retrieve all x coord from mesh2d_node_x.
                # Add them to a frame [0,1,2,3].(mesh2d_nFaces=variable, mesh2d_nMax_face_nodes=4)
                #dst.createDimension( "mesh2d_face_x_bnd")
                #coords of face nodes. Search node Y position based on 
                # mesh2d_face_nodes [0,1,2,3] and retrieve all Y coord from mesh2d_node_y.
                # Add them to a frame [0,1,2,3]. (mesh2d_nFaces=variable, mesh2d_nMax_face_nodes=4)
                #dst.createDimension( "mesh2d_face_y_bnd")

             # copy all file data except for the excluded
            for name, variable in src.variables.items():
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                
                #add variable attributes that were missing.
                if name == "mesh2d":
                    mesh_dict = src[name].__dict__
                    mesh_dict = {k: mesh_dict[k] for k in set(list(mesh_dict.keys())) - set(["coordinates",\
                            "face_face_connectivity" , "boundary_node_connectivity",\
                            "boundary_edge_dimension"])}
                    dst[name].setncatts(mesh_dict)
                    #add edge_coordinates attribute if missing
                    #if(not(hasattr(mesh_dict, 'edge_coordinates'))):
                    #    dst[name].edge_coordinates = "mesh2d_face_x mesh2d_face_y" 

                else:
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                
                if name not in ["mesh2d_edge_nodes","mesh2d_face_nodes","mesh2d_face_x","mesh2_face_y","mesh2d", "mesh2d_node_x",\
                        "mesh2d_node_y","time", "mesh2d_edge_faces","mesh2d_edge_nodes","mesh2d_edge_type","mesh2d_edge_x",\
                        "mesh2d_edge_y", "mesh2d_face_nodes", "mesh2d_face_x","mesh2d_face_x_bnd",\
                        "mesh2d_face_y","mesh2d_face_y_bnd",\
                        "dryFlagElement"]:
                    if("mesh2d_nFaces" in dst[name].dimensions):
                        if(not(hasattr(dst[name], 'mesh'))):
                            dst[name].mesh = "mesh2d"
                        if(not(hasattr(dst[name], 'location'))):
                            dst[name].location = "face"
                        if(not(hasattr(dst[name], 'grid_mapping'))):
                            dst[name].grid_mapping = "projected_coordinate_system"
                        if(not(hasattr(dst[name], 'coordinates'))):
                            dst[name].coordinates = "mesh2d_face_x mesh2d_face_y"
                    
                    elif("mesh2d_nNodes" in dst[name].dimensions):
                        if(not(hasattr(dst[name], 'mesh'))):
                            dst[name].mesh = "mesh2d"
                        if(not(hasattr(dst[name], 'location'))):
                            dst[name].location = "node"
                        if(not(hasattr(dst[name], 'grid_mapping'))):
                            dst[name].grid_mapping = "projected_coordinate_system"
                        if(not(hasattr(dst[name], 'coordinates'))):
                            dst[name].coordinates = "mesh2d_node_x mesh2d_node_y"                        

                    else:
                        pass
                else:
                    pass

                #add the data
                dst[name][:] = src[name][:]
            
            # #if mesh2d is missing
            # if("mesh2d" not in list(src.variables.keys())):
            #     x = dst.createVariable("mesh2d", np.int, ())
            #     dst["mesh2d"].cf_role = "mesh_topology"
            #     dst["mesh2d"].long_name = "Topology data of 2D mes"
            #     dst["mesh2d"].topology_dimension = int(2)
            #     dst["mesh2d"].node_coordinates = "mesh2d_node_x mesh2d_node_y"
            #     dst["mesh2d"].node_dimension = "mesh2d_nNodes"
            #     dst["mesh2d"].max_face_nodes_dimension = "mesh2d_nMax_face_nodes"
            #     dst["mesh2d"].edge_node_connectivity = "mesh2d_edge_nodes"
            #     dst["mesh2d"].edge_dimension = "mesh2d_nEdges"
            #     dst["mesh2d"].edge_coordinates = "mesh2d_edge_x mesh2d_edge_y"
            #     dst["mesh2d"].face_node_connectivity = "mesh2d_face_nodes"
            #     dst["mesh2d"].face_dimension = "mesh2d_nFaces"
            #     dst["mesh2d"].edge_face_connectivity = "mesh2d_edge_faces"
            #     dst["mesh2d"].face_coordinates = "mesh2d_face_x mesh2d_face_y"
            #     dst["mesh2d"].layer_dimension = "mesh2d_nLayers"
            #     dst["mesh2d"].interface_dimension = "mesh2d_nInterfaces"
            #     dst["mesh2d"].vertical_dimensions = "mesh2d_nLayers: mesh2d_nInterfaces (padding: none)"


    def read_ugrid(
        self,
        fn : Union[str, xr.DataArray, xr.Dataset],
        crs: Union["pyproj.CRS", str] = None,
        bounds: Optional[tuple] = None,
        as_geom : Optional[bool] = False,
    ):
        """Method to read a UGRID schematization from file.
        Parameters
        ----------
        fn:
            Dictionary describing region of interest.
            Currently supported format is {'ugrid': 'path/to/ugrid_file'} 
        bounds : array-like of floats, optional
            (xmin, ymin, xmax, ymax) bounding box to clip an extract of the grid. Should be in the same crs as grid.
            By default None (no clipping).
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.       
        as_geom : boolean
            choose this ugrid file as the geometry for the data to be mapped to
        """
        #TODO: temporary implementation outside of HYDROMT-core
        
        #check 
        
        try: 
            mesh = xu.open_dataset(
                fn,
            )
        except:
            #TODO: NOt correctly implemented yet how to handle old UGRID formats
            self.logger.debug("Trying to corrected ugrid format of file.")
            
            #correct old Ugrid format to new
            cor_mesh_fn = self.correct_to_current_ugrid_format(
                fn = fn,
            )
            mesh = xu.open_dataset(
                    cor_mesh_fn,
            )

            # try:
            #     mesh = xu.open_dataset(
            #         fn = cor_mesh_fn
            #     )
            # except:
            #     self.logger.debug("Not correctly implemented yet how to handle old UGRID formats or altering formats")
            #     return

        #correct mesh connections if necessary
        mesh = self.correct_mesh_connections(
                fn = fn, 
                mesh = mesh
        ) 

        #Select region of interest if applicable
        if(bounds != None):
            mesh = mesh.ugrid.sel(x=slice(bounds[0], bounds[2]), y=slice(bounds[1], bounds[3]))
        
        #Set loaded file as basis for mesh
        if(as_geom == True):
           self.set_as_base_mesh(mesh=mesh, crs=crs)
        
        return(mesh)


    def read_grid(self,
        fn: Union[str, xr.Dataset],
        bounds: Optional[tuple] = None,
        crs: Union["pyproj.CRS", str] = None,
        name: str = None,
        as_geom : Optional[bool] = False,
    ):
        """Method to read a GRID schematization from file.
        Parameters
        ----------
        region : dict
            Dictionary describing region of interest.
            Currently supported format is {'grid': 'path/to/grid_file'}
        bounds : array-like of floats, optional
            (xmin, ymin, xmax, ymax) bounding box to clip an extract of the grid. Should be in the same crs as grid.
            By default None (no clipping).
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.        
        as_geom : boolean
            choose this ugrid file as the geometry for the data to be mapped to
        """
        #TODO: temporary implementation outside of HYDROMT-core
        import xugrid as xu
        import rioxarray
        
        raster = rioxarray.open_rasterio(
                fn
        )
        
        #check if crs can be added
        if(crs != None):
            raster.rio.write_crs(crs, inplace=True)
            self.logger.debug(f"CRS manually added {raster.rio.crs}")
        else:
            if(hasattr(raster,"crs")):
                self.logger.debug(f"CRS used from file {raster.rio.crs}")
            else:
                self.logger.debug(f"No CRS supplied with file {fn}")

        #if extension is NetCDF
        var_list = []
        if(os.path.splitext(os.path.basename(fn))[1] == ".nc"):
            ignored_variables = set(["lat","layer","lon","time"])
            if(not isinstance(raster,xr.core.dataarray.DataArray)):
                bools = [not(x in ignored_variables) for x in list(raster.keys())]
                requested_variables = list(raster.keys())
                for cur_var in requested_variables:
                    var_list.append(xu.UgridDataArray.from_structured(raster[cur_var]))
                mesh = xu.merge(var_list, compat="identical")
            else:
                mesh = xu.merge([xu.UgridDataArray.from_structured(raster)], compat = "identical")
                print(mesh)
                mesh = mesh[next(iter(mesh))].fillna(-999.0)
        else:
            #convert to mesh
            dataarray = xu.UgridDataArray.from_structured(
                raster
            )
        
            #retrieve band, rename to file name and drop nulls
            if(name == None):
                variable_name = os.path.splitext(os.path.basename(fn))[0]
            else:
                variable_name = name
            
            dataarray = dataarray.sel(band=1).rename(
                    variable_name
                )

            var_list.append(dataarray)
            mesh = xu.merge(var_list, compat="identical")
            mesh = mesh.drop_vars("band")
            mesh = mesh.where(mesh.notnull(), drop = True)

        #correct mesh connections if necessary
        mesh = self.correct_mesh_connections(
                fn = fn, 
                mesh = mesh
                )

        #Select region of interest if applicable
        if(bounds != None):
            mesh = mesh.ugrid.sel(x=slice(bounds[0], bounds[2]), y=slice(bounds[1], bounds[3]))

        #If needed set loaded file as basis for mesh
        if(as_geom == True):
           self.set_as_base_mesh(mesh=mesh, crs=crs)           

        return(mesh)

    def read_vector(self,
        fn: Union[str, xr.Dataset],
        datacolumns: Union[list, str],
        bounds: Optional[tuple] = None,
        crs: Union["pyproj.CRS", str] = None,
        as_geom : Optional[bool] = False,
    ):
        """Method to read a VECTOR schematization from file.
        Parameters
        ----------
        region : dict
            Dictionary describing region of interest.
            Currently supported format is {'grid': 'path/to/grid_file'}
        bounds : array-like of floats, optional
            (xmin, ymin, xmax, ymax) bounding box to clip an extract of the grid. Should be in the same crs as grid.
            By default None (no clipping).
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.
        datacolumns : list  of strings
            The data columns that should be included in the UGRID        
        as_geom : boolean
            choose this ugrid file as the geometry for the data to be mapped to
        """
        #TODO: temporary implementation outside of HYDROMT-core
        import xugrid as xu
        import rioxarray
        
        vector = gpd.read_file(
                fn
        )
        
        #check if crs can be added
        if(crs != None):
            vector.set_crs(crs, inplace=True)
            self.logger.debug(f"CRS manually added {vector.crs}")
        else:
            if(hasattr(vector,"crs")):
                self.logger.debug(f"CRS used from file {vector.crs}")
            else:
                self.logger.debug(f"No CRS supplied with file {fn}")

        #extract requested data from shapefile
        var_list = []
        available_columns = list(vector.keys())
        if(all([x in available_columns for x in datacolumns])):
            keep_col = datacolumns + ['geometry']
            vector_sel = vector.drop(columns = vector.columns.difference(keep_col))
            vector_sel_single = vector_sel.explode()
            mesh = xu.UgridDataset.from_geodataframe(vector_sel_single)
            #for data_column in datacolumns:
            #    mesh[data_column] = vector_sel[data_column]

            #mesh = xu.merge(var_list, compat="identical")
            
        else:
            self.logger.debug("Some or all of the datacolumns metioned "+\
                              "do not exists in the vector file.")
            return

        #correct mesh connections if necessary
        mesh = self.correct_mesh_connections(
                fn = fn, 
                mesh = mesh
                )

        #Select region of interest if applicable
        if(bounds != None):
            mesh = mesh.ugrid.sel(x=slice(bounds[0], bounds[2]), y=slice(bounds[1], bounds[3]))

        #If needed set loaded file as basis for mesh
        if(as_geom == True):
           self.set_as_base_mesh(mesh=mesh, crs=crs)           

        return(mesh)

    def add_vector_data(self,
        fn: Union[str, xr.Dataset],
        datacolumns : Union[list, str],
        translation:  Union[str, dict],
        crs: Union["pyproj.CRS", str] = None,
    ):
        """Method to read a vector and impose the data on the present mesh.
        Parameters
        ----------
        fn : WindowsPath
            String describing the path to file.
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.        
        """
        #TODO: temporary implementation outside of HYDROMT-core
        #Code curtousy of Huite Bootsma

        #read the vector data
                #read the raster data
        
        vector = self.read_vector(
            fn = fn,
            crs = crs,
            datacolumns = datacolumns
        )
 
        if(self._mesh2d.ugrid.bounds == vector.ugrid.bounds):
            #same coordinate geometry
            self._mesh2d = xu.merge([self._mesh2d, vector], compat="identical")
        
        else:
            #read polygon file
            polygons_in = gpd.read_file(
                fn
                )
            if(polygons_in.crs == None):
                polygons_in.to_crs({'init' : crs})
        
            if translation.split("_")[0] == "larger":
                #match centroids of UGRID base file with polygon data                

                #Collect face centroids and join them to polygons
                centroids = gpd.GeoDataFrame(
                    geometry=gpd.points_from_xy(*self._mesh2d.ugrid.grid.face_coordinates.T),
                )
                group_joined = gpd.sjoin(centroids, polygons_in, predicate="within", how = "left")
                joined = group_joined.groupby('geometry').first()

                # Add the data to the UGRID data cube as a new variable
                for data_column in datacolumns:
                                      
                    if(joined[data_column].dtype == "string" or\
                        joined[data_column].dtype == "datetime" or\
                        joined[data_column].dtype == "object"):
                        #convert to string and assign categories
                        joined[data_column] = joined[data_column].astype(str)
                        joined[data_column] = joined[data_column].astype('category')
                        
                        #store categories
                        self._categories[data_column] = [[nr,cat]for nr, cat in\
                                         enumerate(joined[data_column].unique())]

                        #convert to integers
                        joined[data_column] = joined[data_column].cat.codes

                    self._mesh2d[data_column] = (self._mesh2d.ugrid.grid.face_dimension, joined[data_column])

            if translation.split("_")[0] == "smaller":
                #match centroids of polygon data with UGRID base file and average or choose most frequent
                polygons_ugrid = self._mesh2d.ugrid.grid.to_shapely(
                    dim = self._mesh2d.ugrid.grid.face_dimension
                )

                #Collect face centroids and join them to polygons
                centroids = gpd.GeoDataFrame(
                    geometry=polygons_in.centroid,
                )
                group_joined = gpd.sjoin(centroids, polygons_ugrid, predicate="within", how = "right")
                
                
                for data_column in datacolumns:
                    #Count number of occurances of value or characterset
                    if(translation.split("_")[1] == "count"):
                    
                        if(set(group_joined[data_column]) > 12):
                            self.logger.debug("More than 12 options available for max occurance count. "+\
                                            "Current limit is set to 12 unique values.")
                            return
                        
                        dfpivot = pd.pivot_table(group_joined,index='PolyID',columns=data_column,aggfunc={data_column:len})
                        dfpivot.columns = dfpivot.columns.droplevel()
                        cols_new = dfpivot.columns.difference(group_joined.columns, sort = False)
                        dfpivot[data_column] = dfpivot[cols_new].idmax(axis1)
                        joined = dfpivot.drop(cols_new, axis=1)
                    
                    #Get the statistics of the joined values           
                    elif(translation.split("_")[1] == "mean"):
                        joined = group_joined.groupby('id_right')[data_column].agg(['mean'])
                    elif(translation.split("_")[1] == "max"):
                        joined = group_joined.groupby('id_right')[data_column].agg(['max'])
                    elif(translation.split("_")[1] == "min"):
                        joined = group_joined.groupby('id_right')[data_column].agg(['min'])
                    elif(translation.split("_")[1] == "sum"):
                        joined = group_joined.groupby('id_right')[data_column].agg(['sum'])
                    else:
                        self.logger.debug("Method for 'smaller' is not defined. "+\
                                        "This should be 'count', 'mean', 'min', 'max' or 'sum'.")
                        return

                    # Add the data to the UGRID data cube as a new variable
                    self._mesh2d[data_column] = (self._mesh2d.ugrid.grid.face_dimension, joined[data_column])    

    def add_raster_data(self,
        fn: Union[str, xr.Dataset],
        translation:  Union[str, dict] = None,
        crs: Union["pyproj.CRS", str] = None,
        name: str = None,
    ):
        """Method to read a GRID schematization from file.
        Parameters
        ----------
        region : dict
            Dictionary describing region of interest.
            Currently supported format is {'grid': 'path/to/grid_file'}
        crs: pyproj.CRS
            The value can be anything accepted by pyproj.CRS.from_user_input(),
             such as an authority string (eg “EPSG:4326”) or a WKT string.        
        as_geom : boolean
            choose this ugrid file as the geometry for the data to be mapped to
        """
        #TODO: temporary implementation outside of HYDROMT-core
        import rasterstats as rs
        import rasterio as rio
        from rasterio.crs import CRS
        #Code curtousy of Huite Bootsma

        #read the raster data
        raster = self.read_grid(
            fn = fn,
            crs = crs,
            name = name,
        )
        
        if(self._mesh2d.ugrid.bounds == raster.ugrid.bounds and translation == None):
            #same coordinate geometry
            self._mesh2d = xu.merge([self._mesh2d, raster], compat="identical")
        else:
            #check if not netcdf file
            if(os.path.splitext(os.path.basename(fn))[1] == ".nc"):
                self.logger.debug("NetCDF raster data adding to a existing grid is not yet implemented.")
            
            #read raster file
            raster_in = rioxarray.open_rasterio(
                fn
            )

            if(raster_in.rio.crs == None):
                raster_in.rio.write_crs(crs, inplace=True)
                
            crs_mesh = next(iter(self._mesh2d.ugrid.crs.values()))
            crs_raster =  raster_in.rio.crs
            
            epsg_mesh = crs_mesh.to_epsg()
            epsg_raster = crs_raster.to_epsg()
            
            if(crs_mesh != crs_raster):
                raster_in=raster_in.rio.reproject(
                    'epsg:'+str(epsg_mesh), 
                    inplace = True
                    )
                
                #raster_in = rio.warp.transform_geom(
                #    CRS.from_epsg(epsg_raster),
                #    CRS.from_epsg(epsg_mesh),
                #    raster_in
                #    )

            if translation.split("_")[0] == "larger":
                self.logger.debug("Intersecting raster data based on mesh grid is not yet implemented.")

            if translation.split("_")[0] == "smaller":
                #use zonal statistics on raster data based on UGRID base file and choose statistic
                polygons_ugrid = self._mesh2d.ugrid.grid.to_shapely(
                    dim = self._mesh2d.ugrid.grid.face_dimension
                )
                polygons_gdf = gpd.GeoDataFrame(index = range(len(polygons_ugrid)), crs = crs, geometry=polygons_ugrid)

                #Apply zonal statistics to raster based on polgyons
                
                #Count number of occurances of cells in polygon
                statistic = None
                aggregation = None

                if(translation.split("_")[1] == "count"):
                    statistic = "count"      
                elif(translation.split("_")[1] == "mean"):
                    statistic = "mean"
                elif(translation.split("_")[1] == "max"):
                    statistic = "max"
                elif(translation.split("_")[1] == "min"):
                    statistic = "min"
                elif(translation.split("_")[1] == "sum"):
                    statistic = "sum"
                elif(translation.split("_")[1] == "fromsquarekm"):
                    statistic = "mean"
                    aggregation = "km2"
                elif(translation.split("_")[1] == "fromhectare"):
                    statistic = "mean"
                    aggregation = "hec"
                elif(translation.split("_")[1] == "fromsquaremeter"):
                    statistic = "mean"
                    aggregation = "m2"
                else:
                    self.logger.debug("Method for 'smaller' is not defined. "+\
                                        "This should be 'count', 'mean', 'min',"+\
                                      "'max', 'sum', 'fromsquarekm', 'fromhectare'"+\
                                      " or 'fromsquaremeter'.")
                    return()

                #Get zonal statistic required via rasterio.clip
                data=[]
                toti = len(polygons_gdf.index.unique()) 
                grid_result = None
                for i, a in enumerate(polygons_gdf.index.unique()):
                    
                    #print progress
                    if (i % 500) == 0:
                        print(" - processed " + str(i) + " of " + str(toti) + " features")
                    
                    grid_result = None
                    gdf_polygon = polygons_gdf[polygons_gdf.index == a]
                    try:
                        da_clip = raster_in.rio.set_spatial_dims(
                                x_dim="x", y_dim="y"
                            ).rio.clip(
                                gdf_polygon["geometry"], all_touched=True
                            )
                    except Exception as e:
                        if type(e)==rioxarray.exceptions.NoDataInBounds:
                            grid_result = np.nan
                            da_clip = None
                        else:
                            stop(e)

                    if(da_clip is not None):
                        if(statistic == "count"):
                            grid_result = da_clip.count(dim=["x", "y"],skipna=True).values[0]
                        elif(statistic == "mean"):
                            grid_result = da_clip.mean(dim=["x", "y"],skipna=True).values[0]
                        elif(statistic == "max"):
                            grid_result = da_clip.max(dim=["x", "y"],skipna=True).values[0]
                        elif(statistic == "min"):
                            grid_result = da_clip.min(dim=["x", "y"],skipna=True).values[0]
                        elif(statistic == "sum"):
                            grid_result = da_clip.sum(dim=["x", "y"],skipna=True).values[0]
                        else:
                            self.logger.error("Method for 'smaller' is not defined. "+\
                                            "This should be 'count', 'mean', 'min', "+\
                                            "'max' or 'sum' : " + str(statistic))
                            return()

                        #aggregate data
                        if(aggregation is not None):
                            if(aggregation == "km2"):
                                grid_result = grid_result * gdf_polygon["geometry"].area.values[0] / 10**6
                            elif(aggregation == "hec"):
                                grid_result = grid_result * gdf_polygon["geometry"].area.values[0] / 10**4
                            elif(aggregation == "m2"):
                                grid_result = grid_result * gdf_polygon["geometry"].area.values[0]
                            else:
                                self.logger.error("Method for aggregation is not defined. "+\
                                                "This should be 'km2', 'hec' or 'm2' : " + str(aggregation))
                                return()

                    #add data together
                    data.append(grid_result)
                    
                # # Get zonal statistic required via rasterstats
                # data_array = pd.DataFrame(rs.zonal_stats(polygons_gdf, raster_in.values,
                #     affine=raster_in.rio.transform(), stats=statistic, nodata=np.nan,
                #     all_touched=False))[statistic]

                # Add the data to the UGRID data cube as a new variable
 
                self._mesh2d[statistic +"_"+ name] = (self._mesh2d.ugrid.grid.face_dimension, data)  
        
        return(self._mesh2d)

    def drop_nans(
        self,
        variables : Union[list, str],
        nan_value : Union[float] = np.nan,
    ):
        if(isinstance(variables,str)):
            variables = [variables]

        if(isinstance(variables,list)):
            for variable in variables:
                if(nan_value == np.nan):
                    self._mesh2d = self._mesh2d.where(uda_ras1.notnull(), drop = True)
                else:
                    self._mesh2d = self._mesh2d.where(self._mesh2d[variable] != nan_value, drop = True)
        else:
            self.logger.debug(f"'variables' should be list or str {variables}")
        
        return()

    def drop_values(
        self,
        variables : Union[list, str],
        conditions : Union[str, list, float],
    ):
        if(isinstance(conditions,str) or isinstance(conditions,float)):
            conditions = [conditions]

        if(isinstance(variables,str)):
            variables = [variables]

        if(not isinstance(conditions,list)):
            self.logger.debug(f"'conditions' should be float, list or str {variables}")

        if(not isinstance(variables,list)):
            self.logger.debug(f"'variables' should be list or str {variables}")
        
        for condition, variable in zip(conditions,variables):
        
            if("<=" in condition):
                try: 
                    value = float(condition.strip("<="))
                except:
                    self.logger.debug(f"'condition' is not valid {condition}")
                    return()
                
                self._mesh2d = self._mesh2d.where(self._mesh2d[variable] <= value, drop = True)


            elif(">=" in condition):
                try: 
                    value = float(condition.strip(">="))
                except:
                    self.logger.debug(f"'condition' is not valid {condition}")
                    return()
                
                self._mesh2d = self._mesh2d.where(self._mesh2d[variable] >= value, drop = True)

            elif(">" in condition):
                try: 
                    value = float(condition.strip(">"))
                except:
                    self.logger.debug(f"'condition' is not valid {condition}")
                    return()
                
                self._mesh2d = self._mesh2d.where(self._mesh2d[variable] > value, drop = True)

            elif("<" in condition):
                try: 
                    value = float(condition.strip("<"))
                except:
                    self.logger.debug(f"'condition' is not valid {condition}")
                    return()
                
                self._mesh2d = self._mesh2d.where(self._mesh2d[variable] > value, drop = True)

            else:
                try: 
                    value = float(condition)
                except:
                    self.logger.debug(f"'condition' is not valid {condition}")
                    return()
                
                self._mesh2d = self._mesh2d.where(self._mesh2d[variable] != value, drop = True)

        return()

    def write_config(self,
        config_name: Optional[str] = None, 
        config_root: Optional[str] = None
    ):
        """Write config to <root/config_fn>"""
        self._assert_write_mode
        if config_name is not None:
            self._config_fn = config_name
        elif self._config_fn is None:
            self._config_fn = self._CONF
        if config_root is None:
            config_root = self.root
        fn = join(config_root, self._config_fn)
        self.logger.info(f"Writing model config to {fn}")
        self._configwrite(fn)

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