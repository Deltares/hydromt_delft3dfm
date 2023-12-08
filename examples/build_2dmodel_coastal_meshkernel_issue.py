# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:29:44 2023

@author: veenstra
"""


import meshkernel
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import numpy as np

is_geographic = True
lon_min,lon_max = 11.82, 12.95
lat_min,lat_max = 45.20, 46.2
if is_geographic:
    crs = 'EPSG:4326'
    x_min, y_min = lon_min, lat_min
    x_max, y_max = lon_max, lat_max
    x_res, y_res = 0.1, 0.1
    projection = meshkernel.ProjectionType.SPHERICAL
else:
    crs = 'EPSG:3857'
    from pyproj import Transformer
    transformer = Transformer.from_crs('EPSG:4326', crs, always_xy=True)
    x_min, y_min = transformer.transform(lon_min, lat_min) # 1315992, 1442201
    x_max, y_max = transformer.transform(lon_max, lat_max) # 5653521, 5812289
    x_res, y_res = 10000, 10000
    projection = meshkernel.ProjectionType.CARTESIAN
    

# Create an instance of MakeGridParameters and set the values
make_grid_parameters = meshkernel.MakeGridParameters(angle=0.0, #TODO: does non-zero result in an orthogonal spherical grid?
                                                     origin_x=x_min,
                                                     origin_y=y_min,
                                                     upper_right_x=x_max, #TODO: angle and upper_right_x cannot be combined: https://github.com/Deltares/MeshKernelPy/issues/74
                                                     upper_right_y=y_max,
                                                     block_size_x=x_res,
                                                     block_size_y=y_res)

mk = meshkernel.MeshKernel(projection=projection)
mk.curvilinear_compute_rectangular_grid_on_extension(make_grid_parameters)
mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d

mesh2d_basegrid = mk.mesh2d_get() #in case of curvi grid: mk.curvilinear_convert_to_mesh2d()
fig, ax = plt.subplots()
mesh2d_basegrid.plot_edges(ax,linewidth=0.8)
dfmt.plot_coastlines(ax=ax, crs=crs)


#select bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

fig, ax = plt.subplots()
data_bathy_sel.elevation.plot()
dfmt.plot_coastlines(ax=ax, crs=4326)

#convert bathy data to GriddedSamples
lon_np = data_bathy_sel.lon.to_numpy()
lat_np = data_bathy_sel.lat.to_numpy()
values_np = data_bathy_sel.elevation.to_numpy().flatten().astype('float') #TODO: astype to avoid "TypeError: incompatible types, c_short_Array_74880 instance instead of LP_c_double instance"
gridded_samples = meshkernel.GriddedSamples(x_coordinates=lon_np,y_coordinates=lat_np,values=values_np) #TODO: does not result in refinement

#also to regular samples, since they have to be converted to 
lon_all_np, lat_all_np = np.meshgrid(lon_np, lat_np)
values_np_rav = values_np.ravel()
if is_geographic:
    x_all_np, y_all_np = lon_all_np, lat_all_np
else:
    x_all_np, y_all_np = transformer.transform(lon_all_np, lat_all_np)
samples = meshkernel.GeometryList(x_coordinates=x_all_np.ravel(),
                                  y_coordinates=y_all_np.ravel(),
                                  values=values_np_rav)

fig, ax = plt.subplots()
ax.scatter(x_all_np, y_all_np, c=values_np_rav)
dfmt.plot_coastlines(ax=ax, crs=crs)

#refinement
mesh_refinement_parameters = meshkernel.MeshRefinementParameters(max_refinement_iterations=2,
                                                                 min_edge_size=1000, #always in meters
                                                                 refinement_type=meshkernel.RefinementType.WAVE_COURANT,
                                                                 connect_hanging_nodes=True,
                                                                 smoothing_iterations=2,
                                                                 max_courant_time=120)
# griddedsamples not possible with coordinate conversion
# mk.mesh2d_refine_based_on_gridded_samples(gridded_samples=gridded_samples,
#                                           mesh_refinement_params=mesh_refinement_parameters,
#                                           use_nodal_refinement=True) #TODO: what does this do?
mk.mesh2d_refine_based_on_samples(samples=samples,
                                  relative_search_radius=1.1,
                                  minimum_num_samples=10,
                                  mesh_refinement_params=mesh_refinement_parameters,
                                  )

mesh2d_refinedgrid = mk.mesh2d_get()
fig, ax = plt.subplots()
mesh2d_refinedgrid.plot_edges(ax,linewidth=0.8)
dfmt.plot_coastlines(ax=ax, crs=crs)
