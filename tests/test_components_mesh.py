# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 02:38:07 2026

@author: veenstra
"""

import os
import numpy as np
from hydromt_delft3dfm import DFlowFMModel


def test_mesh_properties(tmpdir):
    root = os.path.join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=3857,
    )
    
    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=500,
    )
    
    assert mod1.crs.to_epsg() == 3857
    assert mod1.mesh.data.ugrid.crs["mesh2d"].to_epsg() == 3857
    assert mod1.mesh.crs == "EPSG:3857"
    
    mesh_bounds = mod1.mesh.bounds["mesh2d"]
    mesh_bounds_expected = np.array([
        1384000.0,
         5855500.0,
         1394000.0,
         5867000.0,
         ])
    assert np.allclose(mesh_bounds, mesh_bounds_expected)
    
    region_x, region_y = mod1.mesh._region_data.iloc[0].geometry.exterior.coords.xy
    region_x = np.array(region_x)
    region_y = np.array(region_y)
    assert np.isclose(region_x.min(), 1384000.0)
    assert np.isclose(region_x.max(), 1394000.0)
    assert np.isclose(region_y.min(), 5855500.0)
    assert np.isclose(region_y.max(), 5867000.0)
