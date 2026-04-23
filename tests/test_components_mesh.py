# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 02:38:07 2026

@author: veenstra
"""

import os
import numpy as np
import pytest
from hydromt_delft3dfm import DFlowFMModel


@pytest.mark.parametrize("crs", [3857, 4326])
def test_mesh_properties(tmpdir, crs):
    if crs == 3857:
        res = 500
        xmin, ymin, xmax, ymax = 1384000.0, 5855500.0, 1394000.0, 5867000.0
    elif crs == 4326:
        res = 0.05
        xmin, ymin, xmax, ymax = 12.45, 46.45, 12.50, 46.55
    
    root = os.path.join(tmpdir, "dflowfm_example")
    mod1 = DFlowFMModel(
        root=root,
        mode="w",
        crs=crs,
    )
    
    mod1.setup_mesh2d(
        region=dict(bbox=[12.4331, 46.4661, 12.5212, 46.5369]),
        res=res,
    )
    
    assert mod1.crs.to_epsg() == crs
    assert mod1.mesh.data.ugrid.crs["mesh2d"].to_epsg() == crs
    assert mod1.mesh.crs == f"EPSG:{crs}"
    
    mesh_bounds = mod1.mesh.bounds["mesh2d"]
    mesh_bounds_expected = np.array([xmin, ymin, xmax, ymax])
    assert np.allclose(mesh_bounds, mesh_bounds_expected)
    
    region_x, region_y = mod1.mesh._region_data.iloc[0].geometry.exterior.coords.xy
    region_x = np.array(region_x)
    region_y = np.array(region_y)
    assert np.isclose(region_x.min(), xmin)
    assert np.isclose(region_x.max(), xmax)
    assert np.isclose(region_y.min(), ymin)
    assert np.isclose(region_y.max(), ymax)
