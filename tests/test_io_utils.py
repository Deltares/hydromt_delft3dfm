from hydromt_delft3dfm.utils.io_utils import read_branches_gui
import geopandas as gpd
from hydrolib.core.dflowfm.mdu import FMModel

def test_read_branches_gui():
    fmmodel = FMModel()
    fmmodel.filepath = "test2.mdu"
    gdf = gpd.GeoDataFrame()
    gdf["branchid"] = [1,2,3]
    gdf_out = read_branches_gui(gdf, fmmodel)

    assert gdf_out.shape == (3, 4)
    assert gdf_out["branchid"].tolist() == [1,2,3]
    assert (gdf_out["branchtype"] == "river").all()
    assert gdf_out[["manhole_up","manhole_dn"]].isnull().all().all()
