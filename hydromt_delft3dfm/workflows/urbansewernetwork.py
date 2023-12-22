"""Workflows to approximate urban sewer networks from OSM data."""

import logging
from pathlib import Path
from typing import Optional, Union

import networkx as nx
import numpy as np
import pandas as pd
import xarray as xr

# hydromt
from hydromt import DataCatalog
from scipy.optimize import curve_fit

from hydromt_delft3dfm import graph_utils, workflows

logger = logging.getLogger(__name__)

__all__ = [
    "setup_urban_sewer_network_topology_from_osm",
    "setup_urban_sewer_network_bedlevel_from_dem",
    "setup_urban_sewer_network_runpoffarea_from_landuse",
    "setup_rainfall_function_from_stats",
    # calculation workflows
    "select_connected_branches",
    "optimise_pipe_topology",
    "calculate_hydraulic_parameters",
]


def setup_urban_sewer_network_topology_from_osm(
    region,
    highway_types: list[str] = None,
    waterway_types: list[str] = None,
    simplify: bool = True,
    use_connected_only: bool = True,
    logger: logging.Logger = logger,
) -> nx.MultiDiGraph:
    """
    Set up a complete urban sewer network topology OpenStreetMap.

    Include data for waterways and highways for open and closed systems.
    The two systems are connected and outlets are added at intersection points.

    Adds/Updates model layers:

    * **rivers** geom: 1D rivers vector
    * **pipes** geom: 1D pipes vector
    * **branches** geom: 1D branches vector
    * **outlets** geom: 1D outlets vector

    Parameters
    ----------
    region : object
        The geographical region for which the sewer network is to be constructed.
        This could be a bounding box, polygon, or any other representation of a region.
    highway_types : list
        List of highway types (e.g., ["motorway", "primary"]) to include
        from the OpenStreetMap data.
    waterway_types : list
        List of waterway types (e.g., ["river", "stream"]) to include
        from the OpenStreetMap data.
    simplify: bool
        Simplify a graph’s topology by removing interstitial nodes.
        Cautious use because might leave fragmented edges.
        By default True.
    use_connected_only: bool
        Use largest connected components in the network only.
        by default True.
    logger : Logger object
        An instance of a logger to capture logs and messages during the process.
        Useful for debugging and understanding the workflow steps.

    Returns
    -------
    - graph : nx.MultiDiGraph
        The final preprocessed graph representing the urban sewer network.
        It integrates both the open and closed drainage systems and includes outlets at
        the intersection points.
    """
    # input
    if waterway_types is None:
        waterway_types = ["river", "stream", "brook", "canal", "ditch"]
    if highway_types is None:
        highway_types = [
            "motorway",
            "motorway_link",
            "primary",
            "primary_link",
            "secondary",
            "secondary_link",
            "tertiary",
            "tertiary_link",
            "residential",
        ]
    # output
    _required_columns = [
        "id",
        "branchid",
        "node_start",
        "node_end",
        "branchtype",
        "geometry",
    ]

    # 1. Build the graph for waterways (open system)
    logger.info(f"Download waterway {waterway_types} form osm")
    graph_osm_open = workflows.create_graph_from_openstreetmap(
        region=region,
        buffer=1000,
        osm_key="waterway",
        osm_values=waterway_types,
        simplify=simplify,
        logger=logger,
    )

    # 2. Build the graph for highways (closed system)
    logger.info(f"Download highway {highway_types} form osm")
    graph_osm_closed = workflows.create_graph_from_openstreetmap(
        region=region,
        buffer=500,
        osm_key="highway",
        osm_values=highway_types,
        simplify=simplify,
        logger=logger,
    )

    # 3. obtain open and closed systems and outlets by intersecting
    open_system, closed_system, outlets = workflows.intersect_lines(
        graph_utils.graph_to_network(graph_osm_open)[0],
        graph_utils.graph_to_network(graph_osm_closed)[0],
    )
    # add attributes
    open_system["branchid"] = open_system["id"]
    open_system["branchtype"] = "river"
    closed_system["branchid"] = closed_system["id"]
    closed_system["branchtype"] = "pipe"
    outlets["nodeid"] = outlets["id"]
    outlets["nodetype"] = "outlet"
    # trim columns
    branches = pd.concat([open_system, closed_system])[_required_columns]

    # 4. Recreate the graph from the complete system topology
    logger.info("Create graph and update open system, closed system and outlets.")
    graph = workflows.create_graph_from_geodataframe(
        branches,
        is_directed=False,
        logger=logger,
    )

    # TODO set graph

    # 5. Add outlets to graph
    # TODO replaced by super method
    graph = workflows.setup_graph_from_geodataframe(
        graph,  # TODO replace by self.graph
        data_catalog=None,  # TODO replace by self.data_catalog
        vector_fn=outlets,
        variables=["nodeid", "nodetype"],
        max_dist=1.0,
        rename=dict(nodeid="nodeid", nodetype="nodetype"),
        graph_component="nodes",
        logger=logger,
    )

    if use_connected_only is True:
        graph = filter_connected_branches(graph, ["river", "pipe"])

    return graph


def filter_connected_branches(
    graph: nx.Graph, branchtypes: list[str] = ["river", "pipe"]
):
    """Set up the urban sewer network based on connected componants."""
    graph_out = nx.Graph()
    for branchtype in branchtypes:
        _graph = select_connected_branches(graph, branchtype)
        graph_out = nx.compose(graph_out, _graph)
    return graph_out


def setup_urban_sewer_network_bedlevel_from_dem(
    graph: nx.Graph,  # TODO replace by self.graphs
    data_catalog: DataCatalog,  # TODO replace by self.data_catalog
    dem_fn: Union[str, Path, xr.DataArray, xr.Dataset],
    fill_method: Optional[str] = None,
    logger: logging.Logger = logger,
):
    """
    Update the graph with elevation data from a DEM.

    Add elevation to nodes, compute gradient for edges,
    and reverse the direction of edges with negative gradient

    The function samples elevation data from the provided DEM to the nodes of the graph.
    It then computes the gradient of edges based on the elevation difference between
    the nodes and the edge length. If the computed gradient is negative, the direction
    of the edge is reversed.

    Parameters
    ----------
    graph : nx.Graph
        Input graph to be updated.
    data_catalog : DataCatalog
        Data catalog to read and process raster data.
    dem_fn : Union[str, Path, xr.DataArray, xr.Dataset]
        Data catalog key, path to DEM file or DEM xarray data object.
    fill_method : str, optional
        If specified, fills no data values using fill_nodata method.
        Available methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
    logger : logging.Logger, optional
        Logger object to log messages. Default is the global logger.

    Returns
    -------
    nx.Graph
        Graph with updated ["elevtn"] (nodes), ["gradient"] (edges), and potentially
        reversed edge directions.
    """
    logger.info(
        "Update the graph with elevation data from a DEM and compute gradient for edges"
    )

    # Sample DEM
    # TODO replaced by super method
    graph = workflows.setup_graph_from_rasterdataset(
        graph=graph,
        data_catalog=data_catalog,
        raster_fn=dem_fn,
        fill_method=fill_method,
        graph_component="nodes",
        logger=logger,
    )
    # TODO set graph

    return graph


def setup_urban_sewer_network_runpoffarea_from_landuse(
    graph: nx.Graph,  # TODO replace by self.graphs
    data_catalog: DataCatalog,  # TODO replace by self.data_catalog
    landuse_fn: Union[str, Path, xr.DataArray, xr.Dataset],
    landuse_reclass_fn: Union[str, Path, pd.DataFrame],
    runoff_distance: float,
    fill_method: Optional[str] = None,
    logger: logging.Logger = logger,
):
    """
    Update the graph with landuse data.

    Add runoff area ("runoff_area") edges.

    The function samples classifies raster data in `landuse_fn` based on
    "impervious_ratio" in `landuse_reclass_fn`.
    The runoff area is then calculated as the total impervious area within
    the 'runoff_distance' from the network.

    Parameters
    ----------
    graph : nx.Graph
        Input graph to be updated.
    data_catalog : DataCatalog
        Data catalog to read and process raster data.
    landuse_fn : Union[str, Path, xr.DataArray, xr.Dataset]
        Data catalog key, path to landuse file or landuse xarray data object.
    landuse_reclass_fn : Union[str, Path, pd.DataFrame]
        Data catalog key, path to landuse reclass file or landuse reclass dataframe.
        Must contain convertion into "impervious_ratio".
    runoff_distance: float
        Distance within which the runoff will be routed to the network. (m)
    fill_method : str, optional
        If specified, fills no data values using fill_nodata method.
        Available methods are {'linear', 'nearest', 'cubic', 'rio_idw'}.
    logger : logging.Logger, optional
        Logger object to log messages. Default is the global logger.

    Returns
    -------
    nx.Graph
        Graph with updated ["runoff_area"] (edges)
    """
    logger.info(
        "Update the graph with elevation data from a DEM and compute gradient for edges"
    )

    # Sample landuse for mean impervious_ratio within runoff_distance
    # TODO replaced by super method
    graph = workflows.setup_graph_from_raster_reclass(
        graph=graph,
        raster_fn=landuse_fn,
        reclass_table_fn=landuse_reclass_fn,
        reclass_variables=["impervious_ratio"],
        fill_method=fill_method,
        resampling_method="mean",
        graph_component="edges",
        edge_buffer=runoff_distance,
        data_catalog=data_catalog,
        logger=logger,
    )
    # compute runoff area
    runoff_area = {
        (edge_start, edge_end): edge_data["geometry"].length
        * runoff_distance
        * 2
        * edge_data["impervious_ratio"]
        for edge_start, edge_end, edge_data in graph.edges(data=True)
    }

    nx.set_edge_attributes(graph, runoff_area, name="runoff_area")

    return graph


def select_connected_branches(graph, branchtype):
    """Select connected branches based on branchtype."""
    # Truncate to connected branches only
    graph_branchtype = workflows.query_graph_edges_attributes(
        graph, edge_query=f'branchtype=="{branchtype}"'
    )
    graph_branchtype = workflows.get_largest_component(graph_branchtype)

    # Add missing outlet edges
    graph_branchtype = workflows.add_missing_edges_to_subgraph(graph_branchtype, graph)

    return graph_branchtype


def _select_pipes(graph):
    """Select only connected pipe from graph."""
    # Truncate to connected pipes only
    graph_pipes = workflows.query_graph_edges_attributes(
        graph, edge_query='branchtype=="pipe"'
    )
    graph_pipes = workflows.get_largest_component(graph_pipes)

    # Add missing outlet pipes
    graph_pipes = workflows.add_missing_edges_to_subgraph(graph_pipes, graph)
    return graph_pipes


def optimise_pipe_topology(
    graph: nx.Graph, method_for_weight: str = "length", logger=logger
):
    """
    Process the graph by optimising the pipe network topology.

    Include steps to truncate it to connected pipes only, compute a DAG,
    and update the graph using input graph network informations.

    Parameters
    ----------
    graph: nx.Graph
        Graph that contains pipes that needs optimisation.
        Nodes must have attribute "nodetype"; edges must have attribute "branchtype".
    method_for_weight: str
        Choose the method for weight calculation, e.g. "length", "static_loss".
        If "static_loss" is used, the node must have attribute "elevtn".

    Returns
    -------
        nx.Graph: Processed graph in the form of a DAG.

    """
    # FIXME missing edges in between adding them back result in a nondag

    logger.info("Select only connected pipe from graph")
    graph_pipes = _select_pipes(graph)

    logger.info(f"Optimize based on {method_for_weight} method for weight.")
    if method_for_weight == "length":
        graph_pipes_optmised = _optimize_pipe_topology_based_on_length_as_weight(
            graph_pipes
        )
    elif method_for_weight == "static_loss":
        graph_pipes_optmised = _optimize_pipe_topology_based_on_static_loss_as_weight(
            graph_pipes
        )
    else:
        logger.error("method not supported. Use from : length, static_loss)")

    logger.info("compute gradient for pipes")
    graph_pipes_optmised = _compute_gradient_from_elevtn(graph_pipes_optmised)
    return graph_pipes_optmised


def _optimize_pipe_topology_based_on_length_as_weight(graph_pipes):
    """Optimize graph based on length as weight."""
    # dag
    _dag = workflows.setup_dag(
        graph_pipes,
        target_query='nodetype=="outlet"',
        weight="length",  # report="1"
    )
    # set output graph using dag edge directions
    graph_pipes_optmised = workflows.set_directions(graph_pipes.to_undirected(), _dag)

    # recompute geometry based on new edges directions
    graph_pipes_optmised = workflows.update_edge_geometry(graph_pipes_optmised)
    return graph_pipes_optmised


def _optimize_pipe_topology_based_on_static_loss_as_weight(graph_pipes):
    """Optimize graph based on static loss as weight.

    if only elevtn is available at nodes, then the cost along path is used
    need to determin which nodes are sources, which nodes are targets

    Note: not suited for flat area
    """
    # add static loss (based on elevation of nodes) as weight on both directions
    for e in graph_pipes.edges(data=True):
        e[2].update(
            {
                "weight": graph_pipes.nodes[e[1]]["elevtn"]
                - graph_pipes.nodes[e[0]]["elevtn"]
            }
        )  # if upslope -> positive loss, if downslope -> negative loss

    # dag
    _graph = workflows.preprocess_dag(graph_pipes, weight="weight")
    _dag = workflows.setup_dag(
        _graph,
        target_query='nodetype=="outlet"',
        weight="weight",
        algorithm="flowpath",
    )
    _dag = workflows.postprocess_dag(_dag)

    # set output graph using dag edge directions
    graph_pipes_optmised = workflows.set_directions(graph_pipes.to_undirected(), _dag)

    # recompute geometry based on new edges directions
    graph_pipes_optmised = workflows.update_edge_geometry(graph_pipes_optmised)
    return graph_pipes_optmised


def _optimize_pipe_topology_based_on_network_simplex(graph_pipes):
    # G = digraph_to_multidigraph(digraph)
    # H = multidigraph_to_digraph_by_adding_node_at_mid_location(G)
    # # add super target
    # H.add_node(
    #     "T", demand=-1 * sum([n[1] for n in H.nodes(data="demand")])
    # )  # demand must equals to all supply
    # H.add_edge("B", "T")  # Edge from B to T with capacity of 3 and cost of 0
    # H.add_edge("C", "T")  # Edge from B to T with capacity of 3 and cost of 0
    # # Find the minimum cost flow using the network simplex algorithm
    # flow_cost, flow_dict = nx.network_simplex(H)
    # # Print the minimum cost
    # print("The minimum cost of the flow is", flow_cost)
    # # Print the flow values of each edge
    # for u, v in H.edges():
    #     print("The flow on the edge from", u, "to", v, "is", flow_dict[u][v])

    # # Create a list of edge colors based on the flow values
    # edge_colors = []
    # for u, v in H.edges():
    #     if flow_dict[u][v] > 0:
    #         edge_colors.append("green")  # Green color for edges with positive flow
    #     else:
    #         edge_colors.append("white")  # Black color for edges with zero flow

    # # Draw the graph with the nodes, edges, and labels
    # nx.draw(H, with_labels=True, node_color="white", edge_color=edge_colors)

    # # Show the plot
    # plt.show()
    return None


def _compute_gradient_from_elevtn(graph):
    """
    Compute the gradient on edge based on node elevations and edge length.

    Gradient is calculated as:
    gradient = (elevation of source node - elevation of target node) / edge length

    Modifies the input graph in place.
    """
    # TODO: check if need to convert into digraph first

    for e in graph.edges:
        graph.edges[e].update(
            {
                "gradient": (
                    (graph.nodes[e[0]]["elevtn"] - graph.nodes[e[1]]["elevtn"])
                    / graph.edges[e]["length"]
                )
            }
        )
    return graph


def calculate_hydraulic_parameters(
    pipe_network_graph: nx.DiGraph,
    flow_velocity: float = 1,
    pipe_depth: float = 0.5,
    rainfall_depth_function=lambda x: 20,
    rounding_precision: int = 1,
    logger=logger,
):
    """
    Calculate radius for a pipe network.

    Consider various geographic and design assumptions.

    Parameters
    ----------
    pipe_network_graph: networkx.DiGraph
        A directed graph representing the pipe network, where edges have 'length' and
        'runoff_area' attributes.
    flow_velocity : float
        Average flow velocity in the pipes during near-full conditions (m/s).
        This is an empirical value and may vary by region.
        By default 1 m/s is used.
    pipe_depth : float
        Average depth of pipe underground (m).
        Depth is deinfed as the elevation different between the surface and the top of
        the pipe.
        This is an empirical value and may vary by region.
        By default 0.5 m is used.
    rainfall_depth_function : function
        A function that takes concentration time as input and returns rainfall depth
        (mm).
        This can be a constant value for simpler models,
        or a more complex function reflecting design criteria which may vary by region.
        By default, a constant value of 20 mm is returned.
    rounding_precision: int
        Number of decimal places to round to for the estimated dimentions.
        Applies to "diameter", "invlev_up", "invlev_dn".
        By default 1.

    Returns
    -------
    networkx.DiGraph:
        The graph with updated attributes.
        Each edge contains ["runoff_area", "diameter", "invlev_up", "invlev_dn"]
    """
    if not nx.is_directed_acyclic_graph(pipe_network_graph):
        raise ValueError("The input graph must be a directed acyclic graph (DAG)")

    # Calculate pipe radius in topological order
    pipe_network_graph = _calculate_pipe_radius(
        pipe_network_graph, flow_velocity, rainfall_depth_function, logger
    )

    # calculate other hydraulic dimentions
    for s, t in pipe_network_graph.edges():
        # calculate other dimentions
        radius = pipe_network_graph[s][t]["radius"]
        current_pipe_diameter = radius * 2
        invertlevel_up = (
            pipe_network_graph.nodes[s]["elevtn"] - current_pipe_diameter
        ) - pipe_depth
        invertlevel_down = (
            pipe_network_graph.nodes[t]["elevtn"] - current_pipe_diameter
        ) - pipe_depth
        nx.set_edge_attributes(
            pipe_network_graph,
            {
                (s, t): {
                    "diameter": round(current_pipe_diameter, rounding_precision),
                    "invlev_up": round(invertlevel_up, rounding_precision),
                    "invlev_dn": round(invertlevel_down, rounding_precision),
                }
            },
        )

    # FIXME inverstigate why there are still huge pipes ('2889', '4637', 8.3)
    return pipe_network_graph


def __test_calculate_pipe_radius():
    # Define a simple rainfall depth function for testing
    def simple_rainfall_depth_function(concentration_time):
        return 20 * concentration_time

    # Create a small directed graph
    pipe_network_graph = nx.DiGraph()

    # Add edges with 'length' and 'runoff_area' attributes
    pipe_network_graph.add_edge("A", "B", length=100, runoff_area=100)
    pipe_network_graph.add_edge("B", "C", length=100, runoff_area=100)
    pipe_network_graph.add_edge("A", "D", length=100, runoff_area=100)
    pipe_network_graph.add_edge("D", "C", length=100, runoff_area=100)
    pipe_network_graph.add_edge("C", "E", length=100, runoff_area=100)

    # Set the initial radius of all edges to zero
    nx.set_edge_attributes(pipe_network_graph, 0, "radius")

    # Apply the function
    pipe_network_graph = _calculate_pipe_radius(
        pipe_network_graph, 1, simple_rainfall_depth_function
    )

    # Check the results
    for u, v in pipe_network_graph.edges():
        print(f"Edge ({u}, {v}) has radius: {pipe_network_graph[u][v]['radius']}")

    # Visualize the graph
    import matplotlib.pyplot as plt

    nx.draw(
        pipe_network_graph,
        with_labels=True,
        node_color="lightblue",
        node_size=2000,
        font_size=20,
        font_weight="bold",
    )
    plt.show()


def _dfs_longest_path_properties(
    node, visited, reversed_graph, current_path_length=0, longest_path=0
):
    """
    Perform a DFS traversal to find the longest upstream path for a given node.

    Parameters
    ----------
    node: The current node being visited.
    visited: A set of nodes that have already been visited.
    reversed_graph: The graph used for traversal.
    current_path_length: The length of the current path in the traversal.
    longest_path: The length of the longest path found so far.

    Returns
    -------
    The length of the longest upstream path.
    """
    if node in visited:
        return longest_path

    visited.add(node)
    for start_node, target_node in reversed_graph.out_edges(node):
        if target_node in visited:
            continue

        edge_length = reversed_graph[start_node][target_node]["length"]
        new_path_length = current_path_length + edge_length

        longest_path = max(
            longest_path,
            _dfs_longest_path_properties(
                target_node,
                visited.copy(),
                reversed_graph,
                new_path_length,
                longest_path,
            ),
            new_path_length,
        )

    return longest_path


def _dfs_accumulate_properties(
    node,
    visited,
    reversed_graph,
    total_upstream_runoff_area=0,
    total_upstream_volume=0,
):
    if node in visited:
        return (
            total_upstream_runoff_area,
            total_upstream_volume,
        )

    visited.add(node)
    for start_node, target_node in reversed_graph.out_edges(node):
        if target_node in visited:
            continue
        edge_runoff_area = reversed_graph[start_node][target_node]["runoff_area"]
        edge_length = reversed_graph[start_node][target_node]["length"]
        edge_radius = reversed_graph[start_node][target_node]["radius"]
        edge_volume = np.pi * (edge_radius**2) * edge_length

        (
            child_runoff_area,
            child_volume,
        ) = _dfs_accumulate_properties(
            target_node,
            visited.copy(),
            reversed_graph,
            total_upstream_runoff_area + edge_runoff_area,
            total_upstream_volume + edge_volume,
        )

        # Update the totals with the maximum values from all branches
        total_upstream_runoff_area = max(total_upstream_runoff_area, child_runoff_area)
        total_upstream_volume = max(total_upstream_volume, child_volume)

    return (
        total_upstream_runoff_area,
        total_upstream_volume,
    )


def _calculate_pipe_radius(
    pipe_network_graph, flow_velocity, rainfall_depth_function, logger=logger
):
    """
    Calculate the radius of pipes  based on hydraulic calculations.

    This function iterates through the nodes in topological order.
    For each edge, it performs hydraulic calculations to determine the required average
    radius to accommodate the upstream volume.

    Parameters
    ----------
    pipe_network_graph: networkx.DiGraph
        A directed graph representing the pipe network.
    flow_velocity: float
        The flow velocity used in hydraulic calculations.
    rainfall_depth_function: function
        A function that takes concentration time upstream
        and returns the rainfall depth.

    Returns
    -------
    networkx.DiGraph with new properties (radius)

    See Also
    --------
    _dfs_accumulate_properties
    _dfs_longest_path_properties
    """
    nx.set_edge_attributes(pipe_network_graph, 0, "radius")

    reversed_graph = pipe_network_graph.reverse()

    for current_node in nx.topological_sort(pipe_network_graph):
        for start_node, target_node in pipe_network_graph.out_edges(current_node):
            # get accumulated properties upstream
            (
                upstream_runoff_area,
                upstream_volume,
            ) = _dfs_accumulate_properties(start_node, set(), reversed_graph)
            # get longest path properties upstream
            upstream_longestpath = _dfs_longest_path_properties(
                start_node, set(), reversed_graph
            )
            # get downstream edge properties (volume unknown and to be solved)
            runoff_area = reversed_graph[target_node][start_node]["runoff_area"]
            length = reversed_graph[target_node][start_node]["length"]
            # get total properties
            total_runoff_area = upstream_runoff_area + runoff_area
            total_longestpath = upstream_longestpath + length
            concentration_time = total_longestpath / flow_velocity / 3600.0
            rainfall_depth = rainfall_depth_function(concentration_time)
            total_volume = rainfall_depth / 1000.0 * total_runoff_area
            # solve volume and radius
            volume = total_volume - upstream_volume
            radius = np.sqrt(volume / (np.pi * length)) if length > 0 else 0
            # get immediate upstream radius to compare --> no sharp increase
            immediate_upstream_radius = sum(
                [e[2] for e in reversed_graph.out_edges(current_node, data="radius")]
            )
            # assign radius
            reversed_graph[target_node][start_node]["radius"] = radius
            # quality check
            # TODO: there are still very large pipes
            if radius > immediate_upstream_radius and immediate_upstream_radius > 0:
                logger.debug(
                    f"radius {radius} exceeds limits(sum of upstream radius): "
                    + f"Fall back to maximum allowed {immediate_upstream_radius}"
                )
                reversed_graph[target_node][start_node][
                    "radius"
                ] = immediate_upstream_radius

    return reversed_graph.reverse()


def _get_rainfall_from_gpex(
    ds: xr.Dataset,
    rt: int = 2,
):
    """
    Calculate the rainfall depth function from gpex data.

    Temporal downscaling is done by exponential model fitted to
    GEV estimates for a given return period.

    Parameters
    ----------
    ds: xarray.Dataset
        Contain 'gev_estimates' data variable with coodinates 'dur' and 'rt'.
    rt: int
        Return period for which to calculate the rainfall depth (in years).
        Must be among [2 5 10 20 39 50 100 200 500 1000]
        By default 2.

    Returns
    -------
    np.array
        Rainfall intensity as an exponential function
    """
    # Check if 'gev_estimate' in the dataset
    if "gev_estimate" not in ds:
        raise ValueError("Dataset does not contain the required data.")
    durations = ds.coords["dur"].values
    rainfall_rate = ds["gev_estimate"].sel(tr=rt).values / durations

    # Fitting the extremes using the specified distribution and duration
    def exponential_model(duration, rainfall_rate, alpha):
        return rainfall_rate * (duration**alpha)

    # Curve fitting
    params, params_covariance = curve_fit(
        exponential_model, durations, rainfall_rate, p0=[1, 0.5]
    )
    rainfall_rate_max_fitted, alpha_fitted = params

    # plt.figure(figsize=(8, 5))
    # plt.scatter(durations, rainfall_rate, label='Provided Data')
    # plt.plot(durations, exponential_model(durations,
    # rainfall_rate_max_fitted,alpha_fitted),
    # label='Fitted Curve', color='red')
    # # Setting the scale of both axes to logarithmic
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('Duration (hours) [Log Scale]')
    # plt.ylabel('Rainfall rate (mm/hr) [Log Scale]')
    # plt.title('Rainfall rate vs Duration with Fitted Curve (Log-Log Plot)')
    # plt.legend()
    # plt.grid(True, which="both", ls="--")  # Grid lines for log scale
    # plt.show()

    def rainfall_depth_function(tc):
        # compute intensity for tc
        I = exponential_model(tc, rainfall_rate_max_fitted, alpha_fitted)
        # convert to depth
        return I * tc

    return rainfall_depth_function


def setup_rainfall_function_from_stats(
    rainfall_fn=None,
    region=None,
    data_catalog=None,
    assumption_rainfall_intensity: float = 40,
):
    """
    Set up a function to calculate rainfall depth based on concentration time.

    Using rainfall statistics derived from historical data (gpex only) or
    a default assumption of rainfall intensity.

    Parameters
    ----------
    rainfall_fn: Union[str, Path], optional
        The file path to the historical rainfall data or the identifier of the rainfall
        data source.
        Currently, only "gpex" is supported.
        If None, a default assumption is used.
    region: GeoDataFrame or similar, optional
        The geographical region for which the rainfall data is relevant.
        Used to filter the data source.
        To be replaced by self.region
    data_catalog: DataCatalog or similar, optional
        A catalog or repository from which the rainfall data (e.g., "gpex")
        can be retrieved.
        To be replaced by self.data_catalog
    assumption_rainfall_intensity: float, default 40
        The assumed rainfall intensity (in mm/hr) to be used if no historical data is
        specified or available.
        This value is typically based on the capacity of the sewer system in the region.

    Returns
    -------
    function
        A function that calculates the rainfall depth (in mm) based on the
        given concentration time (in hours).
        The function uses either the derived rainfall statistics from the
        "gpex" data or the default assumption.

    See Also
    --------
    _get_rainfall_intensity_assumption_from_gpex
    """
    if rainfall_fn == "gpex":
        logger.info("Get rainfall functions from gpex (return period 2 year)")
        # read gpex data
        ds = data_catalog.get_rasterdataset(
            "gpex",
            bbox=region.to_crs(4326).total_bounds,
            buffer=1,  # ensure there is data
        )
        # get func from data
        calculate_rainfall_depth_given_concentration_time = _get_rainfall_from_gpex(
            ds.isel(lat=0, lon=0)
        )

    else:
        logger.error("NotImplementedError: cannot use rainfall file {rainfall_fn}.")
        logger.debug(f"Fall back to default {assumption_rainfall_intensity} mm/hr.")

        def calculate_rainfall_depth_given_concentration_time(tc):
            return assumption_rainfall_intensity * tc

    return calculate_rainfall_depth_given_concentration_time
