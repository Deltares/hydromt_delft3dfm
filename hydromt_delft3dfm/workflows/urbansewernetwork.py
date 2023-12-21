"""Workflows to approximate urban sewer networks from OSM data."""

import logging
from pathlib import Path
from typing import Optional, Union

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# hydromt
from hydromt import DataCatalog

from hydromt_delft3dfm import graph_utils, workflows
from hydromt.stats import extremes
import hydromt

logger = logging.getLogger(__name__)

__all__ = [
    "setup_urban_sewer_network_topology_from_osm",
    "setup_urban_sewer_network_bedlevel_from_dem",
    # graph additions
    # "setup_graph_nodes",
    # "setup_graph_egdes",
    # graph workflows
    "optimise_pipe_topology",
    "calculate_hydraulic_parameters",
    # TODO to be categorized
    "setup_network_connections_based_on_flowdirections",
    "setup_network_parameters_from_rasters",
    "setup_network_topology_optimization",
    "setup_network_dimentions_from_rainfallstats",
]


def setup_urban_sewer_network_topology_from_osm(
    region,
    highway_types: list[str] = None,
    waterway_types: list[str] = None,
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
    - region : object
        The geographical region for which the sewer network is to be constructed.
        This could be a bounding box, polygon, or any other representation of a region.
    - highway_types : list
        List of highway types (e.g., ["motorway", "primary"]) to include
        from the OpenStreetMap data.
    - waterway_types : list
        List of waterway types (e.g., ["river", "stream"]) to include
        from the OpenStreetMap data.
    - logger : Logger object
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
    _required_columns = None

    # 1. Build the graph for waterways (open system)
    logger.info(f"Download waterway {waterway_types} form osm")
    graph_osm_open = workflows.create_graph_from_openstreetmap(
        region=region,
        buffer=1000,
        osm_key="waterway",
        osm_values=waterway_types,
        logger=logger,
    )

    # 2. Build the graph for highways (closed system)
    logger.info(f"Download highway {highway_types} form osm")
    graph_osm_closed = workflows.create_graph_from_openstreetmap(
        region=region,
        buffer=500,
        osm_key="highway",
        osm_values=highway_types,
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

    # TODO select required columns

    # 4. Recreate the graph from the complete system topology
    logger.info("Create graph and update open system, closed system and outlets.")
    graph = workflows.create_graph_from_geodataframe(
        pd.concat([open_system, closed_system]),
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

    return graph


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
        graph_pipes, target_query='nodetype=="outlet"', weight="length", report="1"
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


def _optimize_pipe_topology_based_on_network_simplex(graph, logger=logger):
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
    area_buffer_distance: float = 30,
    flow_velocity: float = 1,
    rainfall_depth_function=lambda x: 20,
):
    """
    Calculate radius for a pipe network.

    Consider various geographic and design assumptions.

    Parameters
    ----------
    pipe_network_graph: networkx.DiGraph
        A directed graph representing the pipe network, where edges have 'geometry' and
        'length' attributes.
    area_buffer_distance: float
        Estimated distance from the pipe network connected to the drainage system (m).
        Beyond this distance, water is assumed not to drain into the system.
        This is an estimate and may vary by region.
        By default, 30 m is used.
    flow_velocity : float
        Average flow velocity in the pipes during near-full conditions (m/s).
        This is an empirical value and may vary by region.
        By default 1 m/s is used.
    rainfall_depth_function : function
        A function that takes concentration time as input and returns rainfall depth
        (mm).
        This can be a constant value for simpler models,
        or a more complex function reflecting design criteria which may vary by region.
        By default, a constant value of 20 mm is returned.

    Returns
    -------
    networkx.DiGraph: The graph with updated 'radius' attributes for each edge.
    """
    # Compute contributing area for each edge in the pipe network
    nx.set_edge_attributes(
        pipe_network_graph,
        {
            (edge_start, edge_end): edge_data["geometry"]
            .buffer(area_buffer_distance)
            .area
            for edge_start, edge_end, edge_data in pipe_network_graph.edges(data=True)
        },
        name="contributing_area",
    )

    # Initialize pipe radius to zero for all edges
    nx.set_edge_attributes(
        pipe_network_graph,
        {
            (edge_start, edge_end): 0
            for edge_start, edge_end in pipe_network_graph.edges()
        },
        name="radius",
    )

    if not nx.is_directed_acyclic_graph(pipe_network_graph):
        raise ValueError("The input graph must be a directed acyclic graph (DAG)")

    # Calculate pipe radius in topological order
    for current_node in nx.topological_sort(pipe_network_graph):
        for start_node, target_node in pipe_network_graph.out_edges(current_node):
            # Subgraph of all upstream edges to the target node
            upstream_graph = pipe_network_graph.subgraph(
                nx.dfs_tree(pipe_network_graph.reverse(), target_node).nodes
            )

            # Summation of properties for upstream pipes
            total_upstream_area = sum(
                edge_data
                for _, _, edge_data in upstream_graph.edges(data="contributing_area")
            )
            total_upstream_length = sum(
                edge_data for _, _, edge_data in upstream_graph.edges(data="length")
            )

            # Upstream hydraulic calculations
            concentration_time_upstream = total_upstream_length / flow_velocity
            rainfall_depth = rainfall_depth_function(concentration_time_upstream)
            volume_upstream = rainfall_depth * total_upstream_area

            # Average radius to accommodate upstream volume
            avg_radius_needed = (
                np.sqrt(volume_upstream / (np.pi * total_upstream_length))
                if total_upstream_length > 0
                else 0
            )

            # Ensure one pipe is left to calculate the radius
            total_upstream_radius = sum(
                edge_data for _, _, edge_data in upstream_graph.edges(data="radius")
            )
            num_pipes_with_radius = sum(
                1 for _, _, radius in upstream_graph.edges(data="radius") if radius != 0
            )
            assert (
                num_pipes_with_radius == len(upstream_graph.edges()) - 1
            ), "Exactly one pipe should be without a radius"

            # Calculate current pipe radius
            current_pipe_radius = (
                avg_radius_needed * len(upstream_graph.edges()) - total_upstream_radius
            )
            pipe_network_graph[start_node][target_node]["radius"] = current_pipe_radius

    return pipe_network_graph


def setup_network_connections_based_on_flowdirections(
    user_input_graph: nx.graph,
    dem_fn: Union[str, Path],
    logger: logging.Logger = logger,
    **kwargs,
) -> nx.graph:
    """
    Set up flow directions from a Digital Elevation Model (DEM).

    The method performs the following steps:
    1. Retrieves a flow direction raster from a DEM.
    2. Converts the flow direction raster into a graph and a network, using the network
    to represent flow directions.
    3. Snaps nodes from the `user_input_graph` to the flow direction network to ensure
    that the network aligns with natural flow paths.
    4. Identifies and fills in any missing links in the `user_input_graph`, ensuring
    a fully connected network.
    5. Fills in any missing geometry in the `user_input_graph` and converts it into
    a network.

    Parameters
    ----------
    user_input_graph: nx.graph
        The initial graph provided by the user, which will be updated to align
        with the flow direction network.

    dem_fn: Union[str, Path]
        The file path to the digital elevation model (DEM) used to derive the
        flow direction raster. This can be provided as a string or a Path object.

    **kwargs:
        Other keyword arguments that may be required for specific implementations
        of the function.

    Returns
    -------
    nx.graph:
        The processed `user_input_graph`, now an instance of nx.graph class,
        with any missing links and geometry filled in based on flow directions.
        The graph represents the urban drainage network.

    Notes
    -----
    The function reprojects the DEM-derived flow network to the local CRS
    (Coordinate Reference System) of the `user_input_graph` to ensure
    that all spatial data aligns correctly.
    """
    # method implementation goes here
    pass


def setup_network_parameters_from_rasters(
    graph: nx.graph,
    dem_fn,
    landuse_fn,
    water_demand_fn,
    population_fn,
    building_footprint_fn,
    logger: logging.Logger = logger,
    **kwargs,
) -> nx.graph:
    """
    Set up physical parameters of pipes using raster data sources.

    The method performs the following steps:
    1. Adds upstream and downstream street levels to the graph using a
    digital elevation model (DEM).
    2. Adds upstream areas to the graph using land use data.
    3. Adds water demand, population, and building footprint data to the graph using
    corresponding raster data sources.
    4. Approximates the size of the network based on the street size/type using data
    from the graph itself. --> use as weight for the next step.
    5. Approximates the slope of the network based on upstream and downstream street
    levels, the depths, and the length of the edge geometry in the graph.
    --> use as weight for the next step.

    Parameters
    ----------
    graph: nx.graph
        The network graph to be filled with physical parameters.
        This is an instance of the nx.graph class, which wraps around NetworkX Graph.

    dem_fn: Union[str, Path]
        The file path to the digital elevation model (DEM) used to derive the upstream
        and downstream street levels.

    landuse_fn: Union[str, Path]
        The file path to the land use data used to derive the upstream area.

    water_demand_fn: Union[str, Path]
        The file path to the water demand data used to add water demand information
        to the graph.

    population_fn: Union[str, Path]
        The file path to the population data used to add population information
        to the graph.

    building_footprint_fn: Union[str, Path]
        The file path to the building footprint data used to add building footprint
        information to the graph.

    **kwargs:
         Additional keyword arguments for more specific implementations of the function.

    Returns
    -------
    nx.graph:
        The updated graph, an instance of nx.graph class, representing the
        urban drainage network with physical parameters filled in.
    """
    # method implementation goes here
    pass


def setup_network_topology_optimization(
    graph: nx.graph, method: str, logger: logging.Logger = logger, **kwargs
) -> nx.graph:
    """
    Optimise the topology of the urban drainage network represented by the graph.

    The method performs the following steps:
    1. Computes the selected graph metric (defined by 'method') for all nodes
    in the graph.
    2. Removes unnecessary links based on the computed graph metric.
    3. Iterates the above steps as needed until the network is optimized.

    Parameters
    ----------
    graph: nx.graph
        The network graph to be optimized. This is an instance of the nx.graph class,
        which wraps around the NetworkX Graph class.

    method: str
        The method to use for optimization. This determines the graph metric that will
        be computed and used for the optimization process.

    **kwargs:
        Additional keyword arguments for more specific implementations of the function.

    Returns
    -------
    nx.graph:
        The optimized graph, an instance of nx.graph class, representing the
        urban drainage network.

    Notes
    -----
    The betweenness centrality of a node in a graph is a measure of how often it
    appears on the shortest paths between nodes. Nodes with high betweenness centrality
    can be considered 'hubs' in the network. By removing links with low betweenness
    centrality, we can potentially simplify the network topology without significantly
    impacting the overall connectivity.

    The optimal way to optimize network topology will highly depend on the specific
    characteristics of the network, which are often determined by previous steps in the
    network analysis process.

    Other candidate graph metrics that might be useful for network topology optimization
    include:
    - Degree centrality: Measures the number of edges connected to a node. This can
    be useful for identifying nodes that are most connected.
    - Closeness centrality: Measures how close a node is to all other nodes
    in the network. This can be useful for identifying nodes that are centrally located.
    - Eigenvector centrality: Measures a node's influence based on the number of links
    it has to other influential nodes. This can be useful for identifying influential
    nodes in the network.

    If the network's weights have physical meanings, shortest path or
    maximum flow methods can also be considered for optimization.
    These methods, such as the Ford-Fulkerson or Edmonds-Karp algorithm,
    can be used to identify critical paths in the network. This can be particularly
    useful for identifying vulnerabilities in the network or for identifying potential
    areas for network enhancement.

    See Also
    --------
    https://github.com/xldeltares/hydrolib-nowcasting/tree/master

    """
    # method implementation goes here
    pass


def setup_network_dimentions_from_rainfallstats(
    graph: nx.graph,
    rainfall_fn: Union[str, Path],
    rainfall_assumption,
    capacity_assumption=None,
    diameter_range=None,
    velocity_range=None,
    slope_range=None,
    region=None,
    logger: logging.Logger = logger,
) -> nx.graph:
    """
    Update the dimensions of the pipes using historical rainfall data.

    After the graph topology is optimized, this method updates the physical parameters
    of the network, such as the diameter of the pipes. Rainfall statistics are obtained
    from historical rainfall datasets, and assumptions are made about the capacity of
    the sewer network (for example, based on the rainfall return period).
    Alternatively, the problem can be simplified by directly assuming a
    capacity (in mm/hr) for the network, bypassing the need for historical rainfall
    data.

    The dimensions of the pipe (such as diameter) are then approximated based on the
    edge weights used for network optimization, subject to reasonable ranges for
    diameter, velocity, and slope.

    Parameters
    ----------
    graph: nx.graph
        The network graph to be updated with new pipe dimensions. This is an instance
        of the nx.graph class, which wraps around the NetworkX Graph class.

    rainfall_fn: Union[str, Path]
        The file path to the historical rainfall data used to derive rainfall stats.

    rainfall_assumption: float
        An assumption about the capacity of the sewer system based on the rainfall
        return period (in years).
        This is used to derive statistics from the historical rainfall data.

    region: GeoDataFrame
        Geometry file of the bounding area.
        
    capacity_assumption: float, optional
        An assumption about the capacity of the sewer network (in mm/hr).
        If provided, this value is used instead of deriving capacity from historical
        rainfall data.

    diameter_range: tuple, optional
        A tuple specifying the minimum and maximum feasible pipe diameters.
        Used to constrain the approximation of pipe dimensions.

    velocity_range: tuple, optional
        A tuple specifying the minimum and maximum feasible pipe velocities.
        Used to constrain the approximation of pipe dimensions.

    slope_range: tuple, optional
        A tuple specifying the minimum and maximum feasible pipe slopes.
        Used to constrain the approximation of pipe dimensions.

    Returns
    -------
    nx.graph:
        The updated graph, an instance of nx.graph class,
        representing the urban drainage network with new pipe dimensions.

    References
    ----------
    - hydromt.stats.extremes:
    Function used to derive statistics from historical rainfall data.
    """
    # method implementation goes here

    # 1: Get T2 rainfall amount per region
    # draft substeps:
    # - read data from catalog based on region
    # -     * optionally visualize data to get an idea of quality
    # - average precip for region
    # - do T2 calculation for a range of return times & event durations
    
    dc = DataCatalog(logger=logger, data_libs=['deltares_data'])

    # get era5 precipitation data cut down to region area as a data array
    # next line causes ValueError: Invalid pattern: '**' can only be an entire path component
    # possibly related to era5_hourly path; * after {year} (yml line 522)
    precip_da = dc.get_rasterdataset("era5_hourly", geom=region)["precip"]

    # average grid precipitation to reach a timeseries of avg. precipitation for region
    mean_precip = precip_da.mean(dim=["x", "y"])

    # peaks
    # TODO: Implement peak selection that includes event duration; e.g. using moving
    bm_peaks = extremes.get_peaks(precip_da, ev_type="POT", period=f"{rainfall_assumption}Y") 


    # plot for quick visualization:
    fig, ax = plt.subplots(figsize=(10, 3))
    precip_da.to_pands().plot(
        ax=ax, xlabel="time", ylabel="discharge [m3/s]", color=["orange", "green"]
    )
    bm_peaks.to_pandas().plot(
    ax=ax,
    marker="o",
    linestyle="none",
    legend=False,
    color=["darkorange", "darkgreen"],
    markersize=4,
    )

    # fit EV distribution
    precip_da_params = extremes.fit_extremes(bm_peaks, ev_type="BM") #TODO: decide on desired criterium & distribution
    precip_da_params.load()

    
    # TODO: Retrieve return times with extremes.get_return_times(), see https://deltares.github.io/hydromt/latest/_examples/doing_extreme_value_analysis.html


    # 2: Compute storage needed per pipe    
    #draft substeps:
    # - multiply the T2 precipitation by node area for each node
    # - assuming each node has one output node: add upstream nodes for each downstream node

    # 3: Convert storage volumes to pipe diameter
    # - divide volumes by a timestep to arrive at a volume per time that needs to be discharged
    # - divide by pipelength to arrive at pipe diameters

    pass


# func from hybridurb
def update_edges_attributes(
    graph: nx.Graph,
    edges: gpd.GeoDataFrame,
    id_col: str = "id",
) -> nx.Graph():
    """Update the graph by adding new edges attributes specified in edges.

    Only edges with matching ids specified in "id" will be updated.
    """
    # graph df
    _graph_df = nx.to_pandas_edgelist(graph).set_index("id")
    _graph_df["_graph_edge_tuple"] = list(graph.edges)

    # check if edges id in attribute df
    if edges.index.name == id_col:
        edges.index.name = "id"
    elif id_col in edges.columns:
        edges = edges.set_index(id_col)
        edges.index.name = "id"
    else:
        raise ValueError(
            "attributes could not be updated to graph: could not perform join"
        )

    # last item that isnt NA
    graph_df = _graph_df.reindex(
        columns=_graph_df.columns.union(edges.columns, sort=False)
    )
    graph_df.update(edges)

    # add each attribute
    for c in edges.columns:
        dict = {row._graph_edge_tuple: row[c] for i, row in graph_df.iterrows()}
        nx.set_node_attributes(graph, dict, c)

    return graph


# func from hybridurb
def update_nodes_attributes(
    graph: nx.Graph,
    nodes: gpd.GeoDataFrame,
    id_col: str = "id",
) -> nx.Graph():
    """Update the graph by adding new edges attributes specified in edges.

    Only edges with matching ids specified in "id" will be updated.
    """
    # graph df
    _graph_df = pd.DataFrame(graph.nodes(data="id"), columns=["tuple", "id"]).set_index(
        "id"
    )

    # check if edges id in attribute df
    if nodes.index.name == id_col:
        nodes.index.name = "id"
    elif id_col in nodes.columns:
        nodes = nodes.set_index(id_col)
        nodes.index.name = "id"
    else:
        raise ValueError(
            "attributes could not be updated to graph: could not perform join"
        )

    # last item that isnt NA
    _graph_df = _graph_df.reindex(
        columns=_graph_df.columns.union(nodes.columns, sort=False)
    )
    graph_df = pd.concat([_graph_df, nodes]).groupby(level=0).last()
    graph_df = graph_df.loc[_graph_df.index]

    # add each attribute
    for c in nodes.columns:
        dict = {row.tuple: row[c] for i, row in graph_df.iterrows()}
        nx.set_node_attributes(graph, dict, c)

    return graph


def reverse_edges_on_negative_weight(
    graph: Union[nx.MultiDiGraph, nx.DiGraph], weight: str = "gradient"
) -> Union[nx.MultiDiGraph, nx.DiGraph]:
    """Reverse graph edges based on a negative weight attribute.

    Parameters
    ----------
    graph : Union[nx.DiGraph,nx.MultiDiGraph]
        Input graph.
    weight : str
        Name of the edge attribute to check. Default is 'gradient'.

    Returns
    -------
    Union[nx.DiGraph,nx.MultiDiGraph]
        The function modifies the graph in-place and returns it.
    """
    if isinstance(graph, nx.MultiDiGraph):
        # Collect edges that need to be reversed for MultiDiGraph
        edges_to_reverse = [
            (u, v, key, data)
            for u, v, key, data in graph.edges(keys=True, data=True)
            if data.get(weight, 0) < 0
        ]
        for u, v, key, data in edges_to_reverse:
            # Remove original edge
            graph.remove_edge(u, v, key=key)
            # Add reversed edge with preserved attributes
            data[weight] = data[weight] * -1.0
            graph.add_edge(v, u, key=key, **data)

    elif isinstance(graph, nx.DiGraph):
        # Collect edges that need to be reversed for DiGraph
        edges_to_reverse = [
            (u, v, data)
            for u, v, data in graph.edges(data=True)
            if data.get(weight, 0) < 0
        ]
        for u, v, data in edges_to_reverse:
            # Remove original edge
            graph.remove_edge(u, v)
            # Add reversed edge with preserved attributes and reversed weights
            data[weight] = data[weight] * -1.0
            graph.add_edge(v, u, **data)

    else:
        raise ValueError("The graph should be either a DiGraph or a MultiDiGraph.")

    return graph
