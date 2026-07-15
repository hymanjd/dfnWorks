import networkx as nx 
import pandas as pd 
import numpy as np
from scipy.stats import hmean

def make_connection_data_frame(self, G):
    """
    Create a PFLOTRAN unstructured grid connections DataFrame from a NetworkX graph.

    Loops through graph intersection nodes, identifies valid fracture-fracture
    connections, computes the connection area using the harmonic mean of fracture
    apertures, and returns the result as a DataFrame. Only internal connections
    are assigned.

    Parameters
    ----------
    G : networkx.Graph
        Graph containing fracture nodes and intersection nodes. Fracture nodes
        must have an ``aperture`` attribute; intersection nodes must have
        ``x``, ``y``, ``z``, and ``length`` attributes.

    Returns
    -------
    pd.DataFrame
        DataFrame of graph connections with columns
        ``i``, ``j``, ``xc``, ``yc``, ``zc``, and ``area``.
    """
    conns = []
    for u in G.intersections:
        frac_list = list(G.neighbors(u))
        if len(frac_list) == 2:
            frac1 = frac_list[0]
            frac2 = frac_list[1]
            if frac1 != 's' and frac2 != 's' and frac1 != 't' and frac2 != 't': 
                b1 = G.nodes[frac1]['aperture']
                try:
                    b2 = G.nodes[frac2]['aperture']
                    b_hmean = hmean([b1, b2])
                except KeyError:
                    b2 = None 
                    b_hmean = b1
                x = G.nodes[u]['x']
                y = G.nodes[u]['y']
                z = G.nodes[u]['z']
                G.nodes[u]['area'] = G.nodes[u]['length'] * b_hmean
                conns.append((frac1, frac2, x, y, z, G.nodes[u]['area']))
    conns_df = pd.DataFrame(conns, columns=["i", "j", "xc", "yc", "zc", "area"])
    return conns_df     


def add_boundary_nodes_as_conns(self, G, altered_cells_df, conns_df, boundaries_df):
    """
    Add boundary-node connections to an existing connection DataFrame.

    Boundary-to-internal connections from dfnWorks ``.uge`` files connect at
    the cell center, omitting the half-cell segment between the boundary point
    and that center. This function amends the geometry by nudging each boundary
    point 1 % of the way toward its parent cell center, computes the
    corresponding connection areas, matches each boundary row to its added
    boundary cell, and appends the new connections to the existing DataFrame.

    Parameters
    ----------
    G : networkx.Graph
        Graph containing fracture nodes with ``aperture`` attributes.
    altered_cells_df : pd.DataFrame
        Cell DataFrame containing both internal cells and added boundary cells.
        Boundary rows are identified by ``row_type == "boundary"``.
    conns_df : pd.DataFrame
        Existing DataFrame of internal fracture-fracture connections.
    boundaries_df : pd.DataFrame
        Boundary DataFrame containing boundary locations (``x``, ``y``, ``z``),
        ``length``, fracture ``id`` values, and negative boundary IDs
        (``neg_id``).

    Returns
    -------
    pd.DataFrame
        Updated connection DataFrame containing both the original internal
        connections and the new boundary connections, with columns
        ``i``, ``j``, ``xc``, ``yc``, ``zc``, and ``area``.
    """
    boundaries_df_copy = boundaries_df.copy()
    altered_cells_df_copy = altered_cells_df.copy()
    conns_df_copy = conns_df.copy()

    # preserve original geometry for later association
    boundaries_df_copy['x_orig'] = boundaries_df_copy['x']
    boundaries_df_copy['y_orig'] = boundaries_df_copy['y']
    boundaries_df_copy['z_orig'] = boundaries_df_copy['z']
    boundaries_df_copy['length_orig'] = boundaries_df_copy['length']

    for i in range(1, self.num_frac + 1):
        if i in boundaries_df_copy['id'].values:  
            X1 = np.array(altered_cells_df_copy.loc[altered_cells_df_copy['id'] == i, ['x', 'y', 'z']].values[0])

            mask = boundaries_df_copy['id'] == i

            midpoints = boundaries_df_copy.loc[mask, ['x', 'y', 'z']].apply(
                lambda row: np.array(row) + 0.01 * (X1 - np.array(row)),
                axis=1,
                result_type='expand'
            )
            midpoints.columns = ['x', 'y', 'z']
            boundaries_df_copy.loc[mask, 'x'] = midpoints['x'].values
            boundaries_df_copy.loc[mask, 'y'] = midpoints['y'].values
            boundaries_df_copy.loc[mask, 'z'] = midpoints['z'].values

            boundaries_df_copy.loc[mask, "length"] *= G.nodes[i]['aperture']

    boundary_cells = altered_cells_df_copy.loc[
        altered_cells_df_copy['row_type'] == 'boundary',
        ['id', 'parent_cell_id', 'neg_id', 'x', 'y', 'z', 'orig_length']
    ].rename(columns={'id': 'altered_id'})

    merged = boundaries_df_copy.merge(
        boundary_cells,
        left_on=['neg_id', 'x_orig', 'y_orig', 'z_orig', 'length_orig'],
        right_on=['neg_id', 'x', 'y', 'z', 'orig_length'],
        how='left',
        validate='one_to_one'
    )

    new_conns = merged.rename(columns={
        'altered_id': 'i',
        'id': 'j',
        'x_x': 'xc',
        'y_x': 'yc',
        'z_x': 'zc',
        'length': 'area'
    })[['i', 'j', 'xc', 'yc', 'zc', 'area']]

    conns_df_updated = pd.concat([conns_df_copy, new_conns], ignore_index=True)

    return conns_df_updated




def make_cell_data_frame(self, G):
    """
    Create a PFLOTRAN unstructured grid cell DataFrame from a NetworkX graph.

    Assigns volumes and center coordinates to each fracture node, then
    collects those properties into a DataFrame. Only internal cells are
    assigned.

    Parameters
    ----------
    G : networkx.Graph
        Graph containing fracture nodes. Each node at index ``i``
        (1-based) must be addressable via ``G.nodes[i]``.

    Returns
    -------
    pd.DataFrame
        DataFrame of cell properties with columns
        ``id``, ``x``, ``y``, ``z``, and ``volume``.
    """
    self.volume = self.surface_area * self.aperture

    for i in range(1, self.num_frac + 1):
        G.nodes[i]['volume'] = self.volume[i-1]
        G.nodes[i]['center'] = self.centers[i-1]

    cells = []
    for i in range(1, self.num_frac + 1):
        cid = i
        x = G.nodes[i]['center'][0]
        y = G.nodes[i]['center'][1]
        z = G.nodes[i]['center'][2]
        volume = G.nodes[i]['volume']
        cells.append((cid, x, y, z, volume))

    cells_df = pd.DataFrame(cells, columns=["id", "x", "y", "z", "volume"])
    return cells_df

def add_boundary_nodes_as_cells(self, G, cells_df, boundaries_df, omega=1e-6):     
    """
    Add boundary cells to an existing cell DataFrame, conserving total volume.

    Creates one new cell per row in *boundaries_df*, assigns each new cell a
    volume equal to *omega* times the parent cell's original volume, and
    reduces the parent cell volume by the same total amount so that volume is
    conserved exactly.

    Parameters
    ----------
    G : networkx.Graph
        Graph containing fracture nodes (not directly mutated here, but passed
        for interface consistency).
    cells_df : pd.DataFrame
        DataFrame of original internal-only cells with columns
        ``id``, ``x``, ``y``, ``z``, and ``volume``.
    boundaries_df : pd.DataFrame
        Boundary DataFrame containing boundary locations (``x``, ``y``, ``z``),
        ``length``, parent fracture IDs (``id``), and boundary IDs
        (``neg_id``).
    omega : float, optional
        Fraction of the parent cell's original volume assigned to each added
        boundary cell. The parent cell volume is reduced by
        ``omega * n_boundary_cells`` so that the total volume is conserved.
        Default is ``1e-6``.

    Returns
    -------
    pd.DataFrame
        Updated cell DataFrame containing both the original internal cells
        (with reduced volumes) and the newly added boundary cells.
        Additional columns ``row_type``, ``parent_cell_id``,
        ``boundary_index``, and ``orig_length`` are included. Row IDs are
        reassigned as a contiguous 1-based index.

    Raises
    ------
    SystemExit
        If the sum of the modified parent cell volume and its associated
        boundary cell volumes does not equal the original parent cell volume
        within an absolute tolerance of ``1e-9``.
    """
    import numpy as np
    import pandas as pd

    boundaries_df_copy = boundaries_df.copy()
    cells_df_copy = cells_df.copy()
    
    # Preserve metadata on original cells but make all boundarys nan
    cells_df_copy["row_type"] = "internal"
    cells_df_copy["parent_cell_id"] = np.nan
    cells_df_copy["boundary_index"] = np.nan
    cells_df_copy["orig_length"] = np.nan

    # Preserve metadata on boundary rows before modifying length 
    boundaries_df_copy["row_type"] = "boundary"
    boundaries_df_copy["parent_cell_id"] = boundaries_df_copy["id"]
    boundaries_df_copy["boundary_index"] = boundaries_df_copy["neg_id"]
    boundaries_df_copy["orig_length"] = boundaries_df_copy["length"]

    cell_xyz_map = cells_df_copy.set_index("id")[["x", "y", "z"]].to_dict("index")
    cell_volume_map = cells_df_copy.set_index("id")["volume"].to_dict()

    # Identify all cells associated with a boundary and then assign a percentage of volume and reduce parent cell volume
    for i in boundaries_df_copy["id"].unique():
        xyz = cell_xyz_map[i]
        X1 = np.array([xyz["x"], xyz["y"], xyz["z"]], dtype=float)

        mask = boundaries_df_copy["id"] == i
        boundary_xyz = boundaries_df_copy.loc[mask, ["x", "y", "z"]].to_numpy(dtype=float)

        # Determine distance from parent cell center to boundary
        L_values = np.linalg.norm(boundary_xyz - X1, axis=1)

        # Determine original volume of cell
        volume_check_og = float(cell_volume_map[i])

        # Determine how many cells will be added
        n_additions = len(L_values)

        #calculate percentage that will be given to new cells and removed from old
        total_new_fraction = omega * n_additions
        original_fraction = 1.0 - total_new_fraction

        # original cell keeps the remaining percentage
        original_volume = volume_check_og * original_fraction
        cells_df_copy.loc[cells_df_copy["id"] == i, "volume"] = original_volume

        # each new boundary gets omega % of original volume
        boundary_volume = volume_check_og * omega
        boundaries_df_copy.loc[mask, "length"] = boundary_volume

        # check: cell + all boundaries should equal original
        new_total = original_volume + boundaries_df_copy.loc[mask, "length"].sum()

        diff = volume_check_og - new_total
        # print(f"Volume difference for {i} = {diff:.6f}")


        # Exit with error if volume is not preserved
        if not np.isclose(new_total, volume_check_og, atol=1e-9):
            print(f"ERROR: volume not conserved for cell {i}")
            exit()

    altered_cells_df = pd.concat([
        cells_df_copy,
        boundaries_df_copy.rename(columns={'length': 'volume'})
    ]).reset_index(drop=True)

    altered_cells_df['id'] = altered_cells_df.index + 1
    return altered_cells_df

def convert_graph_to_data_frames(self, G, boundaries_df):
    """
    Convert a NetworkX graph and boundary DataFrame into PFLOTRAN-style
    unstructured grid cell and connection DataFrames.

    Orchestrates four sequential steps:

    1. Build the internal cell DataFrame via :func:`make_cell_data_frame`.
    2. Build the internal connection DataFrame via
       :func:`make_connection_data_frame`.
    3. Extend the cell DataFrame with boundary cells via
       :func:`add_boundary_nodes_as_cells`.
    4. Extend the connection DataFrame with boundary connections via
       :func:`add_boundary_nodes_as_conns`.

    Parameters
    ----------
    G : networkx.Graph
        NetworkX graph representation of the discrete fracture network.
    boundaries_df : pd.DataFrame
        Boundary DataFrame used to create both boundary cells and boundary
        connections.

    Returns
    -------
    altered_cells_df : pd.DataFrame
        Updated cell DataFrame with boundary nodes added as cells.
    conns_df_updated : pd.DataFrame
        Updated connection DataFrame with boundary connections appended.
    """

    cells_df = self.make_cell_data_frame(G)
    conns_df = self.make_connection_data_frame(G)  
    altered_cells_df = self.add_boundary_nodes_as_cells(G, cells_df, boundaries_df)
    conns_df_updated = self.add_boundary_nodes_as_conns(G, altered_cells_df, conns_df, boundaries_df)
    return altered_cells_df, conns_df_updated