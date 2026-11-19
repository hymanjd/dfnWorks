import numpy as np
import pandas as pd


def make_connection_data_frame(self, G):
    r"""
    Build the internal connection table from the bipartite graph.

    Each intersection node with exactly two fracture neighbours (and no
    source or target neighbour) becomes one connection. The connection
    aperture :math:`b_A` is chosen so that a PFLOTRAN connection with
    permeability :math:`k = b^2/12` has the same cubic-law hydraulic
    resistance as the two fracture segments it replaces:

    .. math::

        b_A = \frac{l_1/b_1^2 + l_2/b_2^2}{l_1/b_1^3 + l_2/b_2^3},
        \qquad A = w \, b_A

    where :math:`b_1, b_2` are the fracture apertures, :math:`l_1, l_2` are
    the distances from each fracture center to the intersection point, and
    :math:`w` is the intersection length.

    Parameters
    ----------
    self : object
        DFN Class.
    G : networkx.Graph
        Bipartite graph of the DFN. Fracture nodes need ``aperture`` and
        ``center`` attributes (``center`` is set by
        :func:`make_cell_data_frame`). Intersection nodes need ``x``,
        ``y``, ``z`` and ``length`` attributes.

    Returns
    -------
    pandas.DataFrame
        Connection table with columns ``i``, ``j``, ``xc``, ``yc``, ``zc``
        and ``area``, where ``(xc, yc, zc)`` is the intersection point.

    Notes
    -----
    The connection area is also stored on the graph as
    ``G.nodes[u]['area']`` for each intersection node ``u``.

    Logs an error and exits if :func:`make_cell_data_frame` has not been
    run first (fracture nodes have no ``center``), or if :math:`b_A` would
    divide by zero (both distances are zero or an aperture is zero).
    """
    conns = []
    for u in G.intersections:
        frac_list = list(G.neighbors(u))
        if len(frac_list) != 2:
            continue

        frac1, frac2 = frac_list
        if frac1 in ('s', 't') or frac2 in ('s', 't'):
            continue

        for frac in (frac1, frac2):
            if 'center' not in G.nodes[frac]:
                self.print_log(
                    f"Error. Fracture {frac} has no 'center' attribute. "
                    "Run make_cell_data_frame before make_connection_data_frame.",
                    'error')

        b1 = float(G.nodes[frac1]['aperture'])
        b2 = float(G.nodes[frac2]['aperture'])

        X1 = np.asarray(G.nodes[frac1]['center'], dtype=float)
        X2 = np.asarray(G.nodes[frac2]['center'], dtype=float)

        x = float(G.nodes[u]['x'])
        y = float(G.nodes[u]['y'])
        z = float(G.nodes[u]['z'])
        gn = np.array([x, y, z], dtype=float)

        l1 = np.linalg.norm(gn - X1)
        l2 = np.linalg.norm(X2 - gn)

        if b1 <= 0.0 or b2 <= 0.0:
            self.print_log(
                f"Error. Zero or negative aperture on fracture {frac1} ({b1}) or "
                f"{frac2} ({b2}) at intersection {u}.", 'error')

        numerator = l1 / b1**2 + l2 / b2**2
        denominator = l1 / b1**3 + l2 / b2**3
        if denominator == 0.0:
            self.print_log(
                f"Error. Division by zero computing the aperture of intersection {u} "
                f"between fractures {frac1} and {frac2}: both fracture centers are at "
                "the intersection point (l1 = l2 = 0).", 'error')
        b_area = numerator / denominator

        width = float(G.nodes[u]['length'])
        area = width * b_area

        G.nodes[u]['area'] = area
        conns.append((frac1, frac2, x, y, z, area))

    conns_df = pd.DataFrame(conns, columns=["i", "j", "xc", "yc", "zc", "area"])
    return conns_df


def add_boundary_nodes_as_conns(self, G, altered_cells_df, conns_df, boundaries_df):
    """
    Append boundary connections to the internal connection table.

    In dfnWorks ``.uge`` files a boundary connection meets the fracture at
    its cell center, so the half-cell between the boundary point and the
    center is missing. To fix this, each boundary point is moved 1 % of the
    way toward its parent cell center. The connection area is the boundary
    intersection length times the fracture aperture. Each boundary row is
    then matched to its helper cell from :func:`add_boundary_nodes_as_cells`.

    Parameters
    ----------
    self : object
        DFN Class.
    G : networkx.Graph
        Bipartite graph of the DFN. Fracture nodes need an ``aperture``
        attribute.
    altered_cells_df : pandas.DataFrame
        Cell table with fracture cells and boundary helper cells, from
        :func:`add_boundary_nodes_as_cells`. Helper cells have
        ``row_type == "boundary"``.
    conns_df : pandas.DataFrame
        Internal connection table from :func:`make_connection_data_frame`.
    boundaries_df : pandas.DataFrame
        Boundary table from :func:`extract_graph_boundaries_from_dfn`, with
        columns ``neg_id``, ``id``, ``x``, ``y``, ``z`` and ``length``.

    Returns
    -------
    pandas.DataFrame
        Connection table with the internal connections followed by the
        boundary connections. Columns are ``i``, ``j``, ``xc``, ``yc``,
        ``zc`` and ``area``.
    """
    boundaries_df_copy = boundaries_df.copy()
    altered_cells_df_copy = altered_cells_df.copy()
    conns_df_copy = conns_df.copy()

    # keep the original geometry so rows can be matched to helper cells later
    boundaries_df_copy['x_orig'] = boundaries_df_copy['x']
    boundaries_df_copy['y_orig'] = boundaries_df_copy['y']
    boundaries_df_copy['z_orig'] = boundaries_df_copy['z']
    boundaries_df_copy['length_orig'] = boundaries_df_copy['length']

    for i in range(1, self.num_frac + 1):
        if i in boundaries_df_copy['id'].values:
            X1 = np.array(altered_cells_df_copy.loc[altered_cells_df_copy['id'] == i,
                                                    ['x', 'y', 'z']].values[0])

            mask = boundaries_df_copy['id'] == i

            midpoints = boundaries_df_copy.loc[mask, ['x', 'y', 'z']].apply(
                lambda row: np.array(row) + 0.01 * (X1 - np.array(row)),
                axis=1,
                result_type='expand')
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
        validate='one_to_one')

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
    Build the fracture cell table, one cell per fracture.

    The cell volume is fracture surface area times aperture, and the cell
    center is the fracture center. Both are also stored on the graph as
    ``G.nodes[i]['volume']`` and ``G.nodes[i]['center']``.

    Parameters
    ----------
    self : object
        DFN Class. ``surface_area``, ``aperture``, ``centers`` and
        ``num_frac`` must be set.
    G : networkx.Graph
        Bipartite graph of the DFN. Fracture ``i`` (1-based) must be the
        node ``G.nodes[i]``.

    Returns
    -------
    pandas.DataFrame
        Cell table with columns ``id``, ``x``, ``y``, ``z`` and ``volume``.

    Notes
    -----
    This function sets ``self.volume`` on the DFN object.
    """
    self.volume = self.surface_area * self.aperture

    for i in range(1, self.num_frac + 1):
        G.nodes[i]['volume'] = self.volume[i - 1]
        G.nodes[i]['center'] = self.centers[i - 1]

    cells = []
    for i in range(1, self.num_frac + 1):
        x, y, z = G.nodes[i]['center'][:3]
        cells.append((i, x, y, z, G.nodes[i]['volume']))

    cells_df = pd.DataFrame(cells, columns=["id", "x", "y", "z", "volume"])
    return cells_df


def add_boundary_nodes_as_cells(self, G, cells_df, boundaries_df, omega=1e-6):
    """
    Add one helper cell per boundary connection, keeping total volume fixed.

    Each helper cell gets ``omega`` times its parent cell's original volume.
    The parent cell loses the same total, so the parent and its helper
    cells add up to the original volume.

    Parameters
    ----------
    self : object
        DFN Class.
    G : networkx.Graph
        Bipartite graph of the DFN. Not used; kept so all BiGraFV
        functions take the same arguments.
    cells_df : pandas.DataFrame
        Fracture cell table from :func:`make_cell_data_frame`.
    boundaries_df : pandas.DataFrame
        Boundary table from :func:`extract_graph_boundaries_from_dfn`, with
        columns ``neg_id``, ``id``, ``x``, ``y``, ``z`` and ``length``.
    omega : float, optional
        Fraction of the parent cell's original volume given to each helper
        cell. The parent keeps ``1 - omega * n_boundary`` of its volume.
        Default is ``1e-6``.

    Returns
    -------
    pandas.DataFrame
        Cell table with the fracture cells (reduced volumes) followed by the
        helper cells. Extra columns are ``row_type`` (``"internal"`` or
        ``"boundary"``), ``parent_cell_id``, ``boundary_index``, ``neg_id``
        and ``orig_length``. The ``id`` column is renumbered 1, 2, 3, ...

    Notes
    -----
    Logs an error and exits if a parent cell would be left with zero or
    negative volume (``omega * n_boundary >= 1``), or if a parent cell and
    its helper cells do not add up to the original volume (absolute
    tolerance ``1e-9``).
    """
    boundaries_df_copy = boundaries_df.copy()
    cells_df_copy = cells_df.copy()

    # metadata columns; NaN on fracture cells
    cells_df_copy["row_type"] = "internal"
    cells_df_copy["parent_cell_id"] = np.nan
    cells_df_copy["boundary_index"] = np.nan
    cells_df_copy["orig_length"] = np.nan

    # metadata on boundary rows, saved before 'length' is overwritten with volume
    boundaries_df_copy["row_type"] = "boundary"
    boundaries_df_copy["parent_cell_id"] = boundaries_df_copy["id"]
    boundaries_df_copy["boundary_index"] = boundaries_df_copy["neg_id"]
    boundaries_df_copy["orig_length"] = boundaries_df_copy["length"]

    cell_volume_map = cells_df_copy.set_index("id")["volume"].to_dict()

    for i in boundaries_df_copy["id"].unique():
        mask = boundaries_df_copy["id"] == i
        n_additions = int(mask.sum())

        volume_original = float(cell_volume_map[i])

        if omega * n_additions >= 1.0:
            self.print_log(
                f"Error. Cell {i} has {n_additions} boundary connections, so omega = "
                f"{omega} leaves it with {1.0 - omega * n_additions:.3g} of its volume. "
                f"omega must be less than {1.0 / n_additions:.3g}.", 'error')

        # parent keeps what is left after each helper cell takes omega
        parent_volume = volume_original * (1.0 - omega * n_additions)
        cells_df_copy.loc[cells_df_copy["id"] == i, "volume"] = parent_volume

        boundary_volume = volume_original * omega
        boundaries_df_copy.loc[mask, "length"] = boundary_volume

        new_total = parent_volume + boundaries_df_copy.loc[mask, "length"].sum()
        if not np.isclose(new_total, volume_original, atol=1e-9):
            self.print_log(
                f"Error. Volume not conserved for cell {i} when adding boundary cells "
                f"(original {volume_original:.6e}, new {new_total:.6e})", 'error')

    altered_cells_df = pd.concat([
        cells_df_copy,
        boundaries_df_copy.rename(columns={'length': 'volume'})
    ]).reset_index(drop=True)

    altered_cells_df['id'] = altered_cells_df.index + 1
    return altered_cells_df


def convert_graph_to_data_frames(self, G, boundaries_df):
    """
    Build the PFLOTRAN cell and connection tables from the bipartite graph.

    Runs these steps in order:

    1. Fracture cells with :func:`make_cell_data_frame`.
    2. Internal connections with :func:`make_connection_data_frame`.
    3. Boundary helper cells with :func:`add_boundary_nodes_as_cells`.
    4. Boundary connections with :func:`add_boundary_nodes_as_conns`.

    Parameters
    ----------
    self : object
        DFN Class.
    G : networkx.Graph
        Bipartite graph of the DFN.
    boundaries_df : pandas.DataFrame
        Boundary table from :func:`extract_graph_boundaries_from_dfn`.

    Returns
    -------
    altered_cells_df : pandas.DataFrame
        Cell table including the boundary helper cells.
    conns_df_updated : pandas.DataFrame
        Connection table including the boundary connections.
    """
    cells_df = self.make_cell_data_frame(G)
    conns_df = self.make_connection_data_frame(G)
    altered_cells_df = self.add_boundary_nodes_as_cells(G, cells_df, boundaries_df)
    conns_df_updated = self.add_boundary_nodes_as_conns(G, altered_cells_df, conns_df,
                                                        boundaries_df)
    return altered_cells_df, conns_df_updated
