import numpy as np
import pandas as pd

def parse_boundary_ex(self,
    G,
    intersection_list_path
):
    """
    Parse ``intersection_list.dat`` and return a DataFrame of boundary rows.

    Reads each row of the intersection list, identifies boundary connections
    (those with one negative boundary ID and one positive node ID), and
    assembles them into a single DataFrame.

    Parameters
    ----------
    G : networkx.Graph
        NetworkX graph representation of the DFN network. Accepted for
        interface consistency but not directly used in parsing.
    intersection_list_path : str
        Path to the dfnWorks ``intersection_list.dat`` file.

        Expected file format (whitespace-delimited, one row per connection)::

            f1  f2  x  y  z  length

        Boundary rows are identified by one negative boundary ID (``-1``
        through ``-6``) and one positive node ID.

    Returns
    -------
    pd.DataFrame
        DataFrame of boundary connections with columns
        ``neg_id``, ``id``, ``x``, ``y``, ``z``, and ``length``.

        The ``neg_id`` column maps to boundary faces as follows:

        ========  =========
        ``neg_id``  Face
        ========  =========
        ``-1``    top
        ``-2``    bottom
        ``-3``    left_w
        ``-4``    front_n
        ``-5``    right_e
        ``-6``    back_s
        ========  =========

    Notes
    -----
    Lines that are empty, begin with ``f1`` (header), contain fewer than
    six fields, or cannot be parsed as numeric values are silently skipped.
    """

    BOUNDARY_NAME_TO_NEG_ID = {
        "top":     -1,
        "bottom":  -2,
        "left_w":  -3,
        "front_n": -4,
        "right_e": -5,
        "back_s":  -6,
    }

    # Build map dict: neg_id -> list of (node_i, x, y, z, length)
    boundary_rows: dict[int, list] = {neg_id: [] for neg_id in BOUNDARY_NAME_TO_NEG_ID.values()}

    with open(intersection_list_path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.lower().startswith("f1"):
                continue
            parts = s.split()
            if len(parts) < 6:
                continue
            try:
                a      = int(parts[0])
                b      = int(parts[1])
                x      = float(parts[2])
                y      = float(parts[3])
                z      = float(parts[4])
                length = float(parts[5])
            except ValueError:
                continue

            # Boundary rows if one negative boundary ID, one positive node ID
            if a < 0 and b > 0:
                neg_id, node_i = a, b
            elif b < 0 and a > 0:
                neg_id, node_i = b, a
            else:
                continue

            if neg_id in boundary_rows:
                boundary_rows[neg_id].append((node_i, x, y, z, length)) 


    boundaries_df = pd.concat([
        pd.DataFrame(rows, columns=["id", "x", "y", "z", "length"])
        .assign(neg_id=neg_id)
        for neg_id, rows in boundary_rows.items()
    ])[["neg_id", "id", "x", "y", "z", "length"]].reset_index(drop=True)

    return boundaries_df

def extract_graph_boundaries_from_dfn(
    self,
    G,
    intersection_list_path='./dfnGen_output/intersection_list.dat',
    selected_boundary_ids=None,
):
    """
    Extract boundary connections from a dfnWorks intersection list.

    Delegates parsing to :func:`parse_boundary_ex`, then optionally filters
    the result to a user-specified subset of boundary faces.

    Parameters
    ----------
    G : networkx.Graph
        NetworkX graph representation of the DFN network.
    intersection_list_path : str, optional
        Path to the dfnWorks ``intersection_list.dat`` file.
        Default is ``"./dfnGen_output/intersection_list.dat"``.
    selected_boundary_ids : list of int or None, optional
        Negative boundary IDs to retain in the output (e.g. ``[-3, -5]``
        for the left and right faces only). If ``None``, all six boundary
        faces are returned. Default is ``None``.

    Returns
    -------
    pd.DataFrame
        DataFrame of boundary connections with columns
        ``neg_id``, ``id``, ``x``, ``y``, ``z``, and ``length``,
        filtered to *selected_boundary_ids* if provided.

    See Also
    --------
    parse_boundary_ex : Low-level parser for ``intersection_list.dat``.
    """
    boundaries_df = self.parse_boundary_ex(
        G,
        intersection_list_path=intersection_list_path
    )

    # Keep only requested boundaries, if provided
    if selected_boundary_ids is not None:
        boundaries_df = boundaries_df[
            boundaries_df["neg_id"].isin(selected_boundary_ids)
        ].reset_index(drop=True)

    return boundaries_df
