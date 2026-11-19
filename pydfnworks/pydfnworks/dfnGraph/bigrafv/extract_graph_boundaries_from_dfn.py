import pandas as pd

def parse_boundary_ex(self, G, intersection_list_path):
    """
    Read ``intersection_list.dat`` and return its boundary rows.

    A boundary row joins one domain face (negative id ``-1`` to ``-6``) to
    one fracture (positive id).

    Parameters
    ----------
    self : object
        DFN Class.
    G : networkx.Graph
        Bipartite graph of the DFN. Not used; kept so all BiGraFV
        functions take the same arguments.
    intersection_list_path : str
        Path to the dfnWorks ``intersection_list.dat`` file. Each line is
        ``f1 f2 x y z length``, whitespace separated.

    Returns
    -------
    pandas.DataFrame
        Boundary table with columns ``neg_id``, ``id``, ``x``, ``y``, ``z``
        and ``length``, where ``id`` is the fracture id.

    Notes
    -----
    Boundary face ids:

    ==========  ===========
    ``neg_id``  Face
    ==========  ===========
    ``-1``      ``top``
    ``-2``      ``bottom``
    ``-3``      ``left_w``
    ``-4``      ``front_n``
    ``-5``      ``right_e``
    ``-6``      ``back_s``
    ==========  ===========

    Empty lines, the ``f1`` header line, lines with fewer than six fields
    and lines that are not numeric are skipped.
    """

    BOUNDARY_NAME_TO_NEG_ID = {
        "top":     -1,
        "bottom":  -2,
        "left_w":  -3,
        "front_n": -4,
        "right_e": -5,
        "back_s":  -6,
    }

    # neg_id -> list of (node_i, x, y, z, length)
    boundary_rows = {neg_id: [] for neg_id in BOUNDARY_NAME_TO_NEG_ID.values()}

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

    # one flat list in face order (-1 to -6); avoids concatenating empty frames
    all_rows = [(neg_id, *row)
                for neg_id, rows in boundary_rows.items()
                for row in rows]
    boundaries_df = pd.DataFrame(all_rows,
                                 columns=["neg_id", "id", "x", "y", "z", "length"])

    return boundaries_df

def extract_graph_boundaries_from_dfn(
    self,
    G,
    intersection_list_path='./dfnGen_output/intersection_list.dat',
    selected_boundary_ids=None,
):
    """
    Read the boundary connections from a dfnWorks intersection list.

    Calls :func:`parse_boundary_ex`, then keeps only the faces in
    ``selected_boundary_ids`` if it is given.

    Parameters
    ----------
    self : object
        DFN Class.
    G : networkx.Graph
        Bipartite graph of the DFN.
    intersection_list_path : str, optional
        Path to the dfnWorks ``intersection_list.dat`` file.
        Default is ``"./dfnGen_output/intersection_list.dat"``.
    selected_boundary_ids : list of int or None, optional
        Boundary face ids to keep, e.g. ``[-3, -5]`` for the left and
        right faces only. If ``None``, all six faces are kept.
        Default is ``None``.

    Returns
    -------
    pandas.DataFrame
        Boundary table with columns ``neg_id``, ``id``, ``x``, ``y``, ``z``
        and ``length``.

    See Also
    --------
    parse_boundary_ex : Low-level parser for ``intersection_list.dat``.
    """
    boundaries_df = self.parse_boundary_ex(
        G,
        intersection_list_path=intersection_list_path
    )

    if selected_boundary_ids is not None:
        boundaries_df = boundaries_df[
            boundaries_df["neg_id"].isin(selected_boundary_ids)
        ].reset_index(drop=True)

    return boundaries_df
