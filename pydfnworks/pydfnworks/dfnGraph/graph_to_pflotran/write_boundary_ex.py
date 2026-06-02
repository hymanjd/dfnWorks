import numpy as np
import pandas as pd

# Write one .ex file per boundary
def write_ex(filename: str, rows: list):
    """
    Write a ``.ex`` boundary file.

    Parameters
    ----------
    filename : str
        Path (including filename) of the ``.ex`` file to create or overwrite.
    rows : list of tuple
        Each tuple must contain exactly five elements in the order
        ``(node_i, x, y, z, length)`` where:

        - **node_i** (*int*) – Boundary node identifier written as an integer.
        - **x** (*float*) – X-coordinate of the node.
        - **y** (*float*) – Y-coordinate of the node.
        - **z** (*float*) – Z-coordinate of the node.
        - **length** (*float*) – Boundary face connection length value.

    Returns
    -------
    None
        Writes directly to *filename*; no return value.

    Notes
    -----
    File format produced::

        CONNECTIONS    <count>
        <node_id>    <x>    <y>    <z>    <length>

    Coordinates and length are written in 16-digit scientific notation
    (``:.16e``).
    """
    with open(filename, "w") as g:
        g.write(f"CONNECTIONS\t{len(rows)}\n")
        for node_i, x, y, z, length in rows:
            g.write(
                f"{int(node_i)}\t{x:.16e}\t{y:.16e}\t{z:.16e}\t{length:.16e}\n"
            )

def write_boundary_ex_files(cells_df, omega = 0.01):
    """
    Generate one ``.ex`` boundary file for each of the six domain boundaries.

    Reads boundary-cell rows from *cells_df*, shifts each boundary point
    slightly outward from its parent cell centre by a fraction *omega*, then
    delegates file writing to :func:`write_ex`.

    Parameters
    ----------
    cells_df : pd.DataFrame
        Combined DataFrame produced by ``graph_uge`` containing both internal
        fracture-cell rows and boundary rows.  Required columns:

        - ``row_type`` (*str*) – ``"internal"`` for fracture cells; any other
          value flags a boundary row.
        - ``id`` (*int*) – Unique node identifier for the row.
        - ``x``, ``y``, ``z`` (*float*) – Spatial coordinates of the node.
        - ``neg_id`` (*int or NaN*) – Negative boundary identifier for
          boundary rows (``-1`` through ``-6``); ``NaN`` for internal rows.
        - ``parent_cell_id`` (*int*) – ``id`` of the internal cell that owns
          this boundary point.

    omega : float, optional
        Fractional outward shift applied to each boundary point relative to
        the vector from its parent cell centre.  A value of ``0.01`` moves
        the point 1 % farther from the parent centre along that vector.
        Default is ``0.01``.

    Returns
    -------
    None
        Six ``.ex`` files are written to the current working directory; no
        value is returned.

    Output Files
    ------------
    Each filename maps to a ``neg_id`` value:

    ========================  =========
    File                      ``neg_id``
    ========================  =========
    ``boundary_top.ex``       ``-1``
    ``boundary_bottom.ex``    ``-2``
    ``boundary_left_w.ex``    ``-3``
    ``boundary_front_n.ex``   ``-4``
    ``boundary_right_e.ex``   ``-5``
    ``boundary_back_s.ex``    ``-6``
    ========================  =========

    Notes
    -----
    The outward-shift formula applied to each boundary point *P* given its
    parent centre *X1* is:

    .. math::

        P_{\\text{shifted}} = X_1 + (1 + \\omega)(P - X_1)

    Boundary rows whose ``parent_cell_id`` is not found in the internal rows
    are silently skipped.
    """
    BOUNDARY_NAME_TO_NEG_ID = {
        "top":     -1,
        "bottom":  -2,
        "left_w":  -3,
        "front_n": -4,
        "right_e": -5,
        "back_s":  -6,
    }
    cells_df_copy = cells_df.copy()
    # print(f'Cells df copy ^^^^^^^^^^^^^ {cells_df_copy}')
    # Internal rows are original fracture cells
    internal_df = cells_df_copy[cells_df_copy["row_type"] == "internal"].copy()

    # Boundary rows are rows to write
    boundary_df = cells_df_copy[cells_df_copy["neg_id"].notna()].copy()
    boundary_df["neg_id"] = boundary_df["neg_id"].astype(int)
    boundary_df["parent_cell_id"] = boundary_df["parent_cell_id"].astype(int)

    # Map parent cell id to xyz
    parent_xyz_map = internal_df.set_index("id")[["x", "y", "z"]].to_dict("index")

    # Build map: neg_id to list of (node_i, x, y, z, length)
    boundary_rows = {neg_id: [] for neg_id in BOUNDARY_NAME_TO_NEG_ID.values()}

    for _, row in boundary_df.iterrows():
        neg_id = int(row["neg_id"])
        parent_id = int(row["parent_cell_id"])

        if neg_id not in boundary_rows:
            continue
        if parent_id not in parent_xyz_map:
            continue

        # Parent cell center X1
        X1 = np.array([
            parent_xyz_map[parent_id]["x"],
            parent_xyz_map[parent_id]["y"],
            parent_xyz_map[parent_id]["z"],
        ], dtype=float)

        # Current boundary point P
        P = np.array([row["x"], row["y"], row["z"]], dtype=float)
        P_shifted = X1 + (1.0 + omega) * (P - X1)

        boundary_rows[neg_id].append((
            int(row["id"]),
            float(P_shifted[0]),
            float(P_shifted[1]),
            float(P_shifted[2]),
            # 1 / (float(row["orig_length"])*1e-4)
            # 1 / (float(row["orig_length"])*3.464101615137754701e-06)
            1000
        ))

    # print(f'boundary_rows from .ex script ^^^^^^^^^^ {boundary_rows}')
    # Write one file per boundary
    for name, neg_id in BOUNDARY_NAME_TO_NEG_ID.items():
        filename = f"boundary_{name}.ex"
        rows = boundary_rows[neg_id]
        write_ex(filename, rows)
        # print(f"Wrote {len(rows):>6} rows -> {filename}")


def write_boundary_ex(cells_df):
    """
    Function for writing all boundary ``.ex`` files.

    This is a thin wrapper around :func:`write_boundary_ex_files` and is
    the intended external-facing call.  Invoke as::

        write_boundary_ex(cells_df)

    Parameters
    ----------
    cells_df : pd.DataFrame
        DataFrame containing both internal fracture-cell rows and boundary
        rows.  See :func:`write_boundary_ex_files` for the full column
        specification.

    Returns
    -------
    None
        Delegates entirely to :func:`write_boundary_ex_files`; no value is
        returned.

    See Also
    --------
    write_boundary_ex_files : Core implementation with full parameter control.
    write_ex : Low-level ``.ex`` file writer.
    """
    write_boundary_ex_files(cells_df)
