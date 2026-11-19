import numpy as np

from pydfnworks.general.logging import local_print_log


def write_ex(filename, rows):
    """
    Write one PFLOTRAN ``.ex`` boundary file.

    Parameters
    ----------
    filename : str
        Path of the ``.ex`` file to write. An existing file is overwritten.
    rows : list of tuple
        One ``(node_i, x, y, z, area)`` tuple per boundary face, where
        ``node_i`` is the cell id, ``(x, y, z)`` is the face point and
        ``area`` is the face area.

    Returns
    -------
    None

    Notes
    -----
    File format::

        CONNECTIONS    <count>
        <node_i>    <x>    <y>    <z>    <area>

    Values are tab separated. Floats are written in scientific notation
    with 16 digits after the decimal point (``:.16e``).
    """
    with open(filename, "w") as g:
        g.write(f"CONNECTIONS\t{len(rows)}\n")
        for node_i, x, y, z, area in rows:
            g.write(f"{int(node_i)}\t{x:.16e}\t{y:.16e}\t{z:.16e}\t{area:.16e}\n")


def write_boundary_ex_files(cells_df, omega=0.01, area_per_dist=1e-2):
    """
    Write one ``.ex`` boundary file for each of the six domain faces.

    Each boundary helper cell in ``cells_df`` gets one face. The face point
    is the helper cell point moved outward, away from its parent cell
    center, by a fraction ``omega``.

    Parameters
    ----------
    cells_df : pandas.DataFrame
        Cell table from :func:`convert_graph_to_data_frames`, with fracture
        cells and boundary helper cells. Columns used:

        - ``row_type``: ``"internal"`` for fracture cells.
        - ``id``: cell id.
        - ``x``, ``y``, ``z``: cell point.
        - ``neg_id``: boundary face id (``-1`` to ``-6``) for helper cells,
          ``NaN`` for fracture cells.
        - ``parent_cell_id``: id of the fracture cell the helper cell
          belongs to.
    omega : float, optional
        Outward shift as a fraction of the distance from the parent cell
        center. ``0.01`` moves the point 1 % farther out. Default is ``0.01``.
    area_per_dist : float, optional
        Face area per unit distance [m]. Each face gets area
        ``area_per_dist * d``, where ``d`` is the distance from the helper
        cell point to its face point. Default is ``1e-2``.

    Returns
    -------
    None

    Notes
    -----
    Files written to the current working directory:

    ==========================  ==========
    File                        ``neg_id``
    ==========================  ==========
    ``boundary_top.ex``         ``-1``
    ``boundary_bottom.ex``      ``-2``
    ``boundary_left_w.ex``      ``-3``
    ``boundary_front_n.ex``     ``-4``
    ``boundary_right_e.ex``     ``-5``
    ``boundary_back_s.ex``      ``-6``
    ==========================  ==========

    For a helper cell point :math:`P` with parent center :math:`X_1`, the
    face point is

    .. math::

        P_{\\text{shifted}} = X_1 + (1 + \\omega)(P - X_1)

    The face area is :math:`A = R\\,d` with :math:`R` = ``area_per_dist``
    and :math:`d = |P_{\\text{shifted}} - P| = \\omega |P - X_1|`. This
    gives every boundary face the same conductance per unit permeability,
    :math:`A/d = R`. A large fixed area on faces this close to the cell
    makes the boundary connection far too strong: the outflow is quantised
    and solute piles up in the helper cells.

    Helper cells whose ``parent_cell_id`` is not a fracture cell are
    skipped. Logs an error and exits if a face point lands on its helper
    cell (``d = 0``).
    """
    BOUNDARY_NAME_TO_NEG_ID = {
        "top": -1,
        "bottom": -2,
        "left_w": -3,
        "front_n": -4,
        "right_e": -5,
        "back_s": -6,
    }
    cells_df_copy = cells_df.copy()
    internal_df = cells_df_copy[cells_df_copy["row_type"] == "internal"].copy()

    boundary_df = cells_df_copy[cells_df_copy["neg_id"].notna()].copy()
    boundary_df["neg_id"] = boundary_df["neg_id"].astype(int)
    boundary_df["parent_cell_id"] = boundary_df["parent_cell_id"].astype(int)

    parent_xyz_map = internal_df.set_index("id")[["x", "y", "z"]].to_dict("index")

    # neg_id -> list of (node_i, x, y, z, area)
    boundary_rows = {neg_id: [] for neg_id in BOUNDARY_NAME_TO_NEG_ID.values()}

    for _, row in boundary_df.iterrows():
        neg_id = int(row["neg_id"])
        parent_id = int(row["parent_cell_id"])

        if neg_id not in boundary_rows:
            continue
        if parent_id not in parent_xyz_map:
            continue

        X1 = np.array([
            parent_xyz_map[parent_id]["x"],
            parent_xyz_map[parent_id]["y"],
            parent_xyz_map[parent_id]["z"],
        ], dtype=float)

        P = np.array([row["x"], row["y"], row["z"]], dtype=float)
        P_shifted = X1 + (1.0 + omega) * (P - X1)

        dist = float(np.linalg.norm(P_shifted - P))
        if dist == 0.0:
            local_print_log(
                f"Error. Boundary face for cell {int(row['id'])} is at zero distance "
                "from the cell, so its area would be zero. Check that omega > 0 and "
                "that the boundary point is not at the fracture center.", 'error')

        boundary_rows[neg_id].append((
            int(row["id"]),
            float(P_shifted[0]),
            float(P_shifted[1]),
            float(P_shifted[2]),
            area_per_dist * dist,
        ))

    for name, neg_id in BOUNDARY_NAME_TO_NEG_ID.items():
        write_ex(f"boundary_{name}.ex", boundary_rows[neg_id])


def write_boundary_ex(cells_df, omega=0.01, area_per_dist=1e-2):
    """
    Write all six boundary ``.ex`` files.

    Calls :func:`write_boundary_ex_files`.

    Parameters
    ----------
    cells_df : pandas.DataFrame
        Cell table with fracture cells and boundary helper cells. See
        :func:`write_boundary_ex_files` for the columns used.
    omega : float, optional
        Outward shift of each face point. Default is ``0.01``.
    area_per_dist : float, optional
        Face area per unit distance [m]. Default is ``1e-2``.

    Returns
    -------
    None

    See Also
    --------
    write_boundary_ex_files : Does the work.
    write_ex : Writes one ``.ex`` file.
    """
    write_boundary_ex_files(cells_df, omega=omega, area_per_dist=area_per_dist)
