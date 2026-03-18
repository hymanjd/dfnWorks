import numpy as np
import pandas as pd

# Write one .ex file per boundary
def write_ex(filename: str, rows: list):
    with open(filename, "w") as g:
        g.write(f"CONNECTIONS\t{len(rows)}\n")
        for node_i, x, y, z, length in rows:
            g.write(
                f"{int(node_i)}\t{x:.16e}\t{y:.16e}\t{z:.16e}\t{length:.16e}\n"
            )

def write_boundary_ex_files(cells_df, omega = 0.01):
    """
    Takes output from graph_uge to write boundary .ex files based on new cells near boundaries.
    Parameters
    ----------
    cells_df : pd.DataFrame
        DataFrame containing both internal and boundary rows.
    move_fraction/omega : float
        Fractional shift relative to distance from parent cell center.
            0.01   -> 1% farther

    where boundary rows have one negative boundary ID (-1..-6) and one positive node ID.

    Output files:
        boundary_top.ex       (neg_id = -1)
        boundary_bottom.ex    (neg_id = -2)
        boundary_left_w.ex    (neg_id = -3)
        boundary_front_n.ex   (neg_id = -4)
        boundary_right_e.ex   (neg_id = -5)
        boundary_back_s.ex    (neg_id = -6)

    Each .ex file format:
        CONNECTIONS\t<count>
        <node_id>\t<x>\t<y>\t<z>\t<area>
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
    move_fraction = omega
    print(f'Cells df copy ^^^^^^^^^^^^^ {cells_df_copy}')
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

        # make distance (1 + move_fraction) times longer
        P_shifted = X1 + (1.0 + move_fraction) * (P - X1)


        boundary_rows[neg_id].append((
            int(row["id"]),
            float(P_shifted[0]),
            float(P_shifted[1]),
            float(P_shifted[2]),
            float(row["orig_length"])*1e-4
        ))
    print(f'boundary_rows from .ex script ^^^^^^^^^^ {boundary_rows}')
    # Write one file per boundary
    for name, neg_id in BOUNDARY_NAME_TO_NEG_ID.items():
        filename = f"boundary_{name}.ex"
        rows = boundary_rows[neg_id]
        write_ex(filename, rows)
        print(f"Wrote {len(rows):>6} rows -> {filename}")


def write_boundary_ex(cells_df):
    write_boundary_ex_files(cells_df)
