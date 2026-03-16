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

def write_boundary_ex_files(cells_df):
    """
    Takes output from graph_uge to write boundary .ex files based on new cells near boundaries.

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

    # Keep only boundary rows
    boundary_df = cells_df_copy[cells_df_copy["neg_id"].notna()].copy()
    boundary_df["neg_id"] = boundary_df["neg_id"].astype(int)
    # Build map: neg_id -> list of (node_i, x, y, z, length)
    boundary_rows = {neg_id: [] for neg_id in BOUNDARY_NAME_TO_NEG_ID.values()}

    for _, row in boundary_df.iterrows():
        neg_id = int(row["neg_id"])

        if neg_id in boundary_rows:
            boundary_rows[neg_id].append((
                int(row["id"]),
                float(row["x"]),
                float(row["y"]),
                float(row["z"]),
                float(row["orig_length"]) 
            ))

    # Write one file per boundary
    for name, neg_id in BOUNDARY_NAME_TO_NEG_ID.items():
        filename = f"boundary_{name}.ex"
        rows = boundary_rows[neg_id]
        write_ex(filename, rows)
        print(f"Wrote {len(rows):>6} rows -> {filename}")


def write_boundary_ex(cells_df):
    write_boundary_ex_files(cells_df)
