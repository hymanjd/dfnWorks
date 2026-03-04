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

def write_boundary_ex(
    intersection_list_path,
    area_default,
):
    """
    Reads intersection_list.dat and writes one .ex file per boundary face.

    intersection_list.dat format:
        f1 f2 x y z length
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

    # Build map dict: neg_id -> list of (node_i, x, y, z, area)
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

    for name, neg_id in BOUNDARY_NAME_TO_NEG_ID.items():
        filename = f"boundary_{name}.ex"
        rows = boundary_rows[neg_id]
        write_ex(filename, rows)
        print(f"Wrote {len(rows):>6} rows -> {filename}")


    return pd.concat([
        pd.DataFrame(rows, columns=["node_id", "x", "y", "z", "area"])
        .assign(neg_id=neg_id)
        for neg_id, rows in boundary_rows.items()
    ])
    

def extract_graph_boundaries_from_dfn(area_default, intersection_list_path='./dfnGen_output/intersection_list.dat'):
    boundaries_df = write_boundary_ex(
        intersection_list_path=intersection_list_path,
        area_default=area_default
    )
    print(boundaries_df)

    return boundaries_df
