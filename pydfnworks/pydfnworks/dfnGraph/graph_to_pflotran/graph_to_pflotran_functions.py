import numpy as np
import pandas as pd
import os

def extract_boundary_string(fragment: str, key: str) -> str | None:
    """
    Extract kwarg value for patterns like:
      inflow = 'left'
      inflow='left'
      inflow = "left"
    Returns the string value without quotes, or None if not found.
    """
    idx = fragment.find(key)
    if idx < 0:
        return None

    # start at key, find '='
    eq = fragment.find("=", idx)
    if eq < 0:
        return None

    s = fragment[eq + 1 :].lstrip()  # after '='
    if not s:
        return None

    q = s[0]
    if q not in ("'", '"'):
        return None

    endq = s.find(q, 1)
    if endq < 0:
        return None

    return s[1:endq].strip()


def write_boundary_ex(
    driver_path: str,
    intersection_list_path: str,
    inflow_ex_path: str = "inflow_nodes.ex",
    outflow_ex_path: str = "outflow_nodes.ex",
    area_default: float = 99999.99999999999,
):
    """
    Read driver.py to get inflow/outflow strings, then write
    inflow_nodes.ex and outflow_nodes.ex by filtering intersection_list.dat.

    intersection_list.dat format:
      f1 f2 x y z length
    where boundary rows have one negative boundary ID (-1..-6) and one positive node ID.
    """
    
    BOUNDARY_NAME_TO_NEG_ID = {
        "top": -1,
        "bottom": -2,
        "left": -3,
        "front": -4,
        "right": -5,
        "back": -6,
    }

    
    # read driver and find the create_graph line
    with open(driver_path, "r", errors="replace") as f:
        lines = f.readlines()

    call_text = None
    for i, line in enumerate(lines):
        if "DFN.create_graph" in line:
            chunk = line
            j = i + 1
            while ")" not in chunk and j < len(lines):
                chunk += lines[j]
                j += 1
            call_text = chunk
            break

    if call_text is None:
        raise ValueError(f"Could not find a line containing 'DFN.create_graph' in {driver_path}")

    inflow = extract_boundary_string(call_text, "inflow")
    outflow = extract_boundary_string(call_text, "outflow")

    if inflow is None or outflow is None:
        raise ValueError(
            f"Found DFN.create_graph call but could not parse inflow/outflow in {driver_path}.\n"
            f"Call text was:\n{call_text}"
        )

    inflow = inflow.lower()
    outflow = outflow.lower()

    if inflow not in BOUNDARY_NAME_TO_NEG_ID:
        raise ValueError(f"Unsupported inflow='{inflow}'. Must be one of {list(BOUNDARY_NAME_TO_NEG_ID)}")
    if outflow not in BOUNDARY_NAME_TO_NEG_ID:
        raise ValueError(f"Unsupported outflow='{outflow}'. Must be one of {list(BOUNDARY_NAME_TO_NEG_ID)}")

    inflow_neg = BOUNDARY_NAME_TO_NEG_ID[inflow]
    outflow_neg = BOUNDARY_NAME_TO_NEG_ID[outflow]

    # scan intersection_list.dat and collect boundary rows
    inflow_rows = []
    outflow_rows = []

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
                a = int(parts[0])
                b = int(parts[1])
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
            except ValueError:
                continue

            # boundary rows: one negative boundary id and one positive node id
            if a < 0 and b > 0:
                neg_id, node_i = a, b
            elif b < 0 and a > 0:
                neg_id, node_i = b, a
            else:
                continue

            row_out = (node_i, x, y, z, area_default)

            if neg_id == inflow_neg:
                inflow_rows.append(row_out)
            if neg_id == outflow_neg:
                outflow_rows.append(row_out)

    # write .ex files
    def write_boundary_fun(path: str, rows):
        with open(path, "w") as g:
            g.write(f"CONNECTIONS\t{len(rows)}\n")
            for node_i, x, y, z, area in rows:
                g.write(f"{int(node_i)}\t{x: .12E}\t{y: .12E}\t{z: .12E}\t{area: .14g}\n")

    write_boundary_fun(inflow_ex_path, inflow_rows)
    write_boundary_fun(outflow_ex_path, outflow_rows)

    return {
        "driver_path": driver_path,
        "inflow": inflow,
        "outflow": outflow,
        "inflow_neg_id": inflow_neg,
        "outflow_neg_id": outflow_neg,
        "n_inflow": len(inflow_rows),
        "n_outflow": len(outflow_rows),
        "inflow_ex_path": inflow_ex_path,
        "outflow_ex_path": outflow_ex_path,
    }

def parse_uge(path):
    """
    Parse a UGE file with blocks:

        CELLS N
        id  x  y  z  volume

        CONNECTIONS M
        i  j  xc  yc  zc  area

    Returns:
    cells_df : DataFrame[id, x, y, z, volume]
    conns_df : DataFrame[i, j, xc, yc, zc, area]
    """
    cells, conns, mode = [], [], None
    with open(path, "r") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            u = s.upper()
            if u.startswith("CELLS"):
                mode = "cells"
                continue
            if u.startswith("CONNECTIONS"):
                mode = "conns"
                continue
            parts = s.split()
            if mode == "cells" and len(parts) >= 5:
                cid = int(parts[0])
                x, y, z, vol = map(float, parts[1:5])
                cells.append((cid, x, y, z, vol))
            elif mode == "conns" and len(parts) >= 6:
                i, j = map(int, parts[:2])
                xc, yc, zc, area = map(float, parts[2:6])
                conns.append((i, j, xc, yc, zc, area))

    cells_df = pd.DataFrame(cells, columns=["id", "x", "y", "z", "volume"])
    conns_df = pd.DataFrame(conns, columns=["i", "j", "xc", "yc", "zc", "area"])
    return cells_df, conns_df


def parse_boundary_ex(path):
    """
    Parse inflow/outflow .ex file with a CONNECTIONS section.
    Returns empty DataFrame if file is missing.
    """
    if path is None or not os.path.exists(path):
        return pd.DataFrame(columns=["cell", "x", "y", "z", "area"])

    rows, mode = [], None
    with open(path, "r") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.upper().startswith("CONNECTIONS"):
                mode = "conns"
                continue
            if mode == "conns":
                p = s.split()
                if len(p) >= 5:
                    rows.append(
                        (
                            int(p[0]),
                            float(p[1]),
                            float(p[2]),
                            float(p[3]),
                            float(p[4]),
                        )
                    )
    return pd.DataFrame(rows, columns=["cell", "x", "y", "z", "area"])


def add_projection_and_lengths(
    cells_df,
    conns_df,
    tol_rel=1e-12,
):
    """
    For each connection i-j with graph node gn:

        Xi = cell center i
        Xj = cell center j
        gn = original graph node (xc,yc,zc)

    New PFLOTRAN face p is:

        - Normally: orthogonal projection of gn onto segment [Xi, Xj]
        - Special case: if gn lies (within tolerance) on Xi or Xj and
          is colinear with Xi-Xj, then p = midpoint(Xi, Xj).
        - Case 3 near endpoints OR outside segment. This avoids p being equal (or extremely close) to Xi or Xj, which
          can trigger PFLOTRAN "Face ... cannot be projected ..." failures

    We compute and add columns:

        l1    = |Xi - gn|
        l2    = |gn - Xj|
        l3    = |Xi - Xj|        (PFLOTRAN length)
        lij   = l1 + l2
        eps   = l3 / lij         (if lij > 0)
        ahat  = eps * area       (area from graph.uge is b_h * lij)
        offset = |gn - p|        (distance from original node to PFLOTRAN face)

        xp, yp, zp : coordinates of p
    """
    cent = {
        int(r.id): np.array([r.x, r.y, r.z], float)
        for _, r in cells_df.iterrows()
    }

    xp, yp, zp = [], [], []
    l1s, l2s, l3s, lijs, eps_list, ahats, offsets = [], [], [], [], [], [], []

    for _, r in conns_df.iterrows():
        Xi = cent[int(r.i)]
        Xj = cent[int(r.j)]
        gn = np.array([r.xc, r.yc, r.zc], float)

        v = Xj - Xi
        L = np.linalg.norm(v)
        denom = float(v @ v)

        # geometric lengths based on original graph node
        l1 = np.linalg.norm(gn - Xi)
        l2 = np.linalg.norm(Xj - gn)
        l3 = np.linalg.norm(Xj - Xi)
        lij = l1 + l2

        # # default: orthogonal projection of gn onto Xi-Xj
        # if denom == 0.0:
        #     p_proj = Xi.copy()
        #     offset_orth = np.linalg.norm(gn - p_proj)
        # else:
        #     t = ((gn - Xi) @ v) / denom
        #     t = max(0.0, min(1.0, t))  # clamp to [0,1]
        #     p_proj = Xi + t * v
        #     offset_orth = np.linalg.norm(gn - p_proj)

        # # detect special "colinear endpoint" case:
        # # gn almost on Xi or Xj, AND very close to the line
        # is_endpoint = (
        #     L > 0.0
        #     and min(np.linalg.norm(gn - Xi), np.linalg.norm(gn - Xj))
        #     < tol_rel * max(L, 1.0)
        # )
        # is_colinear = offset_orth < tol_rel * max(L, 1.0)

        # if is_endpoint and is_colinear:
        #     # Use midpoint as new face (PFLOTRAN face between cells)
        #     p = 0.5 * (Xi + Xj)
        # else:
        #     # Use true orthogonal projection
        #     p = p_proj
            
            
            
        if denom == 0.0:
            p_proj = Xi.copy()
            offset_orth = np.linalg.norm(gn - p_proj)
        else:
            t = ((gn - Xi) @ v) / denom
            # orthogonal distance to the infinite line (using unclamped t)
            p_proj = Xi + t * v
            offset_orth = np.linalg.norm(gn - p_proj)

            # detect special "colinear endpoint" case (Case 2)
            is_endpoint = (
                L > 0.0
                and min(np.linalg.norm(gn - Xi), np.linalg.norm(gn - Xj))
                < tol_rel * max(L, 1.0)
            )
            is_colinear = offset_orth < tol_rel * max(L, 1.0)

            frac = 0.05      # 5% buffer inside the segment
            eps_t = frac     # treat within 5% of an endpoint as unsafe

            if is_endpoint and is_colinear:
                # midpoint
                p_proj = 0.5 * (Xi + Xj)

            else:
                # Case 3: keep face away from endpoints
                if t <= eps_t:
                    p_proj = Xi + frac * v
                elif t >= 1.0 - eps_t:
                    p_proj = Xi + (1.0 - frac) * v
                else:
                    # Safe interior projection
                    p_proj = Xi + t * v

            
        offset = np.linalg.norm(gn - p_proj)

        # epsilon and ahat using area from graph.uge
        area = float(r.area)
        if lij > 0.0:
            eps = l3 / lij
            ahat = eps * area
        else:
            eps = 0.0
            # ahat = 0
            ahat = area

        xp.append(p_proj[0]); yp.append(p_proj[1]); zp.append(p_proj[2])
        l1s.append(l1); l2s.append(l2); l3s.append(l3)
        lijs.append(lij); eps_list.append(eps); ahats.append(ahat)
        offsets.append(offset)

    out = conns_df.copy()
    out["xp"] = xp
    out["yp"] = yp
    out["zp"] = zp
    out["l1"] = l1s
    out["l2"] = l2s
    out["l3"] = l3s
    out["lij"] = lijs
    out["eps"] = eps_list
    out["ahat"] = ahats
    out["offset"] = offsets
    return out


def write_pflotran_uge(path, cells_df, conns_df):
    """
    Write a UGE-style file with:

        CELLS <Ncells>
        id   x   y   z   volume

        CONNECTIONS <Nconns>
        i    j    xp    yp    zp    ahat_ij

    Where (xp,yp,zp) are the *new* PFLOTRAN faces we computed and ahat_ij
    is the scaled area for that connection.
    """
    with open(path, "w") as f:
        f.write(f"CELLS\t{len(cells_df)}\n")
        for _, r in cells_df.iterrows():
            f.write(
                f"{int(r.id)}\t"
                f"{r.x: .12E}\t{r.y: .12E}\t{r.z: .12E}\t"
                f"{r.volume: .12E}\n"
            )

        f.write(f"CONNECTIONS\t{len(conns_df)}\n")
        for _, r in conns_df.iterrows():
            f.write(
                f"{int(r.i)}\t{int(r.j)}\t"
                f"{r.xp: .12E}\t{r.yp: .12E}\t{r.zp: .12E}\t"
                f"{r.ahat: .16E}\n"
            )


def convert_graph_to_pflotran(
    driver_path,
    in_uge_path,
    out_uge_path,
    inflow_path=None,
    outflow_path=None,
    tol_rel=1e-6,):
    """
    Main function call for converting bipartite dfnWorks graph into PFLOTRAN usable format:

    1. Parse graph UGE.
    2. Compute PFLOTRAN faces and ahat_ij.
    3. Write PFLOTRAN UGE.

    Parameters:
    in_uge_path : str
        Path to input graph.uge
    out_uge_path : str
        Path to output pflotran .uge
    inflow_path, outflow_path : str or None
        Currently unused in the conversion, but parsed for convenience
        and possible future use.
    tol_rel : float
        Relative tolerance for detecting colinear/endpoint case.

    Returns:
        .uge file for PFLOTRAN and optionally cells_df, conns_df : DataFrames
        Cells and connections with added PFLOTRAN geometry/columns.
    """
    # Run boundary condition function
    info = write_boundary_ex(
        driver_path=driver_path,
        intersection_list_path='./dfnGen_output/intersection_list.dat',
        inflow_ex_path=inflow_path,
        outflow_ex_path=outflow_path,
        area_default=99999.99999999999,
    )

    print("\nBoundary extraction complete:")
    print(f"  inflow  = {info['inflow']}  (neg id {info['inflow_neg_id']})")
    print(f"  outflow = {info['outflow']} (neg id {info['outflow_neg_id']})")
    print(f"  inflow connections : {info['n_inflow']}")
    print(f"  outflow connections: {info['n_outflow']}")
    print(f"  wrote {info['inflow_ex_path']}")
    print(f"  wrote {info['outflow_ex_path']}")


    # Now continue writing UGE for PFLOTRAN 
    cells_df, conns_df = parse_uge(in_uge_path)

    # parse boundaries if provided (even if not used in conversion yet)
    _ = parse_boundary_ex(inflow_path)
    _ = parse_boundary_ex(outflow_path)

    conns_df = add_projection_and_lengths(
        cells_df, conns_df, tol_rel=tol_rel
    )
    write_pflotran_uge(out_uge_path, cells_df, conns_df)
    return cells_df, conns_df


#if __name__ == "__main__":
#    # Example standalone usage
#    cells_df, conns_df = convert_graph_to_pflotran(
#        "graph.uge",
#        "pflotran_connection_lengths.uge",
#        inflow_path="inflow_nodes.ex",
#        outflow_path="outflow_nodes.ex",
#    )
#    print("Wrote pflotran_connection_lengths.uge")
