import numpy as np
import pandas as pd


def compute_graph_to_pflotran_geometries(
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
    ahats = []
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
            
        if denom == 0.0:
            p_proj = Xi.copy()
            offset_orth = np.linalg.norm(gn - p_proj)
        else:
            t = ((gn - Xi) @ v) / denom
            # orthogonal distance to the infinite line
            p_proj = Xi + t * v
            offset_orth = np.linalg.norm(gn - p_proj)

            # detect "colinear endpoint" case (Case 2)
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
        ahats.append(ahat)

    return pd.DataFrame({
        "i": conns_df["i"].values,
        "j": conns_df["j"].values,
        "xp": xp,
        "yp": yp,
        "zp": zp,
        "ahat": ahats
    })














