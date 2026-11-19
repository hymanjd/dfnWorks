import numpy as np
import pandas as pd

from pydfnworks.general.logging import local_print_log


def compute_graph_to_pflotran_geometries(cells_df, conns_df, eps_t=0.05, tol_rel=1e-6):
    """
    Compute the PFLOTRAN face point and face area for each connection.

    For each connection ``i``-``j``, the graph node ``gn`` at
    ``(xc, yc, zc)`` is projected onto the segment ``[Xi, Xj]`` between the
    two cell centers. The face point ``p`` is placed as follows:

    1. **Normal case**: ``p`` is the orthogonal projection of ``gn`` onto
       ``[Xi, Xj]``.
    2. **Colinear endpoint case**: if ``gn`` sits on ``Xi`` or ``Xj``
       (within tolerance) and is colinear with the segment, ``p`` is the
       midpoint of ``[Xi, Xj]``.
    3. **Near-endpoint case**: if the projection parameter ``t`` is within
       ``eps_t`` of either end (``t <= eps_t`` or ``t >= 1 - eps_t``),
       ``p`` is moved inward to ``t = eps_t`` or ``t = 1 - eps_t``. A face
       point on a cell center makes PFLOTRAN fail with
       ``"Face ... cannot be projected ..."``.

    The face area is scaled by the ratio of the straight-line distance
    between the cell centers to the path length through the graph node,
    ``ahat = area * |Xi - Xj| / (|Xi - gn| + |gn - Xj|)``.

    Parameters
    ----------
    cells_df : pandas.DataFrame
        Cell table with columns ``id``, ``x``, ``y`` and ``z``.
    conns_df : pandas.DataFrame
        Connection table with columns ``i``, ``j``, ``xc``, ``yc``, ``zc``
        and ``area``.
    eps_t : float, optional
        How close, as a fraction of the segment, the face point may get to
        either cell center. Default is ``0.05``.
    tol_rel : float, optional
        Relative tolerance for the colinear endpoint case. ``gn`` counts as
        on an endpoint if its distance to it is less than
        ``tol_rel * max(L, 1.0)``, where ``L = |Xi - Xj|``.
        Default is ``1e-6``.

    Returns
    -------
    pandas.DataFrame
        One row per connection with columns ``i``, ``j``, ``xp``, ``yp``,
        ``zp`` and ``ahat``.

    Notes
    -----
    Logs an error and exits if ``lij = 0`` for a connection (both cell
    centers and the graph node are the same point), since ``ahat`` would
    divide by zero.

    Quantities computed for each connection:

    .. list-table::
       :header-rows: 1
       :widths: 15 85

       * - Symbol
         - Description
       * - ``Xi``, ``Xj``
         - Centers of cells ``i`` and ``j``, from ``cells_df``.
       * - ``gn``
         - Graph node ``(xc, yc, zc)``, from ``conns_df``.
       * - ``l1``
         - ``|Xi - gn|``
       * - ``l2``
         - ``|gn - Xj|``
       * - ``l3``
         - ``|Xi - Xj|``, the PFLOTRAN connection length.
       * - ``lij``
         - ``l1 + l2``
       * - ``ahat``
         - Scaled face area, ``area * l3 / lij``.
       * - ``xp, yp, zp``
         - The PFLOTRAN face point ``p``.
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

        # path lengths through the original graph node
        l1 = np.linalg.norm(gn - Xi)
        l2 = np.linalg.norm(Xj - gn)
        l3 = L
        lij = l1 + l2

        if denom == 0.0:
            p_proj = Xi.copy()
        else:
            t = ((gn - Xi) @ v) / denom
            # orthogonal distance from gn to the line through Xi and Xj
            offset_orth = np.linalg.norm(gn - (Xi + t * v))

            # Case 2: gn is on an endpoint and colinear with the segment
            is_endpoint = (L > 0.0 and min(l1, l2) < tol_rel * max(L, 1.0))
            is_colinear = offset_orth < tol_rel * max(L, 1.0)

            if is_endpoint and is_colinear:
                p_proj = 0.5 * (Xi + Xj)
            elif t <= eps_t:
                # Case 3: keep the face away from the endpoints
                p_proj = Xi + eps_t * v
            elif t >= 1.0 - eps_t:
                p_proj = Xi + (1.0 - eps_t) * v
            else:
                # Case 1
                p_proj = Xi + t * v

        if lij == 0.0:
            local_print_log(
                f"Error. Division by zero computing the face area of connection "
                f"{int(r.i)}-{int(r.j)}: both cell centers and the graph node are the "
                "same point (l1 + l2 = 0).", 'error')
        ahat = float(r.area) * (l3 / lij)

        xp.append(p_proj[0])
        yp.append(p_proj[1])
        zp.append(p_proj[2])
        ahats.append(ahat)

    return pd.DataFrame({
        "i": conns_df["i"].values,
        "j": conns_df["j"].values,
        "xp": xp,
        "yp": yp,
        "zp": zp,
        "ahat": ahats
    })
