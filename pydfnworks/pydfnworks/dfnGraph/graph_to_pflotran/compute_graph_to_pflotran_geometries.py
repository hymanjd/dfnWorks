import numpy as np
import pandas as pd


def compute_graph_to_pflotran_geometries(
    cells_df,
    conns_df,
    eps_t=0.05,
    tol_rel=1e-6
):
    """
    Compute PFLOTRAN face geometry for each connection in the unstructured grid.

    For each connection i-j, the original graph node ``gn`` at ``(xc, yc, zc)``
    is projected onto the segment ``[Xi, Xj]`` to produce a PFLOTRAN face point
    ``p``. Three cases determine how ``p`` is placed:

    1. **Normal case**: ``p`` is the orthogonal projection of ``gn`` onto
    ``[Xi, Xj]``.
    2. **Colinear endpoint case**: If ``gn`` lies on (or within tolerance of)
    either endpoint ``Xi`` or ``Xj`` and is colinear with the segment, ``p``
    is set to the midpoint of ``[Xi, Xj]`` to avoid degenerate geometry.
    3. **Near-endpoint case**: If the projected parameter ``t`` falls within
    ``eps_t`` of either end of the segment (``t <= eps_t`` or
    ``t >= 1 - eps_t``), ``p`` is clamped inward to avoid coinciding with
    ``Xi`` or ``Xj``, which would trigger PFLOTRAN
    ``"Face ... cannot be projected ..."`` failures.

    Parameters
    ----------
    cells_df : pd.DataFrame
        Cell DataFrame containing cell IDs and centre coordinates.
        Required columns: ``id``, ``x``, ``y``, ``z``.
    conns_df : pd.DataFrame
        Connection DataFrame containing pairs of connected cell IDs, the
        original graph node coordinates, and face area.
        Required columns: ``i``, ``j``, ``xc``, ``yc``, ``zc``, ``area``.
    eps_t : float, optional
        Fractional tolerance controlling how close the projected face point
        ``p`` may be to either endpoint of ``[Xi, Xj]``. If the projection
        parameter ``t`` falls within ``eps_t`` of 0 or 1, ``p`` is clamped
        inward to ``eps_t`` or ``1 - eps_t`` respectively. Default is
        ``0.05``.
    tol_rel : float, optional
        Relative tolerance for detecting the colinear endpoint case (Case 2).
        A graph node is considered coincident with an endpoint if its distance
        to that endpoint is less than ``tol_rel * max(L, 1.0)``, where
        ``L = |Xi - Xj|``. Default is ``1e-6``.

    Returns
    -------
    pd.DataFrame
        One row per connection with columns ``i``, ``j``, ``xp``, ``yp``,
        ``zp``, and ``ahat``.

    Notes
    -----
    The following intermediate quantities are computed per connection:

    .. list-table::
    :header-rows: 1
    :widths: 10 90

    * - Symbol
        - Description
    * - ``Xi``
        - Centre coordinates of cell ``i``, taken from ``cells_df``.
    * - ``Xj``
        - Centre coordinates of cell ``j``, taken from ``cells_df``.
    * - ``gn``
        - Original graph node coordinates ``(xc, yc, zc)`` from ``conns_df``.
    * - ``l1``
        - Distance from ``Xi`` to ``gn``, i.e. ``|Xi - gn|``.
    * - ``l2``
        - Distance from ``gn`` to ``Xj``, i.e. ``|gn - Xj|``.
    * - ``l3``
        - Distance ``|Xi - Xj|``, used as the PFLOTRAN connection length.
    * - ``lij``
        - Sum ``l1 + l2``.
    * - ``eps``
        - Length ratio ``l3 / lij`` (computed when ``lij > 0``).
    * - ``ahat``
        - Scaled face area ``eps * area``.
    * - ``offset``
        - Distance from the original graph node to the face point, ``|gn - p|``.
    * - ``xp, yp, zp``
        - Coordinates of the PFLOTRAN face point ``p``.
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

            if is_endpoint and is_colinear:
                # midpoint
                p_proj = 0.5 * (Xi + Xj)
                t = 0.0
            else:
                # Case 3: keep face away from endpoints
                if t <= eps_t:
                    p_proj = Xi + eps_t * v
                elif t >= 1.0 - eps_t:
                    p_proj = Xi + (1.0 - eps_t) * v
                else:
                    p_proj = Xi + t * v

            
        offset = np.linalg.norm(gn - p_proj)

        # epsilon and ahat using area from graph.uge
        area = float(r.area)
        eps = l3 / lij
        ahat = area * eps

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














