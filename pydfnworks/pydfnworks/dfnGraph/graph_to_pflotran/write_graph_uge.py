import pandas as pd

def write_graph_uge(cells_df, conns_df,filename = "full_mesh_vol_area.uge"):
    """
    Write cell and connection data to a PFLOTRAN unstructured grid (``.uge``) file.

    Parameters
    ----------
    cells_df : pd.DataFrame
        Cell DataFrame with one row per cell.
        Required columns: ``id``, ``x``, ``y``, ``z``, ``volume``.
    conns_df : pd.DataFrame
        Connection DataFrame with one row per connection.
        Required columns: ``i``, ``j``, ``xp``, ``yp``, ``zp``, ``area``.
    filename : str, optional
        Path and name of the ``.uge`` file to create or overwrite.
        Default is ``"full_mesh_vol_area.uge"``.

    Returns
    -------
    None
        Writes directly to *filename*; no value is returned.

    Notes
    -----
    Output file format::

        CELLS       <Ncells>
        <id>  <x>  <y>  <z>  <volume>
        ...
        CONNECTIONS <Nconns>
        <i>  <j>  <xp>  <yp>  <zp>  <area>
        ...

    All floating-point values are written in 16-digit scientific notation
    (``:.16E``). ``xp``, ``yp``, ``zp`` are the PFLOTRAN face-point
    coordinates and ``area`` is the scaled connection area, both computed
    by the face-geometry step prior to calling this function.
    """
    with open(filename, "w") as f:
        f.write(f"CELLS\t{len(cells_df)}\n")
    
        for _, r in cells_df.iterrows():
            f.write(
                f"{int(r.id)}\t"
                f"{r.x:.16E}\t{r.y:.16E}\t{r.z:.16E}\t"
                f"{r.volume:.16E}\n"
            )

        f.write(f"CONNECTIONS\t{len(conns_df)}\n")
        for _, r in conns_df.iterrows():
            f.write(
                f"{int(r.i)}\t{int(r.j)}\t"
                f"{r.xp:.16E}\t{r.yp:.16E}\t{r.zp:.16E}\t"
                f"{r.ahat:.16E}\n"
            )
