def write_graph_uge(cells_df, conns_df, filename="full_mesh_vol_area.uge"):
    """
    Write the cell and connection tables to a PFLOTRAN unstructured grid (``.uge``) file.

    Parameters
    ----------
    cells_df : pandas.DataFrame
        Cell table with columns ``id``, ``x``, ``y``, ``z`` and ``volume``.
    conns_df : pandas.DataFrame
        Connection table from :func:`compute_graph_to_pflotran_geometries`,
        with columns ``i``, ``j``, ``xp``, ``yp``, ``zp`` and ``ahat``.
    filename : str, optional
        Name of the ``.uge`` file to write. An existing file is overwritten.
        Default is ``"full_mesh_vol_area.uge"``.

    Returns
    -------
    None

    Notes
    -----
    File format::

        CELLS       <Ncells>
        <id>  <x>  <y>  <z>  <volume>
        ...
        CONNECTIONS <Nconns>
        <i>  <j>  <xp>  <yp>  <zp>  <ahat>
        ...

    Values are tab separated. Floats are written in scientific notation
    with 16 digits after the decimal point (``:.16E``).
    """
    with open(filename, "w") as f:
        f.write(f"CELLS\t{len(cells_df)}\n")
        for _, r in cells_df.iterrows():
            f.write(f"{int(r.id)}\t"
                    f"{r.x:.16E}\t{r.y:.16E}\t{r.z:.16E}\t"
                    f"{r.volume:.16E}\n")

        f.write(f"CONNECTIONS\t{len(conns_df)}\n")
        for _, r in conns_df.iterrows():
            f.write(f"{int(r.i)}\t{int(r.j)}\t"
                    f"{r.xp:.16E}\t{r.yp:.16E}\t{r.zp:.16E}\t"
                    f"{r.ahat:.16E}\n")
