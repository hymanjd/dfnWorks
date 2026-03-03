import pandas as pd

def write_graph_uge(cells_df, conns_df,filename = "graph.uge"):
    """
    Write a UGE-style file with:

        CELLS <Ncells>
        id   x   y   z   volume

        CONNECTIONS <Nconns>
        i    j    xp    yp    zp    ahat_ij

    Where (xp,yp,zp) are the *new* PFLOTRAN faces we computed and ahat_ij
    is the scaled area for that connection.
    """
    with open(filename, "w") as f:
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