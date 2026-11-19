from .compute_graph_to_pflotran_geometries import compute_graph_to_pflotran_geometries
from .write_graph_uge import write_graph_uge
from .write_boundary_ex import write_boundary_ex


def bigrafv(self,
            G,
            uge_filename="full_mesh_vol_area.uge",
            intersection_list_path="./dfnGen_output/intersection_list.dat",
            selected_boundary_ids=None,
            cells_filename="altered_cells_df.dat",
            properties_filename="dfn_properties.h5",
            boundary_shift=0.01,
            boundary_area_per_dist=1e-2):
    """
    BiGraFV: convert the bipartite graph of a DFN into PFLOTRAN unstructured grid input files.

    Each fracture becomes one PFLOTRAN cell. Each fracture-fracture
    intersection becomes one connection. Every fracture that touches a
    domain boundary gets a small helper cell, and the helper cells are
    written to one ``.ex`` file per domain face.

    The pipeline runs these steps in order:

    1. Read boundary rows from ``intersection_list.dat`` with
       :func:`extract_graph_boundaries_from_dfn`.
    2. Build the cell and connection tables with
       :func:`convert_graph_to_data_frames`.
    3. Place the PFLOTRAN face points and scale the face areas with
       :func:`compute_graph_to_pflotran_geometries`.
    4. Write the six boundary ``.ex`` files with :func:`write_boundary_ex`.
    5. Write the cell table to ``cells_filename``. This file is used to add
       cell attributes to VTK output for visualization.
    6. Write the PFLOTRAN UGE file with :func:`write_graph_uge`.
    7. Write cell permeabilities to ``properties_filename`` with
       :func:`dump_h5_file_for_graph_perm`.

    Parameters
    ----------
    self : object
        DFN Class. Fracture centers, surface areas, apertures and
        permeabilities must be set on it.
    G : networkx.Graph
        Bipartite graph of the DFN, made with
        ``DFN.create_graph("bipartite", ...)``.
    uge_filename : str, optional
        Name of the PFLOTRAN UGE file to write.
        Default is ``"full_mesh_vol_area.uge"``.
    intersection_list_path : str, optional
        Path to the dfnWorks ``intersection_list.dat`` file.
        Default is ``"./dfnGen_output/intersection_list.dat"``.
    selected_boundary_ids : list of int or None, optional
        Boundary face ids to keep, e.g. ``[-3, -5]`` for the left and right
        faces only. If ``None``, all six faces are kept. Default is ``None``.
    cells_filename : str, optional
        Name of the text file the cell table is written to.
        Default is ``"altered_cells_df.dat"``.
    properties_filename : str, optional
        Name of the HDF5 file the cell permeabilities are written to.
        Default is ``"dfn_properties.h5"``.
    boundary_shift : float, optional
        How far each ``.ex`` face point is moved outward past its boundary
        helper cell, as a fraction of the distance from the parent cell
        center. Default is ``0.01``.
    boundary_area_per_dist : float, optional
        Face area per unit distance [m] for the ``.ex`` faces. Each face
        gets area ``boundary_area_per_dist * d``, where ``d`` is the
        distance from the helper cell to its face point. Default is
        ``1e-2``.

    Returns
    -------
    cells_df : pandas.DataFrame
        Cell table, including the boundary helper cells.
    conns_df : pandas.DataFrame
        Connection table with columns ``i``, ``j``, ``xp``, ``yp``, ``zp``
        and ``ahat``.

    Notes
    -----
    All files are written to the current working directory.

    See Also
    --------
    extract_graph_boundaries_from_dfn : Reads boundary rows from ``intersection_list.dat``.
    convert_graph_to_data_frames : Builds the cell and connection tables.
    compute_graph_to_pflotran_geometries : Computes PFLOTRAN face points and areas.
    write_boundary_ex : Writes one ``.ex`` file per domain face.
    write_graph_uge : Writes the PFLOTRAN UGE file.
    dump_h5_file_for_graph_perm : Writes cell permeabilities to HDF5.
    """
    self.print_log("\n--> BiGraFV: Starting")
    boundaries_df = self.extract_graph_boundaries_from_dfn(
        G,
        intersection_list_path=intersection_list_path,
        selected_boundary_ids=selected_boundary_ids)
    self.print_log(f"--> Boundary connections found: {len(boundaries_df)}")

    cells_df, conns_df = self.convert_graph_to_data_frames(G, boundaries_df)
    n_boundary = int((cells_df["row_type"] == "boundary").sum())
    self.print_log(f"--> Cells: {len(cells_df)} "
                   f"({len(cells_df) - n_boundary} fracture, {n_boundary} boundary)")
    self.print_log(f"--> Connections: {len(conns_df)}")

    conns_df = compute_graph_to_pflotran_geometries(cells_df, conns_df)

    write_boundary_ex(cells_df,
                      omega=boundary_shift,
                      area_per_dist=boundary_area_per_dist)
    self.print_log("--> Wrote boundary files boundary_*.ex")

    cells_df.to_csv(cells_filename, sep=" ", index=False)
    self.print_log(f"--> Wrote {cells_filename}")

    write_graph_uge(cells_df, conns_df, uge_filename)
    self.print_log(f"--> Wrote {uge_filename}")

    self.dump_h5_file_for_graph_perm(cells_df, filename=properties_filename)
    self.print_log(f"--> Wrote {properties_filename}")

    self.print_log("--> BiGraFV: Complete\n")
    return cells_df, conns_df
