from .compute_graph_to_pflotran_geometries import compute_graph_to_pflotran_geometries
from .write_graph_uge import write_graph_uge
from .write_boundary_ex import write_boundary_ex,write_boundary_ex_files, write_ex

def graph_to_pflotran(self, G, uge_filename = "full_mesh_vol_area.uge"):
    """
    Convert a dfnWorks-generated DFN into a PFLOTRAN unstructured grid.

    Orchestrates the full pipeline from graph representation to PFLOTRAN input files in seven sequential steps:

    1. Extract boundary information via :func:`extract_graph_boundaries_from_dfn`.
    2. Build cell and connection DataFrames via :func:`convert_graph_to_data_frames`.
    3. Compute orthogonal PFLOTRAN connection face geometry via
       :func:`compute_graph_to_pflotran_geometries`.
    4. Write boundary condition ``.ex`` files via :func:`write_boundary_ex`.
    5. Export the cell DataFrame to ``altered_cells_df.dat`` for debugging.
    6. Write the PFLOTRAN UGE mesh file via :func:`write_graph_uge`.
    7. Export permeability-related graph properties to HDF5 via
       :func:`dump_h5_file_for_graph_perm`.

    Parameters
    ----------
    G : networkx.Graph
        NetworkX graph representation of the DFN network.
    uge_filename : str, optional
        Name of the PFLOTRAN UGE output file.
        Default is ``"full_mesh_vol_area.uge"``.

    Returns
    -------
    None
        All output is written to disk; no value is returned.

    Raises
    ------
    FileNotFoundError
        If required DFN input files (e.g. ``intersection_list.dat``) are
        not found.
    ValueError
        If graph boundary extraction or DataFrame conversion fails due to
        missing graph attributes.

    See Also
    --------
    extract_graph_boundaries_from_dfn : Parses boundary rows from ``intersection_list.dat``.
    convert_graph_to_data_frames : Builds cell and connection DataFrames from the graph.
    compute_graph_to_pflotran_geometries : Computes PFLOTRAN face-point geometry.
    write_boundary_ex : Writes one ``.ex`` file per boundary face.
    write_graph_uge : Writes the PFLOTRAN UGE mesh file.
    dump_h5_file_for_graph_perm : Exports graph permeability properties to HDF5.
    """
    boundaries_df = self.extract_graph_boundaries_from_dfn(
        G,
        intersection_list_path='./dfnGen_output/intersection_list.dat'
        #selected_boundary_ids=[-3, -5]
    )
    cells_df, conns_df = self.convert_graph_to_data_frames(G, boundaries_df)
    conns_df = compute_graph_to_pflotran_geometries(cells_df, conns_df)
    write_boundary_ex(cells_df)
    cells_df.to_csv("altered_cells_df.dat", sep=" ", index=False) #this is used for writing attributes to VTK file for visualization
    write_graph_uge(cells_df, conns_df, uge_filename)
    self.dump_h5_file_for_graph_perm(cells_df,filename = 'dfn_properties.h5')

    return cells_df, conns_df

