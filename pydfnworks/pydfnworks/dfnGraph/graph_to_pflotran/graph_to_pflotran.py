from .compute_graph_to_pflotran_geometries import compute_graph_to_pflotran_geometries
from .write_graph_uge import write_graph_uge
from .extract_graph_boundaries_from_dfn import extract_graph_boundaries_from_dfn
from .write_boundary_ex import write_boundary_ex,write_boundary_ex_files, write_ex

def graph_to_pflotran(self, G, uge_filename = "graph_dfnTest.uge"):
    ## write .ex files for boundaries and make pd dataframe of boundary connections
    boundaries_df = extract_graph_boundaries_from_dfn(area_default = 1e5, intersection_list_path='./dfnGen_output/intersection_list.dat')
    cells_df, conns_df = self.convert_graph_to_data_frames(G, boundaries_df)
    ## fix connections / cross sectional areas 
    conns_df = compute_graph_to_pflotran_geometries(cells_df, conns_df)
    write_boundary_ex(cells_df)
    ## write .uge after reprojecting faces and areas
    write_graph_uge(cells_df, conns_df, "graph_testRun.uge")
    self.dump_h5_file_for_graph_perm(filename = 'permeability.h5')


