from .compute_graph_to_pflotran_geometries import compute_graph_to_pflotran_geometries
from .write_graph_uge import write_graph_uge
from .extract_graph_boundaries_from_dfn import extract_graph_boundaries_from_dfn

def graph_to_pflotran(self, G, uge_filename = "graph_dfnTest.uge"):
    cells_df, conns_df = self.convert_graph_to_data_frames(G)
    ## fix connections / cross sectional areas 
    conns_df = compute_graph_to_pflotran_geometries(cells_df, conns_df)
    ## write .uge after reprojecting faces and areas
    write_graph_uge(cells_df, conns_df, "graph_testRun.uge")
    ## write .ex files for boundaries
    boundaries_df = extract_graph_boundaries_from_dfn(area_default = 1e5, intersection_list_path='./dfnGen_output/intersection_list.dat')
    self.dump_h5_file_for_graph_perm(filename = 'permeability.h5')


