from .compute_graph_to_pflotran_geometries import compute_graph_to_pflotran_geometries
from .write_graph_uge import write_graph_uge

def driver_for_graph_files_for_pflotran(self, G, uge_filename = "graph_dfnTest.uge"):
    cells_df, conns_df = self.convert_graph_to_data_frames(G)
    ## fix connections / cross sectional areas 
    conns_df = compute_graph_to_pflotran_geometries(cells_df, conns_df)
    write_graph_uge(cells_df, conns_df, "graph_testRun.uge")
    # dump_graph_ex(G, boundary = "inflow")
    # dump_graph_ex(G, boundary = "outflow")
    #self.dump_boundary_files()
