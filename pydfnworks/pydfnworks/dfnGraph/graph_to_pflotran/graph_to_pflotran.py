from .compute_graph_to_pflotran_geometries import compute_graph_to_pflotran_geometries
from .write_graph_uge import write_graph_uge
# from .extract_graph_boundaries_from_dfn import extract_graph_boundaries_from_dfn
from .write_boundary_ex import write_boundary_ex,write_boundary_ex_files, write_ex

def graph_to_pflotran(self, G, uge_filename = "full_mesh_vol_area.uge"):
    ## write .ex files for boundaries and make pd dataframe of boundary connections
    boundaries_df = self.extract_graph_boundaries_from_dfn(G, area_default = 1e5, intersection_list_path='./dfnGen_output/intersection_list.dat')
    cells_df, conns_df = self.convert_graph_to_data_frames(G, boundaries_df)
    ## fix connections / cross sectional areas 
    conns_df = compute_graph_to_pflotran_geometries(cells_df, conns_df)
    write_boundary_ex(cells_df)
    ## write .uge after reprojecting faces and areas
    write_graph_uge(cells_df, conns_df, "full_mesh_vol_area.uge")
    self.dump_h5_file_for_graph_perm(cells_df,filename = 'dfn_properties.h5')


