import networkx as nx 
import pandas as pd 
from scipy.stats import hmean
def make_connection_data_frame(self, G):
    conns = []
    for u in G.intersections:
        frac_list = list(G.neighbors(u))
        if len(frac_list) == 2:
            frac1 = frac_list[0]
            frac2 = frac_list[1]
            if frac1 != 's' and frac2 != 's' and frac1 != 't' and frac2 != 't': 
                b1 = G.nodes[frac1]['aperture']
                b2 = G.nodes[frac2]['aperture']
                b_hmean = hmean([b1, b2])
                x = G.nodes[u]['x']
                y = G.nodes[u]['y']
                z = G.nodes[u]['z']
                G.nodes[u]['area'] = G.nodes[u]['length'] * b_hmean
                conns.append((frac1, frac2, x, y, z, G.nodes[u]['area']))
    conns_df = pd.DataFrame(conns, columns=["i", "j", "xc", "yc", "zc", "area"])
    #conns_df = conns_df.sort_values(by=["i", "j"]).reset_index(drop=True)
    return conns_df       

def make_cell_data_frame(self, G):
    self.volume = self.surface_area * self.aperture

    for i in range(1, self.num_frac + 1):
        G.nodes[i]['volume'] = self.volume[i-1]
        G.nodes[i]['center'] = self.centers[i-1]

    cells = []
    for i in range(1, self.num_frac + 1):
        cid = i
        x = G.nodes[i]['center'][0]
        y = G.nodes[i]['center'][1]
        z = G.nodes[i]['center'][2]
        volume = G.nodes[i]['volume']
        cells.append((cid, x, y, z, volume))

    cells_df = pd.DataFrame(cells, columns=["id", "x", "y", "z", "volume"])
    return cells_df

def convert_graph_to_data_frames(self, G):

    cells_df = self.make_cell_data_frame(G)
    conns_df = self.make_connection_data_frame(G)

    return cells_df, conns_df