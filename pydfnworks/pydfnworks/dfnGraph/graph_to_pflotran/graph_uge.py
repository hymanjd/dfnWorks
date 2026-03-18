import networkx as nx 
import pandas as pd 
import numpy as np
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
                try:
                    b2 = G.nodes[frac2]['aperture']
                    b_hmean = hmean([b1, b2])
                except KeyError:
                    b2 = None  # or some default value
                    b_hmean = b1
                # b2 = G.nodes[frac2]['aperture']
                # b_hmean = hmean([b1, b2])
                x = G.nodes[u]['x']
                y = G.nodes[u]['y']
                z = G.nodes[u]['z']
                G.nodes[u]['area'] = G.nodes[u]['length'] * b_hmean
                conns.append((frac1, frac2, x, y, z, G.nodes[u]['area']))
    conns_df = pd.DataFrame(conns, columns=["i", "j", "xc", "yc", "zc", "area"])
    #conns_df = conns_df.sort_values(by=["i", "j"]).reset_index(drop=True)
    return conns_df     


def add_boundary_nodes_as_conns(self, G, altered_cells_df, conns_df, boundaries_df):
    boundaries_df_copy = boundaries_df.copy()
    altered_cells_df_copy = altered_cells_df.copy()
    conns_df_copy = conns_df.copy()
    # turn length into area by multiplying by fracture aperture
    for i in range(1, self.num_frac + 1):
        if i in boundaries_df_copy['id'].values:  
            X1 = np.array(altered_cells_df_copy.loc[altered_cells_df_copy['id'] == i, ['x', 'y', 'z']].values[0])
            # X1_vol = np.array(altered_cells_df_copy.loc[altered_cells_df_copy['id'] == i, ['volume']].values[0])
            # find all boundary rows for this fracture and calculate L for each
            mask = boundaries_df_copy['id'] == i
            L_values = boundaries_df_copy.loc[mask, ['x', 'y', 'z']].apply(
                lambda row: np.linalg.norm(np.array(row) - X1), axis=1
            )

            #P(t) = X1 + t * (X0 - X1) - currently 0.5 is for the halfway point
            total_L = L_values.sum() #unused
            weights = L_values / total_L #unused

            midpoints = boundaries_df_copy.loc[mask, ['x', 'y', 'z']].apply(
                lambda row: X1 + 0.1 * (np.array(row) - X1), axis=1, result_type='expand'
            )
            midpoints.columns = ['x', 'y', 'z']
            boundaries_df_copy.loc[mask, 'x'] = midpoints['x'].values
            boundaries_df_copy.loc[mask, 'y'] = midpoints['y'].values
            boundaries_df_copy.loc[mask, 'z'] = midpoints['z'].values

            boundaries_df_copy.loc[mask, "length"] *= G.nodes[i]['aperture']
            # boundaries_df_copy.rename(columns={'length': 'area'})


    # match neg_ids to get the altered_cells id
    merged = boundaries_df_copy.merge(altered_cells_df_copy[['id', 'neg_id']], on='neg_id')
    # build and concat to conns_df
    conns_df_updated = pd.concat([
        conns_df_copy,
        merged.rename(columns={'id_y': 'i', 'id_x': 'j', 'x': 'xc', 'y': 'yc', 'z': 'zc', 'length': 'area'})[['i', 'j', 'xc', 'yc', 'zc', 'area']]
    ]).reset_index(drop=True)
    print(f'Total area: {conns_df_updated["area"].sum():.12f}')
    return conns_df_updated




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

def add_boundary_nodes_as_cells(self, G, cells_df, boundaries_df, omega=0.01):
    import numpy as np
    import pandas as pd

    boundaries_df_copy = boundaries_df.copy()
    cells_df_copy = cells_df.copy()
    
    # Preserve metadata on original cells but make all boundarys nan
    cells_df_copy["row_type"] = "internal"
    cells_df_copy["parent_cell_id"] = np.nan
    cells_df_copy["boundary_index"] = np.nan
    cells_df_copy["orig_length"] = np.nan

    # Preserve metadata on boundary rows before modifying length 
    boundaries_df_copy["row_type"] = "boundary"
    boundaries_df_copy["parent_cell_id"] = boundaries_df_copy["id"]
    boundaries_df_copy["boundary_index"] = boundaries_df_copy["neg_id"]
    boundaries_df_copy["orig_length"] = boundaries_df_copy["length"]

    move_fraction = omega

    cell_xyz_map = cells_df_copy.set_index("id")[["x", "y", "z"]].to_dict("index")
    cell_volume_map = cells_df_copy.set_index("id")["volume"].to_dict()

    for i in boundaries_df_copy["id"].unique():
    # calculate scaling factor (omega) for dividing cell volume based on coordinate of cell center and boundary

        xyz = cell_xyz_map[i]
        X1 = np.array([xyz["x"], xyz["y"], xyz["z"]], dtype=float)
        # find all boundary rows for this fracture and calculate L for each
        mask = boundaries_df_copy["id"] == i
        boundary_xyz = boundaries_df_copy.loc[mask, ["x", "y", "z"]].to_numpy(dtype=float)

        L_values = np.linalg.norm(boundary_xyz - X1, axis=1)

        volume_check_og = cell_volume_map[i]
        print(volume_check_og)
        # change volumes to share and preserve total volume
        volume_partition = 1 / (len(L_values) + 1)
        new_volume = volume_check_og * (1 - volume_partition)

        cells_df_copy.loc[cells_df_copy["id"] == i, "volume"] = new_volume
        # each boundary gets an equal share of the removed volume
        boundaries_df_copy.loc[mask, "length"] = (
            volume_partition * float(volume_check_og) / len(L_values)
        )
        # check: cell + all boundaries should equal original
        new_total = new_volume + boundaries_df_copy.loc[mask, "length"].sum()
        print(f"Volume difference for {i} = {volume_check_og - new_total:.6f}")

    altered_cells_df = pd.concat([
        cells_df_copy,
        boundaries_df_copy.rename(columns={'length': 'volume'})
    ]).reset_index(drop=True)

    altered_cells_df['id'] = altered_cells_df.index + 1
    print('Altered cells df ######################')
    print(altered_cells_df)
    return altered_cells_df

def convert_graph_to_data_frames(self, G, boundaries_df):

    cells_df = self.make_cell_data_frame(G)
    conns_df = self.make_connection_data_frame(G)
    altered_cells_df = self.add_boundary_nodes_as_cells(G, cells_df, boundaries_df)
    conns_df_updated = self.add_boundary_nodes_as_conns(G, altered_cells_df, conns_df, boundaries_df)
    # return cells_df, conns_df
    # return cells and connections after updating for nodes and connections at the boundaries based on intersection_list.dat
    return altered_cells_df, conns_df_updated