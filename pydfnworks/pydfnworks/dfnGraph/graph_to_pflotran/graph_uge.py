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
    # preserve original geometry for later association
    boundaries_df_copy['x_orig'] = boundaries_df_copy['x']
    boundaries_df_copy['y_orig'] = boundaries_df_copy['y']
    boundaries_df_copy['z_orig'] = boundaries_df_copy['z']
    boundaries_df_copy['length_orig'] = boundaries_df_copy['length']
    # turn length into area by multiplying by fracture aperture
    # print(boundaries_df_copy)
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

            # midpoints = boundaries_df_copy.loc[mask, ['x', 'y', 'z']].apply(
            #     lambda row: np.array(row) + 0.01 * (X1 - np.array(row)),
            #     axis=1,
            #     result_type='expand'
            # )

            midpoints = boundaries_df_copy.loc[mask, ['x', 'y', 'z']].apply(
                lambda row: np.array(row) + 0.01 * (X1 - np.array(row)),
                axis=1,
                result_type='expand'
            )
            midpoints.columns = ['x', 'y', 'z']
            boundaries_df_copy.loc[mask, 'x'] = midpoints['x'].values
            boundaries_df_copy.loc[mask, 'y'] = midpoints['y'].values
            boundaries_df_copy.loc[mask, 'z'] = midpoints['z'].values

            boundaries_df_copy.loc[mask, "length"] *= G.nodes[i]['aperture']
            # print('boundaries after logic')
            # print(boundaries_df_copy)
            # print(G.nodes[i]['aperture'])
            # print(boundaries_df_copy.loc[mask, "length"])
            # print(boundaries_df_copy.loc[mask, "length"] * G.nodes[i]['aperture'])
            # print('exit')
            # exit()
            # boundaries_df_copy.rename(columns={'length': 'area'})
            # print(f'mask!!!!!!!!!!!!!!!!!!! {mask}')
            # print(f'L_values!!!!!!!!!!!!!!!!!!!! {L_values}')
            # print(f'midpoints!!!!!!!!!!!!!!!! {midpoints}')
            # print(f'boundaries_df !!!!!!!!!!!!!!!!! {boundaries_df}')
            # exit()


    # # match neg_ids to get the altered_cells id
    # merged = boundaries_df_copy.merge(altered_cells_df_copy[['id', 'neg_id']], on='neg_id')
    # # build and concat to conns_df
    # conns_df_updated = pd.concat([
    #     conns_df_copy,
    #     merged.rename(columns={'id_y': 'i', 'id_x': 'j', 'x': 'xc', 'y': 'yc', 'z': 'zc', 'length': 'area'})[['i', 'j', 'xc', 'yc', 'zc', 'area']]
    # ]).reset_index(drop=True)
    boundary_cells = altered_cells_df_copy.loc[
        altered_cells_df_copy['row_type'] == 'boundary',
        ['id', 'parent_cell_id', 'neg_id', 'x', 'y', 'z', 'orig_length']
    ].rename(columns={'id': 'altered_id'})

    merged = boundaries_df_copy.merge(
        boundary_cells,
        left_on=['neg_id', 'x_orig', 'y_orig', 'z_orig', 'length_orig'],
        right_on=['neg_id', 'x', 'y', 'z', 'orig_length'],
        how='left',
        validate='one_to_one'
    )

    new_conns = merged.rename(columns={
        'altered_id': 'i',
        'id': 'j',
        'x_x': 'xc',
        'y_x': 'yc',
        'z_x': 'zc',
        'length': 'area'
    })[['i', 'j', 'xc', 'yc', 'zc', 'area']]

    conns_df_updated = pd.concat([conns_df_copy, new_conns], ignore_index=True)

    print(f'Total area: {conns_df_updated["area"].sum():.12f}')
    print(f'merged!!!!!!!!!!!!!!!!!!! {merged}')
    print(f'conns_df!!!!!!!!!!!!!!!!!!!! {conns_df_updated}')
    # exit()
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
        xyz = cell_xyz_map[i]
        X1 = np.array([xyz["x"], xyz["y"], xyz["z"]], dtype=float)

        # find all boundary rows for this fracture and calculate L for each
        mask = boundaries_df_copy["id"] == i
        boundary_xyz = boundaries_df_copy.loc[mask, ["x", "y", "z"]].to_numpy(dtype=float)

        L_values = np.linalg.norm(boundary_xyz - X1, axis=1)

        volume_check_og = float(cell_volume_map[i])
        print(volume_check_og)

        n_additions = len(L_values)
        new_fraction_each = 0.01058 #4rects
        # new_fraction_each = 0.017805 #inflow constant_pressure
        # new_fraction_each = 0.1


        total_new_fraction = new_fraction_each * n_additions
        original_fraction = 1.0 - total_new_fraction

        # original cell keeps the remaining percentage
        original_volume = volume_check_og * original_fraction
        cells_df_copy.loc[cells_df_copy["id"] == i, "volume"] = original_volume

        # each new boundary gets 5% of original volume
        boundary_volume = volume_check_og * new_fraction_each
        boundaries_df_copy.loc[mask, "length"] = boundary_volume

        # check: cell + all boundaries should equal original
        new_total = original_volume + boundaries_df_copy.loc[mask, "length"].sum()

        diff = volume_check_og - new_total
        print(f"Volume difference for {i} = {diff:.6f}")

        # use tolerance instead of !=
        if not np.isclose(new_total, volume_check_og, atol=1e-9):
            print(f"ERROR: volume not conserved for cell {i}")
            exit()






        # # check: cell + all boundaries should equal original
        # new_total = new_volume + boundaries_df_copy.loc[mask, "length"].sum()
        # print(f"Volume difference for {i} = {volume_check_og - new_total:.6f}")
        # if new_total != 0
        #     exit()

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
