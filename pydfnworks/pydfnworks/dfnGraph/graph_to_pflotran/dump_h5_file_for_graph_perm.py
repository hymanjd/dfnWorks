import h5py
import numpy as np
import pandas as pd

def dump_h5_file_for_graph_perm(self,cells_df, filename):
    """ Write permeability values to cell ids and permeability values to dfn_properties.h5 file for pflotran. 

    Parameters
    ----------
        self : object
            DFN Class
        filename : string
            name of f5 file, default is permeability

    Returns
    ---------
        None

    Notes
    ----------
        Hydraulic properties need to attached to the class prior to running this function. Use DFN.assign_hydraulic_properties() to do so. 
    """
    cells_df_copy = cells_df.copy()

    if self.num_frac == 1 and len(self.perm) != 1:
        perm_value = float(np.mean(self.perm))
        perm_map = pd.Series([perm_value], index=[1])
    else:
        perm_map = pd.Series(self.perm, index=np.arange(1, self.num_frac + 1))

    # choose source id for permeability
    source_ids = np.where(
        cells_df_copy["row_type"] == "boundary",
        cells_df_copy["parent_cell_id"],
        cells_df_copy["id"]
    )

    source_ids = pd.Series(source_ids).astype(int)
    perm_array = perm_map.loc[source_ids].to_numpy(dtype=float)

    iarray = cells_df_copy["id"].to_numpy(dtype=int)

    with h5py.File(filename, mode='w') as h5file:
        h5file.create_dataset('Cell Ids', data=iarray)

        h5file.create_dataset('Permeability', data=perm_array)

    ## debug print to check the output from h5 file
    with h5py.File("dfn_properties.h5", "r") as f:
        cell_ids = f["Cell Ids"]
        perm = f["Permeability"]
