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
    print('*' * 80)
    print("--> Dumping h5 file for permeability")
    print('--> Note: This script assumes isotropic permeability')
    print(f'\n--> Opening HDF5 File {filename}')
    # with h5py.File(filename, mode='w') as h5file:
    #     print('--> Allocating cell index array')
    #     print('--> Writing cell indices')
    #     iarray = np.arange(1, self.num_frac + 1) #need to add the number of boundary cells also
    #     #above will turn into this
    #     #iarray = np.arange(1, self.num_cells + 1)
    #     dataset_name = 'Cell Ids'
    #     h5dset = h5file.create_dataset(dataset_name, data=iarray)
    #     print('--> Writing Permeability')
    #     dataset_name = 'Permeability' #this is n fractures long 
    #     h5dset = h5file.create_dataset(dataset_name, data=self.perm)
    cells_df_copy = cells_df.copy()

    # original fracture permeability map
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
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        h5file.create_dataset('Cell Ids', data=iarray)

        print('--> Writing Permeability')
        h5file.create_dataset('Permeability', data=perm_array)

    ## debug print to check the output from h5 file
    # with h5py.File("dfn_properties.h5", "r") as f:
    #     cell_ids = f["Cell Ids"]
    #     perm = f["Permeability"]

    #     print("Cell Ids shape:", cell_ids.shape)
    #     print("Cell Ids size :", cell_ids.size)
    #     print("Cell Ids sample:", cell_ids[:10])

    #     print("Permeability shape:", perm.shape)
    #     print("Permeability size :", perm.size)
    #     print("Permeability sample:", perm[:10])

# def make cells for the perm function above
#h5dump(filename) > perm.out