import h5py
import numpy as np

def dump_h5_file_for_graph_perm(self, filename):
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
    with h5py.File(filename, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        iarray = np.arange(1, self.num_frac + 1)
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print('--> Writing Permeability')
        dataset_name = 'Permeability'
        h5dset = h5file.create_dataset(dataset_name, data=self.perm)