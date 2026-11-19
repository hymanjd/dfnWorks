import h5py
import numpy as np
import pandas as pd


def dump_h5_file_for_graph_perm(self, cells_df, filename):
    """
    Write cell ids and permeabilities to an HDF5 file for PFLOTRAN.

    Each fracture cell gets its fracture's permeability. Each boundary
    helper cell gets the permeability of its parent fracture.

    Parameters
    ----------
    self : object
        DFN Class. ``perm`` and ``num_frac`` must be set.
    cells_df : pandas.DataFrame
        Cell table from :func:`convert_graph_to_data_frames`, with columns
        ``id``, ``row_type`` and ``parent_cell_id``.
    filename : str
        Name of the HDF5 file to write, e.g. ``"dfn_properties.h5"``.

    Returns
    -------
    None

    Notes
    -----
    Hydraulic properties must be set on the DFN object first, with
    ``DFN.assign_hydraulic_properties()``. Permeability is assumed to be
    isotropic. The file has two datasets, ``Cell Ids`` and ``Permeability``.

    Logs an error and exits if ``self.perm`` does not have exactly one value
    per fracture.
    """
    cells_df_copy = cells_df.copy()

    perm = np.asarray(self.perm, dtype=float)
    if perm.ndim != 1 or perm.size != self.num_frac:
        self.print_log(
            f"Error. Expected one permeability value per fracture ({self.num_frac}), "
            f"but self.perm has shape {perm.shape}. Only isotropic permeability is "
            "supported.", 'error')
    perm_map = pd.Series(perm, index=np.arange(1, self.num_frac + 1))

    # boundary helper cells use their parent fracture's permeability
    source_ids = np.where(cells_df_copy["row_type"] == "boundary",
                          cells_df_copy["parent_cell_id"], cells_df_copy["id"])

    source_ids = pd.Series(source_ids).astype(int)
    perm_array = perm_map.loc[source_ids].to_numpy(dtype=float)

    iarray = cells_df_copy["id"].to_numpy(dtype=int)

    with h5py.File(filename, mode='w') as h5file:
        h5file.create_dataset('Cell Ids', data=iarray)
        h5file.create_dataset('Permeability', data=perm_array)
