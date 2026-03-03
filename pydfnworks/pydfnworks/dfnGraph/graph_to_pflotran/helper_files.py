
import numpy as np
import glob 
from scipy.stats import hmean
import h5py

def process_vtk_files(filename):
    block = False
    data = {}
    files = glob.glob(filename)
    files.sort()
    for filename in files:
        time = filename.split("-")[1][:-4]
        time = filename.split(".")[0][-3:]
        print(filename,time)
        with open(filename,"r") as fp:
            for line in fp.readlines():
                # skip over the mesh
                if "POINT_DATA" in line or "CELL_DATA" in line:
                    # Get number of cells
                    num_cells = int(line.split()[-1])
                    #print("Number of cells {0}".format(num_cells))
                
                #for line in fp.readlines():
                #print(line)
                if block:
                    line = line.split()
                    if not "LOOKUP_TABLE" in line:
                        for val in line:
                            data[variable_name][index] = float(val)
                            index += 1
                            if index >= num_cells:
                                block = False
                
                if 'SCALARS' in line and not block:
                    line = line.split()
                    variable_name = line[1]
                    #print("Processing Variable: {0}".format(variable_name))
                    data[variable_name] = np.zeros(num_cells)
                    index = 0
                    block = True
    return data 


def dump_graph_uge(DFN, G):

    DFN.volume = DFN.surface_area * DFN.aperture  
    for i in range(1, DFN.num_frac + 1):
        G.nodes[i]['volume'] = DFN.volume[i-1]
        G.nodes[i]['center'] = DFN.centers[i-1]

    ## clean this up....
    internal_intersections = []
    for u in G.intersections:
        print(G.nodes[u])
        frac_list = list(G.neighbors(u))
        print(frac_list)
        if len(frac_list) == 2:
            frac1 = frac_list[0]
            frac2 = frac_list[1]
            if frac1 != 's' and frac2 != 's' and frac1 != 't' and frac2 != 't': 
                b1 = G.nodes[frac1]['aperture']
                b2 = G.nodes[frac2]['aperture']
                b_hmean = hmean([b1, b2])
                G.nodes[u]['area'] = G.nodes[u]['length'] * b_hmean
                internal_intersections.append(u)

    with open('graph.uge', 'w') as fp:
        fp.write(f'CELLS\t{DFN.num_frac}\n')
        for i in range(1, DFN.num_frac + 1):
            fp.write(f"{i}\t{G.nodes[i]['center'][0]:.16E}\t{G.nodes[i]['center'][1]:.16E}\t{G.nodes[i]['center'][2]:.16E}\t{G.nodes[i]['volume']:.16E}\n")

        fp.write(f'CONNECTIONS\t{len(internal_intersections)}\n')
        for u in G.intersections:
            print(G.nodes[u])
            frac_list = list(G.neighbors(u))
            print(frac_list)
            if len(frac_list) == 2:
                frac1 = frac_list[0]
                frac2 = frac_list[1]
                if frac1 != 's' and frac2 != 's' and frac1 != 't' and frac2 != 't': 
                    b1 = G.nodes[frac1]['aperture']
                    b2 = G.nodes[frac2]['aperture']
                    b_hmean = hmean([b1, b2])
                    G.nodes[u]['area'] = G.nodes[u]['length'] * b_hmean
                    fp.write(f"{frac1}\t{frac2}\t{G.nodes[u]['x']:.16E}\t{G.nodes[u]['y']:.16E}\t{G.nodes[u]['z']:.16E}\t{G.nodes[u]['area']:.16E}\n")

    print("inflow neighbors")
    inflow_edges= list(G.neighbors('s'))
    inflow_nodes = []
    for u in inflow_edges:
        tmp = list(G.neighbors(u))
        inflow_nodes.extend(tmp)

    inflow_nodes = list(set(inflow_nodes))
    inflow_nodes.remove('s')

    eps = 1e-5
    ## Write boundary ex files 
    with open('inflow_nodes.ex','w') as fp:
        fp.write(f'CONNECTIONS\t{len(inflow_nodes)}\n')
        for u in inflow_nodes:
            fp.write(f"{u}\t{G.nodes[u]['center'][0]-eps:16E}\t{G.nodes[u]['center'][1]-eps:16E}\t{G.nodes[u]['center'][2]-eps:16E}\t{1/eps}\n")

    print("outflow neighbors")
    outflow_edges= list(G.neighbors('t'))
    outflow_nodes = []
    for u in outflow_edges:
        tmp = list(G.neighbors(u))
        outflow_nodes.extend(tmp)

    outflow_nodes = list(set(outflow_nodes))
    outflow_nodes.remove('t')
    ## Write boundary ex files 
    with open('outflow_nodes.ex','w') as fp:
        fp.write(f'CONNECTIONS\t{len(outflow_nodes)}\n')
        for u in outflow_nodes:
            fp.write(f"{u}\t{G.nodes[u]['center'][0]-eps:8E}\t{G.nodes[u]['center'][1]-eps:8E}\t{G.nodes[u]['center'][2]-eps:8E}\t{1/eps}\n")


def dump_h5_file_for_perm(DFN, filename = 'permeability.h5'):
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
        iarray = np.arange(1, DFN.num_frac + 1)
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print('--> Writing Permeability')
        dataset_name = 'Permeability'
        h5dset = h5file.create_dataset(dataset_name, data=DFN.perm)
    print("--> Done writing h5 file")
    print('*' * 80)
    print()
