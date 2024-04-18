"""
.. module:: mesh_dfm.py
   :synopsis: meshing driver for conforming DFM  
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>

"""

import os
import sys
import shutil
import glob 

# pydfnworks Modules
from pydfnworks.dfnGen.meshing.mesh_dfn import mesh_dfn_helper as mh
from pydfnworks.dfnGen.meshing.dfm.hex_mesh_dfm import hex_mesh_driver 


def setup_mesh_dfm_directory(jobname, dirname):
    """ Setup working directory for meshing the DFM. 

    Parameters
    ----------------
        jobname : string
            path to DFN working directory 
        dirname : string 
            name of working directory

    Returns
    --------------
        None

    Notes
    -------------
        None 
    """
    path = jobname + os.sep + dirname
    try: 
        os.mkdir(path)
        os.chdir(path)
    except:
        shutil.rmtree(path)
        os.mkdir(path)
        os.chdir(path)


    print(f"--> Working directory is now {os.getcwd()}")
    # Make symbolic links to required files
    try:
        os.symlink(jobname + os.sep + "full_mesh.inp", "full_mesh.inp")
    except:
        error = f"Error. Unable to make symbolic link to full_mesh.inp file for DFM meshing from {jobname}.\nExitting program."
        sys.stderr.write(error)
        sys.exit(1)

    print("--> Setting up DFM meshing directory complete")


def mesh_dfm(self, mesh_type, dirname = "dfm_mesh", allowed_percentage = 1, cleanup = True, l = None):
    """" Creates a conforming mesh of a DFN using a uniform background tetrahedron mesh. The DFN must be meshed using a uniform triangular mesh. (DFN.mesh_network(uniform_mesh = True))

    Parameters
    ------------------
        mesh_type : string 
            Background mesh type. Options are hex, poisson, octree
        dirname : string
            name of working directory. Default : dfm_mesh
        allowed_percentage : float
            Percentage of the mesh allowed to be missing and still continue
        cleanup : bool
            Clean up working directory. If true dep files are moved into subdirectories

    Returns
    ---------------
        None

    Notes
    --------------
        The final mesh is output in exodus format. This requires that LaGriT is built against exodus.  
         
    """

    print('=' * 80)
    print("Creating conforming DFM mesh using LaGriT : Starting")
    print('=' * 80)

    setup_mesh_dfm_directory(self.jobname, dirname)

    if mesh_type == "hex":
        self.create_hex_dfm(allowed_percentage, cleanup)
    elif mesh_type == "octree":
        self.create_octree_dfm(l)
    elif mesh_type == "poisson":
        self.create_poisson_dfm()


    print('=' * 80)
    print("Creating conforming DFM mesh using LaGriT : Complete")
    print('=' * 80)

