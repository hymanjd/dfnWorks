"""
.. module:: hex_mesh_dfm.py
   :synopsis: meshing driver for conforming DFM  using a background hex mesh
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>

"""

import os
import sys
import glob 
import shutil 

# pydfnworks Modules
from pydfnworks.dfnGen.meshing.mesh_dfn import mesh_dfn_helper as mh
from pydfnworks.dfnGen.meshing.dfm.mesh_dfm_lagrit_scripts import dfm_fracture_facets,  dfm_facets, dfm_diagnostics

def create_domain(domain, h):
    """ Gather domain information. 

    Parameters
    ----------
        domain : dict
            Domain size dictionary from DFN object 
        h : float 
            Meshing length scale from DFN object 

    Returns
    -------
        num_points : int 
            Number of points on side of the domain 
        box_domain : dict
            dictionary of domain min/max for x,y,z

    Notes
    ------
        Exits program is too many points are in domain. 
        Assuming that 

    """
    box_domain = {"x0": None, "x0": None,
                  "y0": None, "y1": None, 
                  "z0": None, "z1": None 
                  }

    # Extent of domain
    box_domain['x0'] = - 0.5*domain['x']
    box_domain['x1'] = 0.5*domain['x'] 
    box_domain['y0'] = - 0.5*domain['y'] 
    box_domain['y1'] = 0.5*domain['y'] 
    box_domain['z0'] = - 0.5*domain['z'] 
    box_domain['z1'] = 0.5*domain['z'] 

    # Mesh size in matrix
    l = h/2
    # # Number of points in each direction in matrix
    num_points = domain['x'] / l + 1
    if num_points**3 > 1e8:
        error = f"Error: Too many elements for DFM meshing.\nValue {num_points**3}\nMaximum is 1e8\nExiting Program"
        sys.stderr.write(error)
        sys.exit(1)

    num_points_x = domain['x'] / l + 1
    num_points_y = domain['y'] / l + 1
    num_points_z = domain['z'] / l + 1
    if num_points_x*num_points_y*num_points_z > 1e8:
        error = f"Error. Too many elements for DFM meshing.\nValue {num_points_x*num_points_y*num_points_z }\nMaximum is 1e8\nExiting Program"
        sys.stderr.write(error)
        sys.exit(1)
    return box_domain, num_points_x, num_points_y, num_points_z 

def hex_dfm_driver(num_points_x, num_points_y, num_points_z , num_fracs, h):
    """ This function creates the main lagrit driver script, which calls the other lagrit scripts.

    Parameters
    ----------
        num_points : int 
            Number of points on side of the domain 
        num_fracs : int 
            Number of Fractures in the DFN
        h : float
            meshing length scale 

    Returns
    -------
        None

    Notes
    -----
        None 
    """
    floop = ""
    for ifrac in range(1, num_fracs + 1):
        if ifrac < num_fracs:
            floop += f"facets_f{ifrac}.table &\n"
        else:
            floop += f"facets_f{ifrac}.table &\n"
            floop += "left.table &\n"
            floop += "right.table &\n"
            floop += "front.table &\n"
            floop += "back.table &\n"
            floop += "top.table &\n"
            floop += "bottom.table"
            
    lagrit_script  = f"""#
#   dfm_mesh_fracture_driver.lgi
#   dfm_box_dimensions.mlgi
#   dfm_build_background_mesh.mlgi
#   dfm_extract_fracture_facets.mlgi
#   dfm_extract_facets.mlgi
#
# extract_fracture_facets.mlgi must be customized for the number of fractures in the DFN
#
# This is the dfnWorks DFN mesh
#
define / INPUT / full_mesh.inp
read / avs / INPUT / mo_dfn
cmo / DELATT / mo_dfn / dfield
cmo / DELATT / mo_dfn / b_a
cmo / DELATT / mo_dfn / numbnd
cmo / DELATT / mo_dfn / if_numb
#
# Diagnostics on fracture mesh extents and resolution
#
cmo / printatt / mo_dfn / -xyz- / minmax
quality
quality/edge_min
quality/edge_max
#
# Define a resolution for the background mesh. This assumes the DFN
# triangulation is uniform resolution triangles. No attempt is made
# to adapt the volume mesh resolution to the DFN triangle resolution.
#
define / NPX / {num_points_x}
# define / NPXM1 / {num_points_x - 1}
define / NPY / {num_points_y}
# define / NPYM1 / {num_points_y - 1}
define / NPZ / {num_points_z}
# define / NPZM1 / {num_points_z - 1}
define / VERTEX_CLOSE / {h / 4}
#
define / MO_BACKGROUND / mo_background
infile dfm_box_dimensions.mlgi
infile dfm_build_background_mesh.mlgi
#
# Remove all vertices of the tet mesh that fall withing a circumsphere of a fracture triangle.
#
addmesh / excavate / mo_tmp / MO_BACKGROUND / mo_dfn
cmo / delete / MO_BACKGROUND
#
# Merge the vertices of the excavated tet mesh with the DFN vertices
#
cmo / create / mo_dfm / / / tet
copypts / mo_dfm / mo_tmp
cmo / delete / mo_tmp
compute / distance_field / mo_dfm / mo_dfn / df_vertex
cmo / select / mo_dfm
pset / pdel / attribute / df_vertex / 1 0 0 / le VERTEX_CLOSE
rmpoint / pset get pdel / inclusive
rmpoint / compress
cmo / DELATT / mo_dfm / df_vertex
copypts / mo_dfm / mo_dfn
#
cmo / setatt / mo_dfm / imt / 1 0 0 / 1
cmo / setatt / mo_dfm / itp / 1 0 0 / 0
cmo / select / mo_dfm
connect
cmo / setatt / mo_dfm / itetclr / 1 0 0 / 1
resetpts / itp
quality
#
#compute / signed_distance_field / mo_dfm / mo_dfn / df_sign_dfm_dfn
#
# crush_thin_tets / mo_dfm / 0.25 / 1 0 0 
dump / avs    / dfm_tet_mesh.inp / mo_dfm
dump / lagrit / dfm_tet_mesh.lg  / mo_dfm
dump / exo    / dfm_tet_mesh.exo / mo_dfm

cmo / delete / mo_dfm
cmo / delete / mo_dfn
#
cmo / status / brief
#
infile dfm_extract_fracture_facets.mlgi
infile dfm_diagnostics.mlgi
#
# Delete this !!!! 
# Hardcoded facesets on boundaries for Alex EES17
cmo / select / mo_dfm
extract / surfmesh / 1 0 0 / mo_surf / mo_dfm / external
cmo / addatt / mo_surf / id_side / vint / scalar / nelements
cmo / select / mo_surf
settets / normal
cmo / copyatt / mo_surf mo_surf / id_side itetclr
cmo / printatt / mo_surf / id_side / minmax
cmo / DELATT / mo_surf / itetclr0
cmo / DELATT / mo_surf / idnode0
cmo / DELATT / mo_surf / idelem0
cmo / DELATT / mo_surf / facecol
cmo / DELATT / mo_surf / itetclr1
cmo / DELATT / mo_surf / idface0
#
cmo / copy / mo_tmp / mo_surf
cmo / select / mo_tmp
eltset / e_bottom / id_side / eq / 1
eltset / e_delete / not / e_bottom
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / DELATT / mo_tmp / id_side
dump / avs2 / bottom.table / mo_tmp / 0 0 0 2
cmo / delete / mo_tmp
#
cmo / copy / mo_tmp / mo_surf
cmo / select / mo_tmp
eltset / e_top / id_side / eq / 2
eltset / e_delete / not / e_top
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / DELATT / mo_tmp / id_side
dump / avs2 / top.table / mo_tmp / 0 0 0 2
cmo / delete / mo_tmp
#
cmo / copy / mo_tmp / mo_surf
cmo / select / mo_tmp
eltset / e_right / id_side / eq / 3
eltset / e_delete / not / e_right
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / DELATT / mo_tmp / id_side
dump / avs2 / right.table / mo_tmp / 0 0 0 2
cmo / delete / mo_tmp
#
cmo / copy / mo_tmp / mo_surf
cmo / select / mo_tmp
eltset / e_back / id_side / eq / 4
eltset / e_delete / not / e_back
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / DELATT / mo_tmp / id_side
dump / avs2 / back.table / mo_tmp / 0 0 0 2
cmo / delete / mo_tmp
#
cmo / copy / mo_tmp / mo_surf
cmo / select / mo_tmp
eltset / e_left / id_side / eq / 5
eltset / e_delete / not / e_left
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / DELATT / mo_tmp / id_side
dump / avs2 / left.table / mo_tmp / 0 0 0 2
cmo / delete / mo_tmp
#
cmo / copy / mo_tmp / mo_surf
cmo / select / mo_tmp
eltset / e_front / id_side / eq / 6
eltset / e_delete / not / e_front
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / DELATT / mo_tmp / id_side
dump / avs2 / front.table / mo_tmp / 0 0 0 2
cmo / delete / mo_tmp
#
dump / exo / dfm_tet_mesh_w_fsets.exo / mo_dfm / / / &
     facesets &
"""
    lagrit_script += floop 
    lagrit_script += """
finish
"""

    with open('dfm_mesh_fracture_driver.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush()

    print("Creating dfm_mesh_fracture_driver.lgi file: Complete\n")

def dfm_box(box_domain):    
    """ This function creates the dfm_box_dimensions.mlgi lagrit script.

    Parameters
    ----------
        box_domain : dict
            dictionary of domain min/max for x,y,z
  
    Returns
    -------
        None 

    Notes
    -----
        None 

    """

    lagrit_script = f"""#
# Define a bounding box that surrounds, and is a big bigger, than the DFN
#
define / X0 / {box_domain['x0']}
define / X1 / {box_domain['x1']}
define / Y0 / {box_domain['y0']}
define / Y1 / {box_domain['y1']}
define / Z0 / {box_domain['z0']}
define / Z1 / {box_domain['z1']}

finish
"""
    with open('dfm_box_dimensions.mlgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush()

    print("Creating dfm_box_dimensions.mlgi file: Complete\n")

def dfm_build():
    """ Create the dfm_build_background_mesh.mlgi lagrit script.

    Parameters
    ----------
        None 

    Returns
    -------
        None 

    Notes
    -----
        Needs to be modified to have different NPX, NPY, NPZ 
    """

    lagrit_script = """#
# Build a uniform background point distribution.
#
cmo / create / MO_BACKGROUND / / / tet
createpts / xyz / NPX NPY NPZ / X0 Y0 Z0 / X1 Y1 Z1 / 1 1 1
cmo / setatt / MO_BACKGROUND / imt / 1 0 0 / 1
connect / noadd
cmo / setatt / MO_BACKGROUND / itetclr / 1 0 0 / 1
#
finish
"""
    with open('dfm_build_background_mesh.mlgi', 'w') as fp: 
        fp.write(lagrit_script)
        fp.flush()
    print("Creating dfm_box_dimensions.mlgi file: Complete\n")

def create_dfm():
    """ This function executes the lagrit scripts. 
    
    Parameters
    ----------
        None    

    Returns
    -------
        None
    
    Notes
    -----
        None
 
    """
    # Run LaGriT
    mh.run_lagrit_script(
        "dfm_mesh_fracture_driver.lgi",
        quiet=False)

def cleanup_mesh_dfm_directory():
    """ Clean up working files from meshing the DFM

    Parameters
    ---------------
        None

    Returns
    ----------------
        None

    Notes
    ---------------
        None

    """
    print("--> Cleaning up working directory")
    # clean up LaGrit Scripts
    lagrit_script_dir = "dfm_lagrit_files" 
    try:
        os.mkdir(lagrit_script_dir)
    except:
        shutil.rmtree(lagrit_script_dir)
        os.mkdir(lagrit_script_dir)
    lagrit_scripts = glob.glob("*lgi")
    for filename in lagrit_scripts:
        shutil.copyfile(filename, lagrit_script_dir + os.sep + filename)
        os.remove(filename)

    extra_files = ['dfm_mesh_fracture_driver.lgi.log','dfm_mesh_fracture_driver.lgi.out',
                   'tmp_interpolate.inp']
    for filename in extra_files:
        shutil.copyfile(filename, lagrit_script_dir + os.sep + filename)
        os.remove(filename)

    table_dir = "tables"
    try:
        os.mkdir(table_dir)
    except:
        shutil.rmtree(table_dir)
        os.mkdir(table_dir)

    table_files = glob.glob("*table")
    for filename in table_files:
        shutil.copyfile(filename, table_dir + os.sep + filename)
        os.remove(filename)

    facets_dir = "facets"
    try:
        os.mkdir(facets_dir)
    except:
        shutil.rmtree(facets_dir)
        os.mkdir(facets_dir)

    facet_files = glob.glob("facets*inp")
    for filename in facet_files:
        shutil.copyfile(filename, facets_dir + os.sep + filename)
        os.remove(filename)


    print("--> Cleaning up working directory: Complete")


def check_dfm_mesh(allowed_percentage):
    """ Checks how many elements of the DFN meshing are missinf from the DFM. If the percentage missing is larger than the allowed percentage, then the program exists.

    Parameters
    ----------------
        allowed_percentage : float
            Percentage of the mesh allowed to be missing and still continue

    Returns
    ----------
        None

    Notes
    ----------
        None
    
    """

    print("--> Checking for missing elements")
    if os.path.isfile('missed_cells_full_mesh.inp'):
        print("--> Missing elements have been found.")
        print(f"--> Missing elements are in the file 'missed_cells_full_mesh.inp' if you want to see them.")
        # get number of missed elements in the 
        with open('missed_cells_full_mesh.inp', 'r') as fp:
            line = fp.readline().split()
            missing_num_elems = int(line[1])
        # get the total number of elements

        with open('full_mesh.inp', 'r') as fp:
            line = fp.readline().split()
            total_num_elems = int(line[1])
        # Compute percentage and compare
        missing_percent = 100*(missing_num_elems/total_num_elems)
        print(f"--> Out of {total_num_elems} elements in the DFN there are {missing_num_elems} missing from the DFM.")
        print(f"--> That's {missing_percent:0.2f} percent of the mesh.")

        if  missing_percent > allowed_percentage:
            error = f"*** Error. Missing percent of mesh is larger than tolerance {allowed_percentage} ***\n*** Exitting ***\n "
            sys.stderr.write(error)
            sys.exit(1)
        else:
            print("--> Doesn't seem to bad. Keep Calm and Carry on.")

    # if the file 'missed_cells_full_mesh.inp' does not exists, this means no elements were missed.  
    else:
        print("--> No missinng elements found. ")

def create_hex_dfm(self, allowed_percentage, cleanup):
    box_domain, num_points_x, num_points_y, num_points_z  = create_domain(self.domain, self.h)
    hex_dfm_driver(num_points_x, num_points_y, num_points_z , self.num_frac, self.h)
    dfm_box(box_domain)    
    dfm_build()
    dfm_fracture_facets(self.num_frac)
    dfm_facets()
    dfm_diagnostics(self.h)
    create_dfm()
    check_dfm_mesh(allowed_percentage)
    if cleanup:
        cleanup_mesh_dfm_directory()
