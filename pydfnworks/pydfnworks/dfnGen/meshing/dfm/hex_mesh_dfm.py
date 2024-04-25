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
from pydfnworks.dfnGen.meshing.dfm.mesh_dfm_lagrit_scripts import dfm_fracture_facets,  dfm_facets, dfm_diagnostics, check_dfm_mesh
from pydfnworks.dfnGen.meshing.dfm.dfm_helper_functions import dfm_get_box_domain 

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
    box_domain =  dfm_get_box_domain(domain)

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

def hex_dfm_driver(box_domain, num_points_x, num_points_y, num_points_z , num_fracs, h, psets):
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
"""
    
    if psets:
        eps = h/4
        lagrit_script += f"""
cmo / select / mo_dfm
cmo / printatt / mo_dfm / -xyz- / minmax

pset/ pleft / geom / xyz / 1, 0, 0 /  &
     {box_domain['x0'] - eps} {box_domain['y0']} {box_domain['z0']} / {box_domain['x0'] + eps} {box_domain['y1']} {box_domain['z1']}  / 0,0,0
pset/ pright / geom / xyz / 1, 0, 0 / &
    {box_domain['x1'] - eps} {box_domain['y0']} {box_domain['z0']} / {box_domain['x1'] + eps} {box_domain['y1']} {box_domain['z1']}  / 0,0,0

pset / pfront / geom / xyz / 1, 0, 0 / & 
    {box_domain['x0']} {box_domain['y0'] - eps}  {box_domain['z0']} / {box_domain['x1']}  {box_domain['y0'] + eps}  {box_domain['z1']}  / 0,0,0 
pset / pback / geom / xyz / 1, 0, 0 / & 
    {box_domain['x0']} {box_domain['y1'] - eps}  {box_domain['z0']}  / {box_domain['x1']}  {box_domain['y1'] + eps}  {box_domain['z1']}  / 0,0,0 

pset / pbottom / geom / xyz / 1, 0, 0 / &
    {box_domain['x0']} {box_domain['y0']} {box_domain['z0'] - eps} / {box_domain['x1']}  {box_domain['y1']} {box_domain['z0'] + eps}/ 0,0,0 
pset / ptop / geom / xyz / 1, 0, 0 /  & 
    {box_domain['x0']} {box_domain['y0']} {box_domain['z1'] - eps} / {box_domain['x1']}  {box_domain['y1']} {box_domain['z1'] + eps} / 0,0,0 

# corners of the mesh 1
pset / p_tmp / inter / pleft pbottom
pset / p_corner_lfb / inter / p_tmp pfront 
pset / p_tmp / delete 

pset / pbottom / not / pbottom p_corner_lfb
pset / pleft / not / pleft p_corner_lfb
pset / pfront / not / pfront p_corner_lfb


cmo / addatt / mo_dfm / p_corner_lfb / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_lfb / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_lfb /pset,get,p_corner_lfb / 1

# corners of the mesh 2
pset / p_tmp / inter / pright pbottom
pset / p_corner_rfb / inter / p_tmp pfront 
pset / p_tmp / delete 

pset / pbottom / not / pbottom p_corner_rfb
pset / pright / not / pright p_rfp_corner
pset / pfront / not / pfront p_corner_rfb

cmo / addatt / mo_dfm / p_corner_rfb / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_rfb / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_rfb /pset,get,p_corner_rfb / 1

# corners of the mesh 3
pset / p_tmp / inter / pleft ptop
pset / p_corner_lft / inter / p_tmp pfront 

pset / pbottom / not / pbottom p_corner_lft
pset / pleft / not / pleft p_corner_lft
pset / pfront / not / pfront p_corner_lft
pset / p_tmp / delete 

cmo / addatt / mo_dfm / p_corner_lft / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_lft / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_lft /pset,get,p_corner_lft / 1

# corners of the mesh 4
pset / p_tmp / inter / pright ptop 
pset / p_corner_rft / inter / p_tmp pfront 
pset / p_tmp / delete 

pset / ptop / not / ptop p_corner_rft
pset / pright / not / pright p_corner_rft
pset / pfront / not / pfront p_corner_rft



cmo / addatt / mo_dfm / p_corner_rft / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_rft / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_rft /pset,get,p_corner_rft / 1


### back face 
# corners of the mesh 1
pset / p_tmp / inter / pleft pbottom
pset / p_corner_lbb / inter / p_tmp pback 
pset / p_tmp / delete 


# corners of the mesh 2
pset / p_tmp / inter / pright pbottom
pset / p_corner_rbb / inter / p_tmp pback 
pset / p_tmp / delete 


# corners of the mesh 3
pset / p_tmp / inter / pleft ptop
pset / p_corner_lbt / inter / p_tmp pback 
pset / p_tmp / delete 


# corners of the mesh 4
pset / p_tmp / inter / pright ptop 
pset / p_corner_rbt / inter / p_tmp pback 
pset / p_tmp / delete 

########

cmo / addatt / mo_dfm / p_corner_rbt / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_rbt / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_rbt /pset,get,p_corner_rbt / 1

cmo / addatt / mo_dfm / p_corner_lbt / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_lbt / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_lbt /pset,get,p_corner_lbt / 1


cmo / addatt / mo_dfm / p_corner_lbb / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_lbb / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_lbb /pset,get,p_corner_lbb / 1

cmo / addatt / mo_dfm / p_corner_rbb / vint / scalar / nnodes
cmo/setatt / mo_dfm / p_corner_rbb / 1,0,0 / 0
cmo/setatt / mo_dfm / p_corner_rbb /pset,get,p_corner_rbb / 1

## clean up PSETS TO MESH 
pset / pbottom / not / pbottom p_corner_lbb
pset / pleft / not / pleft p_corner_lbb
pset / pback / not / pback p_corner_lbb

pset / pbottom / not / pbottom p_corner_rbb
pset / pright / not / pright p_corner_rbb
pset / pback / not / pback p_corner_rbb

pset / ptop / not / ptop p_corner_lbt
pset / pleft / not / pleft p_corner_lbt
pset / pback / not / pback p_corner_lbt

pset / ptop / not / ptop p_corner_rbt
pset / pright / not / pright p_corner_rbt
pset / pback / not / pback p_corner_rbt


pset / pbottom / not / pbottom p_corner_lfb
pset / pleft / not / pleft p_corner_lfb
pset / pfront / not / pfront p_corner_lfb 

pset / pbottom / not / pbottom p_corner_rfb
pset / pright / not / pright p_corner_rfb
pset / pfront / not / pfront p_corner_rfb

pset / ptop / not / ptop p_corner_lft
pset / pleft / not / pleft p_corner_lft
pset / pfront / not / pfront p_corner_lft

pset / ptop / not / ptop p_corner_rft
pset / pright / not / pright p_corner_rft
pset / pfront / not / pfront p_corner_rft


### edges ##### 

pset / p_edge_lb / inter / pleft pbottom
pset / pbottom / not / pbottom p_edge_lb
pset / pleft / not / pleft p_edge_lb

pset / p_edge_lt / inter / pleft ptop
pset / ptop / not / ptop p_edge_lt
pset / pleft / not / pleft p_edge_lt

pset / p_edge_rb / inter / pright pbottom
pset / pbottom / not / pbottom p_edge_rb
pset / pright / not / pright p_edge_rb

pset / p_edge_rt / inter / pright ptop 
pset / ptop / not / ptop p_edge_rt
pset / pright / not / pright p_edge_rt

####### 
pset / p_edge_lfr / inter / pleft pfront
pset / pleft / not / pleft p_edge_lfr
pset / pfront / not / pfront p_edge_lfr

pset / p_edge_lba / inter / pleft pback 
pset / pleft / not / pleft p_edge_lba
pset / pback / not / pback p_edge_lba

pset / p_edge_rfr / inter / pright pfront
pset / pright / not / pright p_edge_rfr
pset / pfront / not / pfront p_edge_rfr

pset / p_edge_rba / inter / pright pback 
pset / pright / not / pright p_edge_rba
pset / pback / not / pback p_edge_rba

####### 


pset / p_edge_frb / inter / pfront pbottom
pset / pfront / not / pfront p_edge_frb
pset / pbottom / not / pbottom p_edge_frb

pset / p_edge_bab / inter / pback pbottom
pset / pback / not / pback p_edge_bab
pset / pbottom / not / pbottom p_edge_bab

pset / p_edge_frtop / inter / pfront ptop
pset / pfront / not / pfront p_edge_frtop
pset / ptop / not / ptop p_edge_frtop

pset / p_edge_btop / inter /  pback ptop
pset / pback / not / pback p_edge_btop
pset / ptop / not / ptop p_edge_btop

####### 

cmo / addatt / mo_dfm / right / vint / scalar / nnodes
cmo/setatt / mo_dfm / right / 1,0,0 / 0
cmo/setatt / mo_dfm / right /pset,get,pright / 1

cmo / addatt / mo_dfm / back / vint / scalar / nnodes
cmo/setatt / mo_dfm / back / 1,0,0 / 0
cmo/setatt / mo_dfm / back /pset,get,pback / 1


cmo / addatt / mo_dfm / left / vint / scalar / nnodes
cmo/setatt / mo_dfm / left / 1,0,0 / 0
cmo/setatt / mo_dfm / left /pset,get,pleft / 1

cmo / addatt / mo_dfm / top / vint / scalar / nnodes
cmo/setatt / mo_dfm / top / 1,0,0 / 0
cmo/setatt / mo_dfm / top /pset,get,ptop / 1

cmo / addatt / mo_dfm / bottom / vint / scalar / nnodes
cmo/setatt / mo_dfm / bottom / 1,0,0 / 0
cmo/setatt / mo_dfm / bottom /pset,get,pbottom / 1

cmo / addatt / mo_dfm / front / vint / scalar / nnodes
cmo/setatt / mo_dfm / front / 1,0,0 / 0
cmo/setatt / mo_dfm / front /pset,get,pfront / 1

dump / dfm_tet_w_psets.inp / mo_dfm
dump / exo / dfm_tet_mesh_w_p_f_sets.exo / mo_dfm / psets / / &
     facesets &
"""
        lagrit_script += floop 
        lagrit_script += """
finish
"""
    else: ## no psets
        lagrit_script += """
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


def create_hex_dfm(self, allowed_percentage, psets):

    box_domain, num_points_x, num_points_y, num_points_z  = create_domain(self.domain, self.h)
    hex_dfm_driver(box_domain, num_points_x, num_points_y, num_points_z , self.num_frac, self.h, psets)
    dfm_box(box_domain)    
    dfm_build()
    dfm_fracture_facets(self.num_frac)
    dfm_facets()
    dfm_diagnostics(self.h)
    create_dfm()
    check_dfm_mesh(allowed_percentage)
