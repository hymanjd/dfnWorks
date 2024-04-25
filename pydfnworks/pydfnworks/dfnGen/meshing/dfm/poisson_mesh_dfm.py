
import os
import sys
import shutil
import glob 
import numpy as np 

# pydfnworks Modules
from pydfnworks.dfnGen.meshing.mesh_dfn import mesh_dfn_helper as mh
from pydfnworks.dfnGen.meshing.dfm.mesh_dfm_lagrit_scripts import dfm_fracture_facets,  dfm_facets, dfm_diagnostics, check_dfm_mesh
from pydfnworks.dfnGen.meshing.dfm.dfm_helper_functions import dfm_get_box_domain 

def create_user_resolution_function():
    lagrit_script = """

read / avs / full_mesh.inp / MO_INTERSECTIONS

define / MO_H_FIELD / mo_poi_h_field

compute / distance_field / MO_H_FIELD / MO_INTERSECTIONS / dfield

math/multiply/MO_H_FIELD/h_field_att/1,0,0/MO_H_FIELD/dfield/SLOPE/

math/add/MO_H_FIELD/h_field_att/1,0,0/MO_H_FIELD/h_field_att/PARAM_B/

math / floor /   MO_H_FIELD / h_field_att / 1 0 0 / &
                 MO_H_FIELD / h_field_att / H_SCALE

math / ceiling / MO_H_FIELD / h_field_att / 1 0 0 / &
                 MO_H_FIELD / h_field_att / MAX_H_SCALE

# dump / avs / h_field_out.inp / MO_H_FIELD

cmo / delete / MO_INTERSECTIONS
finish


"""
    with open('user_resolution.mlgi','w') as fp:
        fp.write(lagrit_script)
        fp.flush()


def poisson_dfm_driver(num_fracs, h, max_resolution, slope, intercept, box_domain, psets):
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

##################### Poisson 3D ##################### 

define / POI_XMIN / {box_domain['x0']}
define / POI_YMIN / {box_domain['y0']} 
define / POI_ZMIN / {box_domain['z0']}
define / POI_XMAX / {box_domain['x1']}
define / POI_YMAX / {box_domain['y1']}
define / POI_ZMAX / {box_domain['z1']}
define / H_SCALE / {h/2}
define / MAX_H_SCALE / {max_resolution} 
# Slope
define / SLOPE / {slope}
# Intersect
define / PARAM_B / {intercept} 


createpts / poisson_disk / 3d_box / MO_BACKGROUND / H_SCALE / & 
    connect / user_resolution.mlgi 

cmo / status / brief 

define / VERTEX_CLOSE / {h / 4}

# extract_fracture_facets.mlgi must be customized for the number of fractures in the DFN
#
# This is the dfnWorks DFN mesh

define / INPUT / full_mesh.inp
read / avs / INPUT / mo_dfn
cmo / DELATT / mo_dfn / dfield
cmo / DELATT / mo_dfn / b_a
cmo / DELATT / mo_dfn / numbnd
cmo / DELATT / mo_dfn / if_numb

# Diagnostics on fracture mesh extents and resolution

cmo / printatt / mo_dfn / -xyz- / minmax
quality
quality/edge_min
quality/edge_max

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

# cmo / addatt / mo_dfm / p_corner_rfb / vint / scalar / nnodes
# cmo/setatt / mo_dfm / p_corner_rfb / 1,0,0 / 0
# cmo/setatt / mo_dfm / p_corner_rfb /pset,get,p_corner_rfb / 1

# corners of the mesh 3
pset / p_tmp / inter / pleft ptop
pset / p_corner_lft / inter / p_tmp pfront 

pset / pbottom / not / pbottom p_corner_lft
pset / pleft / not / pleft p_corner_lft
pset / pfront / not / pfront p_corner_lft
pset / p_tmp / delete 

# cmo / addatt / mo_dfm / p_corner_lft / vint / scalar / nnodes
# cmo/setatt / mo_dfm / p_corner_lft / 1,0,0 / 0
# cmo/setatt / mo_dfm / p_corner_lft /pset,get,p_corner_lft / 1

# corners of the mesh 4
pset / p_tmp / inter / pright ptop 
pset / p_corner_rft / inter / p_tmp pfront 
pset / p_tmp / delete 

pset / ptop / not / ptop p_corner_rft
pset / pright / not / pright p_corner_rft
pset / pfront / not / pfront p_corner_rft


# cmo / addatt / mo_dfm / p_corner_rft / vint / scalar / nnodes
# cmo/setatt / mo_dfm / p_corner_rft / 1,0,0 / 0
# cmo/setatt / mo_dfm / p_corner_rft /pset,get,p_corner_rft / 1


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


cmo /status / brief 
dump / dfm_tet_w_psets.inp / mo_dfm
dump / exo / dfm_poisson.exo / mo_dfm / psets / / &
     facesets &
"""
        lagrit_script += floop 
        lagrit_script += """
finish
"""
    else: ## no psets
        lagrit_script += """
cmo /status / brief

cmo / delete / mo_dfn 
cmo / delete / mo_df_work
cmo / delete / mo_surf  

cmo /status / 

dump / exo / dfm_poisson.exo / mo_dfm / / / &
     facesets &
"""
        lagrit_script += floop 
        lagrit_script += """
finish
"""

    with open('dfm_poisson3D_mesh.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush()

    print("Creating dfm_poisson3D_mesh.lgi file: Complete\n")

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
        "dfm_poisson3D_mesh.lgi",
        quiet=False)
    
def create_poisson_dfm(self, allowed_percentage, psets, min_dist, max_dist, max_resolution_factor, uniform_mesh):

    slope, intercept = mh.compute_mesh_slope_and_intercept(
        self.h, min_dist, max_dist, max_resolution_factor, uniform_mesh)

    create_user_resolution_function()
    dfm_fracture_facets(self.num_frac)
    dfm_facets(self.h)
    dfm_diagnostics(self.h)
    box_domain = dfm_get_box_domain(self.domain)
    poisson_dfm_driver(self.num_frac, self.h, max_resolution_factor, slope, intercept, box_domain, psets)
    create_dfm() 
    check_dfm_mesh(allowed_percentage)
