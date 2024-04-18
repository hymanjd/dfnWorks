import os 
import sys
def dfm_fracture_facets(num_frac):
    """ This function creates the dfm_extract_fracture_facets.mlgi lagrit script.

    Parameters
    ----------
        num_frac : int 
            Number of fractures in the DFN
    
    Returns
    -------
        None

    Notes
    -----
        None 
    """
    floop1 = ""
    floop2 = ""
    for ifrac in range(1,num_frac+1):
        floop1 += f"""
define / FRAC_ID / {ifrac}
define / FRAC_FILE_OUT / facets_f{ifrac}.inp
define / FRAC_TABLE_OUT / facets_f{ifrac}.table
#
infile dfm_extract_facets.mlgi
        """
        if ifrac == 1:
            floop2 += f"""
read / avs / facets_f{ifrac}.inp / mo_merge
cmo / setatt / mo_merge / itetclr / 1 0 0 / {ifrac}
        """
        else:
            floop2 += f"""
read / avs / facets_f{ifrac}.inp / mo
cmo / setatt / mo / itetclr / 1 0 0 / {ifrac}
addmesh / merge / mo_merge / mo_merge / mo
cmo / delete / mo
        """
    lagrit_script = """#
define / INPUT / full_mesh.inp
define / MO_ONE_FRAC / mo_tmp_one_fracture
#
read / avs / dfm_tet_mesh.inp / mo_dfm
#
cmo / create / mo_merge
cmo / status / brief
read / avs / INPUT / mo_dfn
cmo / status / brief
""" + floop1 + floop2 + """
dump / avs / facets_merged.inp / mo_merge
cmo / addatt / mo_merge / id_frac / vint / scalar / nelements
cmo / copyatt / mo_merge / mo_merge / id_frac / itetclr
dump / avs / facets_merged.table / mo_merge / 0 0 0 2
cmo / delete / mo_merge

finish
"""
    with open('dfm_extract_fracture_facets.mlgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush()
    print("Creating dfm_extract_fracture_facets.mlgi file: Complete\n")

def dfm_facets():
    """ This function creates the dfm_extract_facets.mlgi lagrit script.

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

    lagrit_script = f"""#
cmo / copy / MO_ONE_FRAC / mo_dfn
cmo / select / MO_ONE_FRAC
rmmat / FRAC_ID / element / exclusive
rmpoint / compress
resetpts / itp
cmo / status / brief
#
compute / signed_distance_field / mo_dfm / MO_ONE_FRAC / dfield_sign

cmo / delete / MO_ONE_FRAC
#
cmo / copy / mo_df_work / mo_dfm

cmo / DELATT / mo_dfm / dfield_sign

cmo / select / mo_df_work
pset / p1 / attribute / dfield_sign / 1 0 0 / gt 0.0
pset / p2 / attribute / dfield_sign / 1 0 0 / lt 0.0
eltset / e1 / inclusive / pset get p1
eltset / e2 / inclusive / pset get p2
cmo / setatt / mo_df_work / itetclr / eltset get e1 / 1
cmo / setatt / mo_df_work / itetclr / eltset get e2 / 2
resetpts / itp
extract / surfmesh / 1 0 0 / MO_ONE_FRAC_EXTRACT / mo_df_work
#
cmo / select / MO_ONE_FRAC_EXTRACT
eltset / edel / idelem0 / eq / 0
rmpoint / element / eltset get edel
rmpoint / compress
pset / pdel / attribute / dfield_sign / 1 0 0 / gt / 1.e-9
rmpoint / pset get pdel / inclusive
rmpoint / compress
#
# idelem0, idelem1 are element numbers
# idface0, idface1 are the face numbers
#
cmo / DELATT / MO_ONE_FRAC_EXTRACT / itetclr0
cmo / DELATT / MO_ONE_FRAC_EXTRACT / itetclr1
cmo / DELATT / MO_ONE_FRAC_EXTRACT / facecol
#
# Don't keep both sides of the fracture face information.
#
cmo / DELATT / MO_ONE_FRAC_EXTRACT / idelem0
cmo / DELATT / MO_ONE_FRAC_EXTRACT / idface0
#
dump / avs2 / FRAC_FILE_OUT  / MO_ONE_FRAC_EXTRACT
dump / avs2 / FRAC_TABLE_OUT / MO_ONE_FRAC_EXTRACT  / 0 0 0 2
#
cmo / delete / MO_ONE_FRAC_EXTRACT
#
cmo / status / brief
#
finish
"""
    with open('dfm_extract_facets.mlgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush()

    print("Creating dfm_extract_facets.mlgi file: Complete\n")


def dfm_diagnostics(h):
    """
    
    """
    eps_offset = 0.1*h
    lagrit_script = f"""

# Figure out which cells (tringles) from DFN full_mesh.inp were not reproduced
# in the DFM tet fracture faces (facets_f1.inp, facets_f2.inp, etc).
#
read / avs / full_mesh.inp / mo_full
#
read / avs / facets_merged.inp / mo_merge
#
# If the above file exists the next lines can be removed.
#
# Interpolate does not work well on coincident 2D triangulations. C'est la vie.
# To work around this turn the facets into prism volumes by giving them a small
# negative and positive offset and then combine to make prisms. Then you have volume
# cells to interpolate from.
#
#++++++++++++++++++++++++++++++++++++
# EPS_OFFSET  should be set to ~0.1h
#
define / EPS_OFFSET_1  / {-1*eps_offset}
define / EPS_OFFSET_2  /  {eps_offset}
#++++++++++++++++++++++++++++++++++++
offsetsurf / mo_offset_1 / mo_merge / EPS_OFFSET_1
cmo / setatt / mo_offset_1 / imt / 1 0 0 / 1
offsetsurf / mo_offset_2 / mo_merge / EPS_OFFSET_2
cmo / setatt / mo_offset_2 / imt / 1 0 0 / 2
addmesh / merge / mo_offset_1_2 / mo_offset_1 / mo_offset_2
pset / p_bottom / attribute / imt / 1 0 0 / eq / 1
pset / p_top    / attribute / imt / 1 0 0 / eq / 2

extrude / mo_extrude / mo_offset_1_2 / interp / 0 / &
        pset,get,p_bottom / pset,get,p_top

cmo / delete / mo_merge
cmo / delete / mo_offset_1
cmo / delete / mo_offset_2
cmo / delete / mo_offset_1_2
cmo / select / mo_extrude
quality

cmo / addatt / mo_full / mat_interp / vint / scalar / nelements
cmo / setatt / mo_full / mat_interp / 1 0 0 / 2
cmo / setatt / mo_extrude / itetclr / 1 0 0 / 1
interpolate / map / mo_full mat_interp / 1 0 0 / &
                    mo_extrude itetclr
dump / avs / tmp_interpolate.inp / mo_full
cmo / delete / mo_extrude
cmo / select / mo_full
eltset / edelete / mat_interp / eq / 1

cmo / addatt / mo_full / volume / e_area
math / sum / mo_full / area_sum / 1,0,0 / mo_full / e_area

rmpoint / element /  eltset get edelete
rmpoint / compress
# Note: If there are no missed cells, this will return:
# RMPOINT: new point count is            0                                        
# RMPOINT: new element count is          0                                        

cmo / status / brief

cmo / addatt / mo_full / volume / e_area
math / sum / mo_full / area_sum / 1,0,0 / mo_full / e_area
# Note: If there are no missed cells, this MO will be empty and this
# command will return:
# 0 element attribute: e_area
# FATAL ERROR: SUM unable to begin.
# error in command : math/sum/mo_full/area_sum/1,0,0/mo_full/e_area
#
# The attributes that are output in this file could be cleaned up so
# extra unnecessary information is not included.
cmo / DELATT / mo_full / e_area
cmo / DELATT / mo_full / mat_interp
#
# NOTE: If there are no missed cells, mo_full will be an empty (#nodes=0) MO
# No file will be written and LaGriT message will be:
# WARNING: dumpavs             
# WARNING: nnodes=0 nelements = 0
# WARNING: No output
dump / avs / missed_cells_full_mesh.inp / mo_full

cmo / delete / mo_full

finish

"""
    with open('dfm_diagnostics.mlgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush()

    print("Creating dfm_diagonstics.mlgi file: Complete\n")



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
