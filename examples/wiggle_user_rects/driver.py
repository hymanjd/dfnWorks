#"""
#   :synopsis: Driver run file for the wiggle example
#   :version: 2.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#"""
#
# Four user-defined rectangles, two of which FRAM rejects where they are placed:
#
#   fracture 3 crosses fracture 1 only 0.04 from where fracture 2 crosses it,
#              closer than h = 0.1 (intersections too close);
#   fracture 4 pokes only 0.03 through fracture 1, so the intersection runs
#              closer than h to fracture 4's own edge.
#
# With wiggle on, DFNGen shifts each rejected fracture by a few h (smallest
# working shift first) until FRAM accepts it. A shift is only used if the
# fracture keeps every intersection it had. The shifts applied are written to
# output/dfnGen_output/wiggled_fractures.dat. Set wiggle to False to see the
# two fractures rejected instead.

from pydfnworks import *
import os

src_path = os.getcwd()
jobname = src_path + "/output"

DFN = DFNWORKS(jobname, ncpu=4)

DFN.params['domainSize']['value'] = [1.0, 1.0, 1.0]
DFN.params['h']['value'] = 0.1
DFN.params['ignoreBoundaryFaces']['value'] = True
DFN.params['wiggle']['value'] = True
DFN.params['wiggleMaxSteps']['value'] = 5

# 1: vertical fracture through the middle of the domain
DFN.add_user_fract(shape='rect',
                   radii=1.0,
                   aspect_ratio=.65,
                   translation=[0, 0, 0],
                   normal_vector=[1, 0, 0],
                   permeability=1.0e-12)

# 2: horizontal fracture crossing fracture 1
DFN.add_user_fract(shape='rect',
                   radii=.6,
                   translation=[0.4, 0, 0.2],
                   normal_vector=[0, 0, 1],
                   permeability=1.0e-12)

# 3: horizontal fracture 0.04 above fracture 2 (rejected without wiggle)
DFN.add_user_fract(shape='rect',
                   radii=.6,
                   translation=[0.4, 0, 0.24],
                   normal_vector=[0, 0, 1],
                   permeability=1.0e-12)

# 4: small fracture whose edge ends 0.03 past fracture 1
#    (rejected without wiggle)
DFN.add_user_fract(shape='rect',
                   radii=.2,
                   translation=[0.17, 0, -0.3],
                   normal_vector=[0, 0, 1],
                   permeability=1.0e-12)

DFN.make_working_directory(delete=True)
DFN.check_input()
DFN.print_domain_parameters()

DFN.create_network()
DFN.mesh_network()
