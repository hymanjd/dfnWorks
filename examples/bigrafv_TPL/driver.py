#"""
#   :synopsis: Driver run file for BiGraFV TPL example
#   :version: 1.0
#   :maintainer: Collin Sutton
#.. moduleauthor:: Collin Sutton <collin.sutton@mines.edu>
#"""

from pydfnworks import *
import os

src_path = os.getcwd()
jobname = src_path + "/output"
# pflotran_graph_pressure.in solves for steady-state pressure and writes a checkpoint
# pflotran_graph_tracer.in restarts from that checkpoint and runs a tracer
dfnFlow_file = src_path + '/pflotran_graph_pressure.in'
restart_file = src_path + '/pflotran_graph_tracer.in'

DFN = DFNWORKS(jobname, dfnFlow_file=dfnFlow_file, ncpu=8)

DFN.params['domainSize']['value'] = [100, 50, 50]
DFN.params['h']['value'] = 0.2
DFN.params['domainSizeIncrease']['value'] = [5, 5, 5]
DFN.params['keepOnlyLargestCluster']['value'] = True
DFN.params['ignoreBoundaryFaces']['value'] = False
DFN.params['stopCondition']['value'] = 1
DFN.params['boundaryFaces']['value'] = [1, 1, 0, 0, 0, 0]
DFN.params['seed']['value'] = 10

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        alpha=1.8,
                        min_radius=2.0,
                        max_radius=15.0,
                        orientation_distribution="fisher",
                        kappa=1.0,
                        theta=0.0,
                        phi=0.0,
                        p32=1,
                        hy_variable='aperture',
                        hy_function='correlated',
                        number_of_points=8,
                        hy_params={
                            "alpha": 10**-5,
                            "beta": 0.5
                        })

DFN.make_working_directory(delete=True)
DFN.check_input()
DFN.print_domain_parameters()

DFN.create_network()

# Full DFN (optional): mesh the network and run PFLOTRAN on the mesh.
# Q_dfn is the reference flow rate for the bias correction below.
# BiGraFV then writes over these files. Uncomment to run.
# DFN.mesh_network()
# DFN.lagrit2pflotran()
# DFN.pflotran(restart=True, restart_file=restart_file)
# DFN.parse_pflotran_vtk_python()
# DFN.pflotran_cleanup()
# Q_dfn = DFN.effective_perm(inflow_pressure=2e6,
#                            outflow_pressure=1e6,
#                            boundary_file="boundary_left_w.ex",
#                            direction="x",
#                            darcy_vel_file="darcyvel_001.dat")

# BiGraFV - run PFLOTRAN on the bipartite graph
G = DFN.create_graph("bipartite", inflow="left", outflow="right")
DFN.bigrafv(G)
DFN.pflotran(restart=True, restart_file=restart_file)
DFN.pflotran_cleanup()
Q_graph = DFN.effective_perm(inflow_pressure=2e6,
                             outflow_pressure=1e6,
                             boundary_file="boundary_left_w.ex",
                             direction="x",
                             darcy_vel_file="darcyvel_001.dat")

# Bias correction (optional): scale the permeabilities by alpha = Q_dfn / Q_graph
# and rerun BiGraFV with the tracer so its flow rate matches the full DFN.
# Requires the full DFN run above. Uncomment to run.
# alpha = Q_dfn / Q_graph
# DFN.go_home()
# DFN.perm = DFN.perm * alpha
# DFN.bigrafv(G)
# DFN.pflotran(restart=True, restart_file=restart_file)
# DFN.pflotran_cleanup()
# Q_corrected = DFN.effective_perm(inflow_pressure=2e6,
#                                  outflow_pressure=1e6,
#                                  boundary_file="boundary_left_w.ex",
#                                  direction="x",
#                                  darcy_vel_file="darcyvel_001.dat")
# print(f"--> alpha = {alpha:.8e}")
# print(f"--> Relative error corrected = {(Q_corrected - Q_dfn) / Q_dfn:.8e}")
