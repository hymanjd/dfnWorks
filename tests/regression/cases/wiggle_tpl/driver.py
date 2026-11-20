# Regression case: wiggle option with a stochastic (Poisson) fracture family.
# Same family as examples/TPL, generation and meshing only. A fracture FRAM
# rejects is shifted by a few h where it is before being re-translated.

from pydfnworks import *
import os

jobname = os.getcwd() + "/output"

DFN = DFNWORKS(jobname, ncpu=4)

DFN.params['domainSize']['value'] = [15, 15, 15]
DFN.params['h']['value'] = 0.1
DFN.params['domainSizeIncrease']['value'] = [0.5, 0.5, 0.5]
DFN.params['keepOnlyLargestCluster']['value'] = True
DFN.params['ignoreBoundaryFaces']['value'] = False
DFN.params['boundaryFaces']['value'] = [1, 1, 0, 0, 0, 0]
DFN.params['seed']['value'] = 2
DFN.params['wiggle']['value'] = True

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        alpha=1.8,
                        min_radius=1.0,
                        max_radius=10.0,
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
DFN.create_network()
DFN.mesh_network(min_dist=1, max_dist=5, max_resolution_factor=10)
