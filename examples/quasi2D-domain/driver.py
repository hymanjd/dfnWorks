#"""
#   :synopsis: run file for TPL example 
#   :version: 1.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#"""

from pydfnworks import * 

DFN = create_dfn()
DFN.make_working_directory(delete = True)
DFN.check_input()
DFN.create_network()

# DFN.output_report()
DFN.mesh_network(visual_mode = True)

exit()
DFN.dfn_flow()
DFN.dfn_trans()
