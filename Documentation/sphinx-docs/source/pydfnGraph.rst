.. _dfnWorks-python-chapter-dfnGraph:

pydfnworks: dfnGraph
========================================

DFN Class functions used in graph analysis and pipe-network simulations

General Graph Functions
^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.dfnGraph.dfn2graph
    :members: create_graph, dump_json_graph, load_json_graph, plot_graph, dump_fractures, add_fracture_source, add_fracture_target

.. automodule:: pydfnworks.dfnGraph.pruning
    :members:  k_shortest_paths_backbone, greedy_edge_disjoint





Graph-Based Flow and Transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.dfnGraph.graph_flow
    :members: run_graph_flow

.. automodule:: pydfnworks.dfnGraph.graph_transport
    :members: run_graph_transport


BiGraFV
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``DFN.bigrafv()`` converts the bipartite graph of a DFN into PFLOTRAN unstructured grid input files (``.uge``, boundary ``.ex`` files and ``dfn_properties.h5``).

.. automodule:: pydfnworks.dfnGraph.bigrafv.bigrafv
    :members: bigrafv

.. automodule:: pydfnworks.dfnGraph.bigrafv.extract_graph_boundaries_from_dfn
    :members: extract_graph_boundaries_from_dfn, parse_boundary_ex

.. automodule:: pydfnworks.dfnGraph.bigrafv.graph_uge
    :members: convert_graph_to_data_frames, make_cell_data_frame, make_connection_data_frame, add_boundary_nodes_as_cells, add_boundary_nodes_as_conns

.. automodule:: pydfnworks.dfnGraph.bigrafv.compute_graph_to_pflotran_geometries
    :members: compute_graph_to_pflotran_geometries

.. automodule:: pydfnworks.dfnGraph.bigrafv.write_graph_uge
    :members: write_graph_uge

.. automodule:: pydfnworks.dfnGraph.bigrafv.write_boundary_ex
    :members: write_boundary_ex, write_boundary_ex_files, write_ex

.. automodule:: pydfnworks.dfnGraph.bigrafv.dump_h5_file_for_graph_perm
    :members: dump_h5_file_for_graph_perm
