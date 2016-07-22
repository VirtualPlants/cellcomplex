import openalea.mesh

from openalea.mesh import PropertyTopomesh
from openalea.mesh.property_topomesh_io import read_ply_property_topomesh
from openalea.mesh.property_topomesh_analysis import compute_topomesh_property

from openalea.deploy.shared_data import shared_data
dirname = shared_data(openalea.mesh)
filename = dirname + "/p194-t4_L1_topomesh.ply"

topomesh = read_ply_property_topomesh(filename)
world.add(topomesh, "topomesh")

import numpy as np
from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh
from openalea.container import array_dict

compute_topomesh_property(topomesh,'epidermis',1)
compute_topomesh_property(topomesh,'cells',1)

edge_n_cells = np.array(map(len,topomesh.wisp_property('cells',1).values()))
edge_is_cell_edge = edge_n_cells>2
edge_is_cell_edge = edge_is_cell_edge | (edge_n_cells>1)&(topomesh.wisp_property('epidermis',1).values())

edge_is_cell_edge = array_dict(edge_is_cell_edge,list(topomesh.wisps(1)))

mesh,_ = topomesh_to_triangular_mesh(topomesh,1)

edges_to_remove = []
for eid in mesh.edges:
    if not edge_is_cell_edge[eid]:
        edges_to_remove += [eid]
for eid in edges_to_remove:
    del mesh.edges[eid]
    if mesh.edge_data.has_key(eid):
        del mesh.edge_data[eid]
        
world.add(mesh)



