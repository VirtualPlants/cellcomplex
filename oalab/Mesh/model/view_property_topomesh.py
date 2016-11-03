import numpy as np
import openalea.mesh

from openalea.mesh import PropertyTopomesh
from openalea.mesh.property_topomesh_io import read_ply_property_topomesh
from openalea.mesh.property_topomesh_analysis import compute_topomesh_property, compute_topomesh_vertex_property_from_faces
from openalea.mesh.triangular_mesh import topomesh_to_triangular_mesh

from openalea.container import array_dict


from openalea.deploy.shared_data import shared_data
dirname = shared_data(openalea.mesh)
filename = dirname + "/p194-t4_L1_topomesh.ply"

topomesh = read_ply_property_topomesh(filename)
world.add(topomesh, "topomesh")

compute_topomesh_property(topomesh,'normal',2,normal_method='orientation')
compute_topomesh_vertex_property_from_faces(topomesh,'normal',neighborhood=3,adjacency_sigma=1.2)
compute_topomesh_property(topomesh,'mean_curvature',2)
compute_topomesh_vertex_property_from_faces(topomesh,'mean_curvature',neighborhood=3,adjacency_sigma=1.2)



