import numpy as np
import openalea.mesh

from openalea.mesh import PropertyTopomesh
from openalea.mesh.property_topomesh_io import read_ply_property_topomesh
from openalea.mesh.property_topomesh_analysis import compute_topomesh_property, compute_topomesh_vertex_property_from_faces, compute_topomesh_cell_property_from_faces

from openalea.deploy.shared_data import shared_data
dirname = shared_data(openalea.mesh)
filename = dirname + "/p194-t4_L1_topomesh.ply"

topomesh = read_ply_property_topomesh(filename)

compute_topomesh_property(topomesh,'normal',2,normal_method='orientation')
compute_topomesh_vertex_property_from_faces(topomesh,'normal',neighborhood=3,adjacency_sigma=1.2)
compute_topomesh_property(topomesh,'mean_curvature',2)
compute_topomesh_vertex_property_from_faces(topomesh,'mean_curvature',neighborhood=3,adjacency_sigma=1.2)
compute_topomesh_cell_property_from_faces(topomesh,'normal',aggregate='mean')
compute_topomesh_cell_property_from_faces(topomesh,'mean_curvature',aggregate='mean')
compute_topomesh_cell_property_from_faces(topomesh,'area',aggregate='sum')

world.add(topomesh, "topomesh")

