# -*- coding: utf-8 -*-
# -*- python -*-
#
#       PropertyTopomesh
#
#       Copyright 2014-2016 INRIA - CIRAD - INRA
#
#       File author(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       File contributor(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenaleaLab Website : http://virtualplants.github.io/
#
###############################################################################

import numpy as np
from scipy import ndimage as nd
from scipy.cluster.vq import vq

from openalea.container import array_dict, PropertyTopomesh
from openalea.cellcomplex.property_topomesh.utils.array_tools import array_unique
from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property, is_triangular
from copy import deepcopy

from time import time

tetra_triangle_edge_list  = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
tetra_triangle_list  = np.array([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
triangle_edge_list  = np.array([[1, 2],[0, 2],[0, 1]])

quad_edge_list  = np.array([[0,1],[1,2],[2,3],[3,0]])

def tetrahedra_topomesh(tetrahedra, positions, **kwargs):

    tetrahedra = np.array(tetrahedra)
    positions = array_dict(positions)

    tetrahedra_triangles = array_unique(np.concatenate(np.sort(tetrahedra[:,tetra_triangle_list])))

    tetrahedra_triangle_edges = tetrahedra_triangles[:,triangle_edge_list]
    tetrahedra_triangle_vectors = positions.values(tetrahedra_triangle_edges[...,1]) - positions.values(tetrahedra_triangle_edges[...,0])
    tetrahedra_triangle_lengths = np.linalg.norm(tetrahedra_triangle_vectors,axis=2)
    tetrahedra_triangle_perimeters = tetrahedra_triangle_lengths.sum(axis=1)

    tetrahedra_edges = array_unique(np.concatenate(tetrahedra_triangles[:,triangle_edge_list],axis=0))

    start_time = time()
    print "--> Generating tetrahedra topomesh"
    triangle_edges = np.concatenate(tetrahedra_triangles[:,triangle_edge_list],axis=0)
    triangle_edge_matching = vq(triangle_edges,tetrahedra_edges)[0]

    tetrahedra_faces = np.concatenate(np.sort(tetrahedra[:,tetra_triangle_list]))
    tetrahedra_triangle_matching = vq(tetrahedra_faces,tetrahedra_triangles)[0]

    tetrahedra_topomesh = PropertyTopomesh(3)
    for c in np.unique(tetrahedra_triangles):
        tetrahedra_topomesh.add_wisp(0,c)
    for e in tetrahedra_edges:
        eid = tetrahedra_topomesh.add_wisp(1)
        for pid in e:
            tetrahedra_topomesh.link(1,eid,pid)
    for t in tetrahedra_triangles:
        fid = tetrahedra_topomesh.add_wisp(2)
        for eid in triangle_edge_matching[3*fid:3*fid+3]:
            tetrahedra_topomesh.link(2,fid,eid)
    for t in tetrahedra:
        cid = tetrahedra_topomesh.add_wisp(3)
        for fid in tetrahedra_triangle_matching[4*cid:4*cid+4]:
            tetrahedra_topomesh.link(3,cid,fid)
    tetrahedra_topomesh.update_wisp_property('barycenter',0,positions.values(np.unique(tetrahedra_triangles)),keys=np.unique(tetrahedra_triangles))        

    end_time = time()
    print "<-- Generating tetrahedra topomesh [",end_time-start_time,"s]"

    return tetrahedra_topomesh

def triangle_topomesh(triangles, positions, **kwargs):

    triangles = np.array(triangles)
    positions = array_dict(positions)

    edges = array_unique(np.sort(np.concatenate(triangles[:,triangle_edge_list],axis=0)))

    triangle_edges = np.sort(np.concatenate(triangles[:,triangle_edge_list]))

    start_time = time()
    print "--> Generating triangle topomesh"

    triangle_edge_matching = vq(triangle_edges,edges)[0]

    triangle_topomesh = PropertyTopomesh(3)
    for c in np.unique(triangles):
        triangle_topomesh.add_wisp(0,c)
    for e in edges:
        eid = triangle_topomesh.add_wisp(1)
        for pid in e:
            triangle_topomesh.link(1,eid,pid)
    for t in triangles:
        fid = triangle_topomesh.add_wisp(2)
        for eid in triangle_edge_matching[3*fid:3*fid+3]:
            triangle_topomesh.link(2,fid,eid)
    triangle_topomesh.add_wisp(3,0)
    for fid in triangle_topomesh.wisps(2):
        triangle_topomesh.link(3,0,fid)
    triangle_topomesh.update_wisp_property('barycenter',0,positions.values(np.unique(triangles)),keys=np.unique(triangles))    

    end_time = time()
    print "<-- Generating triangle topomesh [",end_time-start_time,"s]"

    return triangle_topomesh


def quad_topomesh(quads, positions, faces_as_cells=False, **kwargs):
    quads = np.array(quads)
    positions = array_dict(positions)

    edges = array_unique(np.sort(np.concatenate(quads[:,quad_edge_list],axis=0)))

    quad_edges = np.sort(np.concatenate(quads[:,quad_edge_list]))

    start_time = time()
    print "--> Generating quad topomesh"

    quad_edge_matching = vq(quad_edges,edges)[0]

    quad_topomesh = PropertyTopomesh(3)
    for c in np.unique(quads):
        quad_topomesh.add_wisp(0,c)
    for e in edges:
        eid = quad_topomesh.add_wisp(1)
        for pid in e:
            quad_topomesh.link(1,eid,pid)
    for q in quads:
        fid = quad_topomesh.add_wisp(2)
        for eid in quad_edge_matching[4*fid:4*fid+4]:
            quad_topomesh.link(2,fid,eid)
    if not faces_as_cells:
        quad_topomesh.add_wisp(3,0)
        for fid in quad_topomesh.wisps(2):
            quad_topomesh.link(3,0,fid)
    else:
        for fid in quad_topomesh.wisps(2):
            quad_topomesh.add_wisp(3,fid)
            quad_topomesh.link(3,fid,fid)

    quad_topomesh.update_wisp_property('barycenter',0,positions.values(np.unique(quads)),keys=np.unique(quads))    

    end_time = time()
    print "<-- Generating quad topomesh [",end_time-start_time,"s]"

    return quad_topomesh

def poly_topomesh(polys, positions, faces_as_cells=False, **kwargs):
    polys = np.array(polys)
    positions = array_dict(positions)

    poly_lengths = np.array(map(len,polys))
    poly_edge_list = [np.transpose([np.arange(l),(np.arange(l)+1)%l]) for l in poly_lengths]

    edges = array_unique(np.sort(np.concatenate([np.array(p)[l] for p,l in zip(polys,poly_edge_list)],axis=0)))

    poly_edges = np.sort(np.concatenate([np.array(p)[l] for p,l in zip(polys,poly_edge_list)],axis=0))

    start_time = time()
    print "--> Generating poly topomesh"

    poly_edge_matching = vq(poly_edges,edges)[0]

    poly_topomesh = PropertyTopomesh(3)
    for c in np.unique(polys):
        poly_topomesh.add_wisp(0,c)
    for e in edges:
        eid = poly_topomesh.add_wisp(1)
        for pid in e:
            poly_topomesh.link(1,eid,pid)
    total_poly_length = 0
    for q,l in zip(polys,poly_lengths):
        fid = poly_topomesh.add_wisp(2)
        for eid in poly_edge_matching[total_poly_length:total_poly_length+l]:
            poly_topomesh.link(2,fid,eid)
        total_poly_length += l
    if not faces_as_cells:
        poly_topomesh.add_wisp(3,0)
        for fid in poly_topomesh.wisps(2):
            poly_topomesh.link(3,0,fid)
    else:
        for fid in poly_topomesh.wisps(2):
            poly_topomesh.add_wisp(3,fid)
            poly_topomesh.link(3,fid,fid)
    poly_topomesh.update_wisp_property('barycenter',0,positions.values(np.unique(polys)),keys=np.unique(polys))    

    end_time = time()
    print "<-- Generating poly topomesh [",end_time-start_time,"s]"

    return poly_topomesh


def edge_topomesh(edges, positions, **kwargs):

    positions = array_dict(positions)

    start_time = time()
    print "--> Generating edge topomesh"

    edge_topomesh = PropertyTopomesh(3)
    for c in np.unique(edges):
        edge_topomesh.add_wisp(0,c)
    for e in edges:
        eid = edge_topomesh.add_wisp(1)
        for pid in e:
            edge_topomesh.link(1,eid,pid)
    edge_topomesh.update_wisp_property('barycenter',0,positions.values(np.unique(edges)),keys=np.unique(edges))    

    end_time = time()
    print "<-- Generating edge topomesh [",end_time-start_time,"s]"

    return edge_topomesh

def vertex_topomesh(positions, **kwargs):
    
    positions = array_dict(positions)

    start_time = time()
    print "--> Generating vertex topomesh"

    vertex_topomesh = PropertyTopomesh(3)
    for c in positions.keys():
        vertex_topomesh.add_wisp(0,c)
    vertex_topomesh.update_wisp_property('barycenter',0,positions.values(positions.keys()),keys=positions.keys())    

    end_time = time()
    print "<-- Generating vertex topomesh [",end_time-start_time,"s]"

    return vertex_topomesh


def dual_topomesh(topomesh,degree=2,vertex_positions='barycenter'):

    dual_topomesh = PropertyTopomesh(topomesh.degree())
    
    if degree == 2:
        for d in xrange(3):
            if d<2:
                dual_topomesh._regions[d] = deepcopy(topomesh._borders[2-d])
            if d>0:
                dual_topomesh._borders[d] = deepcopy(topomesh._regions[2-d])
        
        dual_topomesh._borders[3] = dict(zip(dual_topomesh._borders[2].keys(),[[w] for w in dual_topomesh._borders[2].keys()]))
        dual_topomesh._regions[2] = dict(zip(dual_topomesh._borders[2].keys(),[[w] for w in dual_topomesh._borders[2].keys()]))
    
        edges_to_remove = [e for e in dual_topomesh.wisps(1) if len(list(dual_topomesh.borders(1,e)))<2]
        faces_to_remove = [f for f in dual_topomesh.wisps(2) if np.any([e in edges_to_remove for e in dual_topomesh.borders(2,f)])]
        cells_to_remove = faces_to_remove
    
    for e in edges_to_remove:
        dual_topomesh.remove_wisp(1,e)
    for f in faces_to_remove:
        dual_topomesh.remove_wisp(2,f)
    for c in cells_to_remove:
        dual_topomesh.remove_wisp(3,c)
    
    if 'voronoi' in vertex_positions:
        assert is_triangular(topomesh)
        if degree==2:
            from openalea.cellcomplex.property_topomesh.utils.geometry_tools import triangle_geometric_features
            compute_topomesh_property(topomesh,'vertices',2)
            triangles = topomesh.wisp_property('vertices',2).values(list(dual_topomesh.wisps(0)))
            positions = topomesh.wisp_property('barycenter',0)
            if vertex_positions == 'projected_voronoi':
                centers = triangle_geometric_features(triangles,positions,features=['projected_circumscribed_circle_center'])[:,0]
            else:
                centers = triangle_geometric_features(triangles,positions,features=['circumscribed_circle_center'])[:,0]
            dual_positions = array_dict(centers,list(dual_topomesh.wisps(0)))
    else:
        compute_topomesh_property(topomesh,'barycenter',degree)
        dual_positions = array_dict(topomesh.wisp_property('barycenter',degree).values(list(dual_topomesh.wisps(0))),list(dual_topomesh.wisps(0)))
    
   
    dual_topomesh.update_wisp_property('barycenter',0,dual_positions)
    
    return dual_topomesh
