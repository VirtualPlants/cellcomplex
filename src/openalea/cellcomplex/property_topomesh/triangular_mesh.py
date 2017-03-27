# -*- coding: utf-8 -*-
# -*- python -*-
#
#       TriangularMesh
#
#       Copyright 2015-2016 INRIA - CIRAD - INRA
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
from openalea.container import array_dict

from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property, is_triangular
from openalea.cellcomplex.property_topomesh.property_topomesh_extraction import star_interface_topomesh
from openalea.cellcomplex.triangular_mesh import TriangularMesh

from time import time
from copy import deepcopy

def topomesh_to_triangular_mesh(input_topomesh, degree=3, wids=None, coef=1.0, mesh_center=None, epidermis=False, cell_edges=False, property_name=None, property_degree=None):

    topomesh = deepcopy(input_topomesh)
    if not is_triangular(topomesh) and (degree>1):
        topomesh = star_interface_topomesh(topomesh)




    start_time = time()
    print "--> Creating triangular mesh"

    triangular_mesh = TriangularMesh()

    if property_name is not None:
        if property_degree is None:
            property_degree = degree
        try:
            if not topomesh.has_wisp_property(property_name,degree=property_degree,is_computed=True):
                compute_topomesh_property(topomesh,property_name,degree=property_degree)
            assert len(topomesh.wisp_property(property_name,property_degree).keys()) == topomesh.nb_wisps(property_degree)
        except:
            property_name = None

    if mesh_center is None:
        mesh_center = topomesh.wisp_property('barycenter',0).values().mean(axis=0)
    else:
        mesh_center = np.array(mesh_center)


    considered_start_time = time()
    if wids is None:
        wids = list(topomesh.wisps(degree))

    considered_elements = {}
    for d in xrange(4):
        if d == degree:
            considered_elements[d] = list(wids)
        elif d<degree:
            considered_elements[d] = list(np.unique(np.concatenate([list(topomesh.borders(degree,w,degree-d)) for w in wids])))
        else:
            considered_elements[d] = list(np.unique(np.concatenate([list(topomesh.regions(degree,w,d-degree)) for w in wids])))

    if degree>1:
        compute_topomesh_property(topomesh,'vertices',degree=2)
        compute_topomesh_property(topomesh,'triangles',degree=3)
        compute_topomesh_property(topomesh,'cells',degree=2)

        cell_triangles = np.concatenate(topomesh.wisp_property('triangles',3).values(considered_elements[3])).astype(int)
        # cell_triangles = [t for t in cell_triangles if t in considered_elements[2]]

    considered_end_time = time()
    print "  --> Filtering out mesh elements    [",considered_end_time - considered_start_time,"s]"


    if degree == 3:
        if property_name is not None:
            property_data = topomesh.wisp_property(property_name,property_degree).values(considered_elements[3])
        else:
            property_data = np.array(considered_elements[3])


        vertices_positions = []
        triangle_vertices = []
        triangle_topomesh_cells = []
        vertices_topomesh_vertices = []
        # triangle_topomesh_triangles = []

        if property_data.ndim == 1 or property_degree<3:
            for c in considered_elements[3]:
                if len(list(topomesh.borders(3,c,3)))>0:
                    cell_center = topomesh.wisp_property('barycenter',0).values(list(topomesh.borders(3,c,3))).mean(axis=0)
                    cell_vertices_position = cell_center + coef*(topomesh.wisp_property('barycenter',0).values(list(topomesh.borders(3,c,3)))-cell_center) - mesh_center
                    cell_vertices_index = array_dict(len(vertices_positions) + np.arange(len(list(topomesh.borders(3,c,3)))),list(topomesh.borders(3,c,3)))
                    vertices_positions += list(cell_vertices_position)
                    vertices_topomesh_vertices += list(topomesh.borders(3,c,3))
                    triangle_vertices += list(cell_vertices_index.values(topomesh.wisp_property('vertices',2).values(topomesh.wisp_property('triangles',3)[c])))
                    triangle_topomesh_cells += list(c*np.ones_like(topomesh.wisp_property('triangles',3)[c]))
                    # triangle_topomesh_triangles += topomesh.wisp_property('triangles',3)[c]
            print len(cell_triangles),"Cell triangles"
            print len(triangle_vertices),"Triangle vertices"
            vertices_positions = array_dict(vertices_positions,np.arange(len(vertices_positions)))
            vertices_topomesh_vertices = array_dict(vertices_topomesh_vertices,np.arange(len(vertices_positions)))
            if epidermis:
                compute_topomesh_property(topomesh,'epidermis',2)
                epidermis_triangles = topomesh.wisp_property('epidermis',2).values(cell_triangles)
                triangle_vertices = array_dict(np.array(triangle_vertices)[epidermis_triangles],np.arange(len(cell_triangles[epidermis_triangles])))
                triangle_topomesh_cells = array_dict(np.array(triangle_topomesh_cells)[epidermis_triangles],np.arange(len(cell_triangles[epidermis_triangles])))
                triangle_topomesh_triangles = array_dict(cell_triangles[epidermis_triangles],np.arange(len(cell_triangles[epidermis_triangles])))
            else:
                triangle_vertices = array_dict(triangle_vertices,np.arange(len(cell_triangles)))
                triangle_topomesh_cells = array_dict(triangle_topomesh_cells,np.arange(len(cell_triangles)))
                triangle_topomesh_triangles = array_dict(cell_triangles,np.arange(len(cell_triangles)))
            edge_topomesh_edges = {}

            triangular_mesh.points = vertices_positions.to_dict()
            triangular_mesh.triangles = triangle_vertices.to_dict()
            if property_name is not None:  
                if property_degree == 2:
                    triangle_topomesh_triangle_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(triangle_topomesh_triangles.values()),triangle_topomesh_triangles.keys())
                    triangular_mesh.triangle_data = triangle_topomesh_triangle_property.to_dict()
                elif property_degree == 0:
                    vertex_topomesh_vertex_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(vertices_topomesh_vertices.values()),vertices_topomesh_vertices.keys())
                    triangular_mesh.point_data = vertex_topomesh_vertex_property.to_dict()
                elif property_degree == 3:
                    triangle_topomesh_cell_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(triangle_topomesh_cells.values()),triangle_topomesh_cells.keys())
                    triangular_mesh.triangle_data = triangle_topomesh_cell_property.to_dict()
            else:
                triangular_mesh.triangle_data = triangle_topomesh_cells.to_dict()
        else:
            for c in considered_elements[3]:
                if len(list(topomesh.borders(3,c,3)))>0:
                    cell_center = topomesh.wisp_property('barycenter',0).values(list(topomesh.borders(3,c,3))).mean(axis=0)
                    vertices_positions += [cell_center]

            vertices_positions = array_dict(vertices_positions,np.array([c for c in topomesh.wisps(3) if len(list(topomesh.borders(3,c,3)))>0]))
            vertices_topomesh_vertices = {}
            edge_topomesh_edges = {}
            triangle_topomesh_triangles = {}
            triangle_topomesh_cells = {}

            cell_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(vertices_positions.keys()),vertices_positions.keys())

            triangular_mesh.points = vertices_positions.to_dict()
            triangular_mesh.point_data = cell_property
            triangular_mesh.triangles = {}

    elif degree == 2:
        vertices_positions = []
        triangle_vertices = []
        vertices_topomesh_vertices = []
        for t in cell_triangles:
            triangle_center = topomesh.wisp_property('barycenter',0).values(list(topomesh.borders(2,t,2))).mean(axis=0)
            triangle_vertices_position = triangle_center + coef*(topomesh.wisp_property('barycenter',0).values(list(topomesh.borders(2,t,2)))-triangle_center) - mesh_center
            triangle_vertices_index = array_dict(len(vertices_positions) + np.arange(3),list(topomesh.borders(2,t,2)))
            vertices_positions += list(triangle_vertices_position)
            vertices_topomesh_vertices += list(topomesh.borders(2,t,2))
            triangle_vertices += list([triangle_vertices_index.values(topomesh.wisp_property('vertices',2)[t])])
        vertices_positions = array_dict(vertices_positions,np.arange(len(vertices_positions)))
        vertices_topomesh_vertices = array_dict(vertices_topomesh_vertices,np.arange(len(vertices_positions)))
        triangle_topomesh_cells = np.concatenate([c*np.ones_like([t for t in topomesh.wisp_property('triangles',3)[c] if t in considered_elements[2]]) for c in considered_elements[3]]).astype(int)
        if epidermis:
            compute_topomesh_property(topomesh,'epidermis',2)
            epidermis_triangles = topomesh.wisp_property('epidermis',2).values(cell_triangles)
            triangle_vertices = array_dict(np.array(triangle_vertices)[epidermis_triangles],np.arange(len(cell_triangles[epidermis_triangles])))
            triangle_topomesh_cells = array_dict(np.array(triangle_topomesh_cells)[epidermis_triangles],np.arange(len(cell_triangles[epidermis_triangles])))
            triangle_topomesh_triangles = array_dict(cell_triangles[epidermis_triangles],np.arange(len(cell_triangles[epidermis_triangles])))
        else:
            triangle_vertices = array_dict(triangle_vertices,np.arange(len(cell_triangles)))
            triangle_topomesh_cells = array_dict(triangle_topomesh_cells,np.arange(len(cell_triangles)))
            triangle_topomesh_triangles = array_dict(cell_triangles,np.arange(len(cell_triangles)))
        edge_topomesh_edges = {}

        triangular_mesh.points = vertices_positions.to_dict()
        triangular_mesh.triangles = triangle_vertices.to_dict()
        if property_name is not None:
            if property_degree == 2:
                triangle_topomesh_triangle_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(triangle_topomesh_triangles.values()),triangle_topomesh_triangles.keys())
                triangular_mesh.triangle_data = triangle_topomesh_triangle_property.to_dict()
            elif property_degree == 0:
                vertex_topomesh_vertex_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(vertices_topomesh_vertices.values()),vertices_topomesh_vertices.keys())
                triangular_mesh.point_data = vertex_topomesh_vertex_property.to_dict()
            elif property_degree == 3:
                triangle_topomesh_cell_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(triangle_topomesh_cells.values()),triangle_topomesh_cells.keys())
                triangular_mesh.triangle_data = triangle_topomesh_cell_property.to_dict()
        else:
            triangular_mesh.triangle_data = triangle_topomesh_cells.to_dict()

    elif degree == 1:
        compute_topomesh_property(topomesh,'vertices',degree=1)
        vertices_positions = array_dict(topomesh.wisp_property('barycenter',0).values(considered_elements[0]),considered_elements[0])
        edge_vertices = array_dict(topomesh.wisp_property('vertices',1).values(considered_elements[1]),considered_elements[1])
        triangular_mesh.points = vertices_positions.to_dict()
        triangular_mesh.edges = edge_vertices.to_dict()

        if property_name is not None:
            if property_degree == 1:
                edge_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(considered_elements[1]),considered_elements[1])
                triangular_mesh.edge_data = edge_property.to_dict()
        triangle_topomesh_cells = {}
        triangle_topomesh_triangles = {}
        edge_topomesh_edges = dict(zip(triangular_mesh.edges.keys(),triangular_mesh.edges.keys()))
        vertices_topomesh_vertices = {}

        if cell_edges:
            compute_topomesh_property(topomesh,'epidermis',1)
            compute_topomesh_property(topomesh,'cells',1)
    
            edge_n_cells = np.array(map(len,topomesh.wisp_property('cells',1).values()))
            edge_is_cell_edge = edge_n_cells>2
            edge_is_cell_edge = edge_is_cell_edge | (edge_n_cells>1)&(topomesh.wisp_property('epidermis',1).values())
            edge_is_cell_edge = array_dict(edge_is_cell_edge,list(topomesh.wisps(1)))

            edges_to_remove = []
            for eid in triangular_mesh.edges:
                if not edge_is_cell_edge[eid]:
                    edges_to_remove += [eid]
            for eid in edges_to_remove:
                del triangular_mesh.edges[eid]
                del edge_topomesh_edges[eid]
                if triangular_mesh.edge_data.has_key(eid):
                    del triangular_mesh.edge_data[eid]

    elif degree == 0:
        vertices_positions = array_dict(topomesh.wisp_property('barycenter',0).values(considered_elements[0]),considered_elements[0])
        triangular_mesh.points = vertices_positions.to_dict()

        if property_name is not None:
            if property_degree == 0:
                vertex_property = array_dict(topomesh.wisp_property(property_name,property_degree).values(considered_elements[0]),considered_elements[0])
                triangular_mesh.point_data = vertex_property.to_dict()
        triangle_topomesh_cells = {}
        triangle_topomesh_triangles = {}
        edge_topomesh_edges = {}
        vertices_topomesh_vertices = dict(zip(triangular_mesh.points.keys(),triangular_mesh.points.keys()))

    mesh_element_matching = {}
    mesh_element_matching[0] = vertices_topomesh_vertices
    mesh_element_matching[1] = edge_topomesh_edges
    mesh_element_matching[2] = triangle_topomesh_triangles
    mesh_element_matching[3] = triangle_topomesh_cells

    end_time = time()
    print "<-- Creating triangular mesh [",end_time - start_time,"s]"

    return triangular_mesh, mesh_element_matching