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

from copy import deepcopy

from openalea.container import array_dict

from openalea.cellcomplex.property_topomesh import PropertyTopomesh
from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property, is_triangular

from openalea.cellcomplex.property_topomesh.utils.geometry_tools import triangle_geometric_features


def epidermis_topomesh(topomesh,cells=None):
    
    compute_topomesh_property(topomesh,'epidermis',3)
    compute_topomesh_property(topomesh,'epidermis',1)
    compute_topomesh_property(topomesh,'epidermis',0)
    
    if cells is None:
        epidermis_topomesh = deepcopy(topomesh)
    else:
        faces = np.array(np.unique(np.concatenate([np.array(list(topomesh.borders(3,c))) for c in cells])),int)
        edges = np.array(np.unique(np.concatenate([np.array(list(topomesh.borders(2,t))) for t in faces])),int)
        vertices = np.array(np.unique(np.concatenate([np.array(list(topomesh.borders(1,e))) for e in edges])),int)
        epidermis_topomesh = PropertyTopomesh(3)
        vertices_to_pids = {}
        for v in vertices:
            pid = epidermis_topomesh.add_wisp(0,v)
            vertices_to_pids[v] = pid
        edges_to_eids = {}
        for e in edges:
            eid = epidermis_topomesh.add_wisp(1,e)
            edges_to_eids[e] = eid
            for v in topomesh.borders(1,e):
                epidermis_topomesh.link(1,eid,vertices_to_pids[v])
        faces_to_fids = {}
        for f in faces:
            fid = epidermis_topomesh.add_wisp(2,f)
            faces_to_fids[f] = fid
            for e in topomesh.borders(2,f):
                epidermis_topomesh.link(2,fid,edges_to_eids[e])
        for c in cells:
            cid = epidermis_topomesh.add_wisp(3,c)
            for f in topomesh.borders(3,c):
                epidermis_topomesh.link(3,cid,faces_to_fids[f])
    
    vertices_to_remove = []
    for v in epidermis_topomesh.wisps(0):
        if not topomesh.wisp_property('epidermis',0)[v]:
            vertices_to_remove.append(v)
    for v in vertices_to_remove:
        epidermis_topomesh.remove_wisp(0,v)
    edges_to_remove = []
    for e in epidermis_topomesh.wisps(1):
        if not topomesh.wisp_property('epidermis',1)[e]:
            edges_to_remove.append(e)
    for e in edges_to_remove:
        epidermis_topomesh.remove_wisp(1,e)
    faces_to_remove = []
    for f in epidermis_topomesh.wisps(2):
        if not topomesh.wisp_property('epidermis',2)[f]:
            faces_to_remove.append(f)
    for f in faces_to_remove:
        epidermis_topomesh.remove_wisp(2,f)
    cells_to_remove = []
    for c in epidermis_topomesh.wisps(3):
        if not topomesh.wisp_property('epidermis',3)[c]:
            cells_to_remove.append(c)
    for c in cells_to_remove:
        epidermis_topomesh.remove_wisp(3,c)
    epidermis_topomesh.update_wisp_property('barycenter',0,topomesh.wisp_property('barycenter',0).values(list(epidermis_topomesh.wisps(0))),keys=np.array(list(epidermis_topomesh.wisps(0))))
    return epidermis_topomesh


def cut_surface_topomesh(input_topomesh, z_cut=0, below=True):

    topomesh = deepcopy(input_topomesh)    

    compute_topomesh_property(topomesh,'vertices',2)

    if below:
        triangle_below = array_dict(np.all(topomesh.wisp_property('barycenter',0).values(topomesh.wisp_property('vertices',2).values())[...,2] < z_cut,axis=1),list(topomesh.wisps(2)))
    else:
        triangle_below = array_dict(np.all(topomesh.wisp_property('barycenter',0).values(topomesh.wisp_property('vertices',2).values())[...,2] > z_cut,axis=1),list(topomesh.wisps(2)))
    topomesh.update_wisp_property('below',2,triangle_below)

    triangles_to_remove = [t for t in topomesh.wisps(2) if triangle_below[t]]
    for t in triangles_to_remove:
        topomesh.remove_wisp(2,t)

    topomesh = clean_topomesh(topomesh)

    compute_topomesh_property(topomesh,'triangles',1)
    compute_topomesh_property(topomesh,'vertices',1)
    compute_topomesh_property(topomesh,'length',1)

    topomesh.update_wisp_property('boundary',1,array_dict((np.array(map(len,topomesh.wisp_property('triangles',1).values()))==1).astype(int),list(topomesh.wisps(1))))

    boundary_edges = np.array(list(topomesh.wisps(1)))[topomesh.wisp_property('boundary',1).values()==1]
    boundary_vertices = np.unique(topomesh.wisp_property('vertices',1).values(boundary_edges))

    # z_offset = topomesh.wisp_property('barycenter',0).values()[:,2].std()/8.
    z_offset = np.percentile(topomesh.wisp_property('length',1).values(),10)
    iso_z_positions = np.array([np.concatenate([topomesh.wisp_property('barycenter',0)[v][:2],[z_cut+(1-2*below)*z_offset]]) if v in boundary_vertices else  topomesh.wisp_property('barycenter',0)[v] for v in topomesh.wisps(0)])
    topomesh.update_wisp_property('barycenter',0,array_dict(iso_z_positions,list(topomesh.wisps(0))))

    return topomesh


def clean_topomesh(input_topomesh, clean_properties=False, degree=2):

    topomesh = deepcopy(input_topomesh)

    cells_to_remove = [w for w in topomesh.wisps(3) if topomesh.nb_borders(3,w)==0]
    for w in cells_to_remove:
        topomesh.remove_wisp(3,w)

    if degree > 2:
        triangles_to_remove = [w for w in topomesh.wisps(2) if topomesh.nb_regions(2,w)==0]
        for w in triangles_to_remove:
            topomesh.remove_wisp(2,w)

    edges_to_remove = [w for w in topomesh.wisps(1) if topomesh.nb_regions(1,w)==0]
    for w in edges_to_remove:
        topomesh.remove_wisp(1,w)
        
    vertices_to_remove = [w for w in topomesh.wisps(0) if topomesh.nb_regions(0,w)==0]
    for w in vertices_to_remove:
        topomesh.remove_wisp(0,w)

    if clean_properties:
        for degree in xrange(4):
            for property_name in topomesh.wisp_property_names(degree):
                topomesh.update_wisp_property(property_name,degree,dict(zip(list(topomesh.wisps(degree)),topomesh.wisp_property(property_name,degree).values(list(topomesh.wisps(degree))))))

    return topomesh
    

def cell_topomesh(input_topomesh, cells=None, copy_properties=False):
    from time import time
    start_time = time()

    #topomesh = PropertyTopomesh(topomesh=input_topomesh)
    topomesh = PropertyTopomesh(3)

    if cells is None:
        cells = set(topomesh.wisps(3))
    else:
        cells = set(cells)

    faces = set()
    for c in cells:
        topomesh._borders[3][c] = deepcopy(input_topomesh._borders[3][c])
        faces = faces.union(set(topomesh._borders[3][c]))

    edges = set()
    for f in faces:
        topomesh._borders[2][f] = deepcopy(input_topomesh._borders[2][f])
        topomesh._regions[2][f] = deepcopy(list(set(input_topomesh._regions[2][f]).intersection(cells)))
        edges = edges.union(set(topomesh._borders[2][f]))

    vertices = set()
    for e in edges:
        topomesh._borders[1][e] = deepcopy(input_topomesh._borders[1][e])
        topomesh._regions[1][e] = deepcopy(list(set(input_topomesh._regions[1][e]).intersection(faces)))
        vertices = vertices.union(set(topomesh._borders[1][e]))

    for v in vertices:
        topomesh._regions[0][v] = deepcopy(list(set(input_topomesh._regions[0][v]).intersection(edges)))


    if copy_properties:
        for degree in xrange(4):
            for property_name in input_topomesh.wisp_property_names(degree):
                input_property = input_topomesh.wisp_property(property_name,degree)
                property_keys = np.intersect1d(input_property.keys(),list(topomesh.wisps(degree)))
                if len(property_keys)>0:
                    topomesh.update_wisp_property(property_name,degree,array_dict(input_topomesh.wisp_property(property_name,degree).values(property_keys),property_keys))

    topomesh.update_wisp_property('barycenter',0,array_dict(input_topomesh.wisp_property('barycenter',0).values(list(vertices)),list(vertices)))

    # cells_to_remove = [c for c in topomesh.wisps(3) if not c in cells]
    # cells_to_remove = list(set(topomesh.wisps(3)).difference(set(cells)))

    # for c in cells_to_remove:
    #     topomesh.remove_wisp(3,c)

    # faces_to_remove = [w for w in topomesh.wisps(2) if topomesh.nb_regions(2,w)==0]
    # for w in faces_to_remove:
    #     topomesh.remove_wisp(2,w)

    # edges_to_remove = [w for w in topomesh.wisps(1) if topomesh.nb_regions(1,w)==0]
    # for w in edges_to_remove:
    #     topomesh.remove_wisp(1,w)
        
    # vertices_to_remove = [w for w in topomesh.wisps(0) if topomesh.nb_regions(0,w)==0]
    # for w in vertices_to_remove:
    #     topomesh.remove_wisp(0,w)

    end_time = time()
    #end_time = time()
    print "<-- Extracting cell topomesh     [",end_time-start_time,"s]"

    return topomesh



def star_interface_topomesh(topomesh, inner_interfaces=True, verbose=False):
    from time import time

    triangle_edge_list  = np.array([[1, 2],[0, 2],[0, 1]])

    triangular_topomesh = PropertyTopomesh(3)
    triangle_vertex_positions = {}

    for v in topomesh.wisps(0):
        triangular_topomesh.add_wisp(0,v)
        triangle_vertex_positions[v] = topomesh.wisp_property('barycenter',0)[v]

    for e in topomesh.wisps(1):
        triangular_topomesh.add_wisp(1,e)
        for v in topomesh.borders(1,e):
            triangular_topomesh.link(1,e,v)

    for c in topomesh.wisps(3):
        triangular_topomesh.add_wisp(3,c)


    compute_topomesh_property(topomesh,'regions',2)
    compute_topomesh_property(topomesh,'edges',2)
    compute_topomesh_property(topomesh,'vertices',1)
    compute_topomesh_property(topomesh,'vertices',2)

    face_centers = {} 
    face_triangles = {}

    start_time = time()
    print "--> Triangulating Interfaces"
    for interface in topomesh.wisps(2):
        if interface%100 == 0:
            interface_start_time = time()

        if topomesh.nb_borders(2,interface)>0:
            interface_cells = topomesh.wisp_property('regions',2)[interface]
            interface_edges = topomesh.wisp_property('vertices',1).values(topomesh.wisp_property('edges',2)[interface])
            interface_vertices = np.unique(interface_edges)

            if (len(interface_vertices)>2) and (inner_interfaces or (len(interface_cells) == 1)):

                interface_positions = array_dict(topomesh.wisp_property('barycenter',0).values(interface_vertices),interface_vertices)
                interface_center = interface_positions.values().mean(axis=0)

                center_pid = triangular_topomesh.add_wisp(0)
                triangle_vertex_positions[center_pid] = interface_center

                face_centers[interface] = center_pid
                face_triangles[interface] = []

                vertex_center_edges = {}
                for v in interface_vertices:
                    eid = triangular_topomesh.add_wisp(1)
                    triangular_topomesh.link(1,eid,v)
                    triangular_topomesh.link(1,eid,center_pid)
                    vertex_center_edges[v] = eid

                for e in topomesh.borders(2,interface):
                    fid = triangular_topomesh.add_wisp(2)
                    face_triangles[interface] += [fid]
                    triangular_topomesh.link(2,fid,e)
                    for v in topomesh.borders(1,e):
                        triangular_topomesh.link(2,fid,vertex_center_edges[v])
                    for cid in interface_cells:
                        triangular_topomesh.link(3,cid,fid)

            if verbose:
                if interface%100 == 0:
                    interface_end_time = time()
                    # print "  --> Interface ",interface," / ",topomesh.nb_wisps(2),' ',interface_cells,'     [',interface_end_time-interface_start_time,'s]'
                    print "  --> Interface ",interface," / ",topomesh.nb_wisps(2),'     [',(interface_end_time-interface_start_time),'s]'
    end_time = time()

    # for property_name in topomesh.wisp_property_names(0):
    #     try:
    #         center_property = [topomesh.wisp_property(property_name,0).values(topomesh.wisp_property('vertices',2)[f]).mean(axis=0) for f in face_centers.keys()]
    #     except:
    #         center_property = [topomesh.wisp_property(property_name,0).values(topomesh.wisp_property('vertices',2)[f])[0] for f in face_centers.keys()]
    #     vertex_property = array_dict(list(topomesh.wisp_property(property_name,0).values(list(topomesh.wisps(0))))+center_property,list(topomesh.wisps(0))+[face_centers[f] for f in face_centers.keys()])
    #     triangular_topomesh.update_wisp_property(property_name,0,vertex_property)

    for property_name in topomesh.wisp_property_names(2):
        triangle_faces = np.concatenate([[f for t in face_triangles[f]] for f in face_triangles.keys()])
        triangle_keys = np.concatenate([face_triangles[f] for f in face_triangles.keys()])
        triangle_property = array_dict(topomesh.wisp_property(property_name,2).values(triangle_faces),triangle_keys)
        triangular_topomesh.update_wisp_property(property_name,2,triangle_property)


    for property_name in topomesh.wisp_property_names(3):
        triangular_topomesh.update_wisp_property(property_name,3,topomesh.wisp_property(property_name,3))


    print "--> Triangulating Interfaces  [",end_time-start_time,"s]"
    triangular_topomesh.update_wisp_property('barycenter',degree=0,values=triangle_vertex_positions)
    return triangular_topomesh


def surface_dual_topomesh(topomesh, vertex_placement='center', exterior_vertex=1, exterior_distance=0, face_positions=None, vertices_as_cells=None):
    if vertices_as_cells is None:
        vertices_as_cells = topomesh.nb_wisps(3) != 1

    dual_topomesh = PropertyTopomesh(3)

    for f in topomesh.wisps(2):
        dual_topomesh.add_wisp(0,f)

    for e in topomesh.wisps(1):
        dual_topomesh.add_wisp(1,e)
        for f in topomesh.regions(1,e):
            dual_topomesh.link(1,e,f)

    for v in topomesh.wisps(0):
        if v != exterior_vertex:
            dual_topomesh.add_wisp(2,v)
            for e in topomesh.regions(0,v):
                dual_topomesh.link(2,v,e)

            if vertices_as_cells:
                dual_topomesh.add_wisp(3,v)
                dual_topomesh.link(3,v,v)
            else:
                for c in topomesh.regions(0,v,3):
                    if not dual_topomesh.has_wisp(3,c):
                        dual_topomesh.add_wisp(3,c)
                    dual_topomesh.link(3,c,v)

    for degree in [0,1,2]:
        for property_name in topomesh.wisp_property_names(degree):
            dual_topomesh.update_wisp_property(property_name,2-degree,topomesh.wisp_property(property_name,degree))

    if vertices_as_cells:
        for property_name in topomesh.wisp_property_names(0):
            dual_topomesh.update_wisp_property(property_name,3,topomesh.wisp_property(property_name,0))
    else:
        for property_name in topomesh.wisp_property_names(3):
            dual_topomesh.update_wisp_property(property_name,3,topomesh.wisp_property(property_name,3))


    compute_topomesh_property(topomesh,'barycenter',2)
    face_centers = topomesh.wisp_property('barycenter',2)
    positions = topomesh.wisp_property('barycenter',0)

    if is_triangular(topomesh):
        compute_topomesh_property(topomesh,'vertices',2)
        triangle_ids = np.array(list(topomesh.wisps(2)))
        triangles = topomesh.wisp_property('vertices',2).values(triangle_ids)

        if exterior_vertex is not None:
            exterior_triangle_ids = triangle_ids[np.where(np.any(triangles==exterior_vertex,axis=1))]
            exterior_triangles = triangles[np.where(np.any(triangles==exterior_vertex,axis=1))]
            exterior_edges = [t[t!=exterior_vertex] for t in exterior_triangles]
            exterior_edge_ids = np.concatenate([[e for e in topomesh.borders(2,t) if not exterior_vertex in topomesh.borders(1,e)] for t in exterior_triangle_ids])
            exterior_edge_triangle_ids = np.concatenate([[t for t in topomesh.regions(1,e) if not exterior_vertex in topomesh.borders(2,t,2)] for e in exterior_edge_ids])

            triangle_ids = triangle_ids[np.where(np.all(triangles!=exterior_vertex,axis=1))]
            triangles = triangles[np.where(np.all(triangles!=exterior_vertex,axis=1))]
            
        else:
            exterior_triangles = []
            exterior_triangle_ids = []
            exterior_edges = []
            exterior_edge_triangle_ids = []
            exterior_centers = []

        if vertex_placement == 'voronoi':
            centers = triangle_geometric_features(triangles,positions,['circumscribed_circle_center'])[:,0]
        elif vertex_placement == 'projected_voronoi':
            centers = triangle_geometric_features(triangles,positions,['projected_circumscribed_circle_center'])[:,0]
        elif vertex_placement == 'center':
            centers = triangle_geometric_features(triangles,positions,['barycenter'])[:,0]

        if exterior_vertex is not None:
            exterior_edge_centers = positions.values(exterior_edges).mean(axis=1)
            exterior_face_centers = face_centers.values(exterior_edge_triangle_ids)
            exterior_face_dual_centers = array_dict(centers,triangle_ids).values(exterior_edge_triangle_ids)

            exterior_edge_vectors = np.sign(np.einsum("...ij,...ij->...i",exterior_face_centers-exterior_edge_centers,exterior_face_dual_centers-exterior_edge_centers))[:,np.newaxis]*(exterior_edge_centers-exterior_face_dual_centers)
            print exterior_edge_vectors
            exterior_edge_vectors = exterior_edge_vectors/np.linalg.norm(exterior_edge_vectors,axis=1)[:,np.newaxis]
            exterior_dual_distance = np.linalg.norm(exterior_face_dual_centers-exterior_edge_centers,axis=1)

            #exterior_centers = face_centers.values(exterior_edge_triangle_ids) + exterior_coef*(positions.values(exterior_edges).mean(axis=1)-face_centers.values(exterior_edge_triangle_ids))
            exterior_centers = exterior_face_dual_centers + np.maximum(exterior_distance-exterior_dual_distance,0)[:,np.newaxis]*exterior_edge_vectors

        print list(triangle_ids),list(exterior_triangle_ids)
        print list(topomesh.wisps(2))

        dual_topomesh.update_wisp_property('barycenter',0,array_dict(list(centers)+list(exterior_centers),list(triangle_ids)+list(exterior_triangle_ids)))
    else:
        dual_topomesh.update_wisp_property('barycenter',0,face_centers)

    edges_to_remove = [e for e in dual_topomesh.wisps(1) if dual_topomesh.nb_borders(1,e) != 2]
    for e in edges_to_remove:
        dual_topomesh.remove_wisp(1,e)

    return clean_topomesh(dual_topomesh)


def triangulation_add_exterior(triangulation_topomesh):
    assert is_triangular(triangulation_topomesh)

    print np.min(list(triangulation_topomesh.wisps(0)))

    positions = triangulation_topomesh.wisp_property('barycenter',0)

    if not triangulation_topomesh.has_wisp(0,1):
        exterior_vertex = triangulation_topomesh.add_wisp(0,1)
    else:
        exterior_vertex = triangulation_topomesh.add_wisp(0,np.max(list(triangulation_topomesh.wisps(0)))+1)

    positions[exterior_vertex] = np.zeros(3)

    boundary_edges = [e for e in triangulation_topomesh.wisps(1) if triangulation_topomesh.nb_regions(1,e)==1]
    boundary_vertices = np.unique([list(triangulation_topomesh.borders(1,e)) for e in boundary_edges])

    exterior_edges = {}
    for v in boundary_vertices:
        e = triangulation_topomesh.add_wisp(1)
        triangulation_topomesh.link(1,e,v)
        triangulation_topomesh.link(1,e,exterior_vertex)
        exterior_edges[v] = e 

    for e in boundary_edges:
        t = triangulation_topomesh.add_wisp(2)
        for v in triangulation_topomesh.borders(1,e):
            triangulation_topomesh.link(2,t,exterior_edges[v])
        triangulation_topomesh.link(2,t,e)

    triangulation_topomesh.update_wisp_property('barycenter',0,positions)

    compute_topomesh_property(triangulation_topomesh,'triangles',0)
    compute_topomesh_property(triangulation_topomesh,'vertices',1)
    compute_topomesh_property(triangulation_topomesh,'regions',1)
    compute_topomesh_property(triangulation_topomesh,'triangles',1)
    compute_topomesh_property(triangulation_topomesh,'cells',1)
    compute_topomesh_property(triangulation_topomesh,'vertices',2)
    compute_topomesh_property(triangulation_topomesh,'cells',2)
    compute_topomesh_property(triangulation_topomesh,'regions',2)
    compute_topomesh_property(triangulation_topomesh,'vertices',3)
    compute_topomesh_property(triangulation_topomesh,'edges',3)
    compute_topomesh_property(triangulation_topomesh,'triangles',3)
    compute_topomesh_property(triangulation_topomesh,'vertices',3)
    compute_topomesh_property(triangulation_topomesh,'epidermis',2)

    return exterior_vertex


def triangulation_remove_exterior(triangulation_topomesh, exterior_vertex=1):
    assert is_triangular(triangulation_topomesh)
        
    triangles_to_remove = []
    for t in triangulation_topomesh.wisps(2):
        if exterior_vertex in triangulation_topomesh.borders(2,t,2):
            triangles_to_remove.append(t)
    edges_to_remove = []
    for e in triangulation_topomesh.wisps(1):
        if exterior_vertex in triangulation_topomesh.borders(1,e):
            edges_to_remove.append(e)
    triangulation_topomesh.remove_wisp(0,exterior_vertex)
    for e in edges_to_remove:
        triangulation_topomesh.remove_wisp(1,e)
    for t in triangles_to_remove:
        triangulation_topomesh.remove_wisp(2,t)

    triangulation_topomesh.update_wisp_property('barycenter',0,triangulation_topomesh.wisp_property('barycenter',0).values(list(triangulation_topomesh.wisps(0))),list(triangulation_topomesh.wisps(0)))   

    compute_topomesh_property(triangulation_topomesh,'triangles',0)
    compute_topomesh_property(triangulation_topomesh,'vertices',1)
    compute_topomesh_property(triangulation_topomesh,'regions',1)
    compute_topomesh_property(triangulation_topomesh,'triangles',1)
    compute_topomesh_property(triangulation_topomesh,'cells',1)
    compute_topomesh_property(triangulation_topomesh,'vertices',2)
    compute_topomesh_property(triangulation_topomesh,'cells',2)
    compute_topomesh_property(triangulation_topomesh,'regions',2)
    compute_topomesh_property(triangulation_topomesh,'vertices',3)
    compute_topomesh_property(triangulation_topomesh,'edges',3)
    compute_topomesh_property(triangulation_topomesh,'triangles',3)
    compute_topomesh_property(triangulation_topomesh,'vertices',3)
    compute_topomesh_property(triangulation_topomesh,'epidermis',2)


def topomesh_connected_components(topomesh,degree=2):
    if topomesh.nb_wisps(3)>0:
        new_component = np.max(list(topomesh.wisps(3)))+1
    else:
        new_component = 0
    considered_fids = set()
    component_cells = []
    for fid in topomesh.wisps(2):
        if not fid in considered_fids:
            c = new_component
            topomesh.add_wisp(3,c)
            topomesh.link(3,c,fid)
            considered_fids |= {fid}
            neighbor_fids = list(topomesh.border_neighbors(2,fid))
            while(len(neighbor_fids) > 0):
                n_fid = neighbor_fids.pop()
                if not n_fid in considered_fids:
                    topomesh.link(3,c,n_fid)
                    considered_fids |= {n_fid}
                    neighbor_fids += list(set(list(topomesh.border_neighbors(2,n_fid))).difference(considered_fids))
            component_cells += [c]
            new_component += 1

    component_meshes = [cell_topomesh(topomesh,cells=[c],copy_properties=True) for c in component_cells]
    component_size = [component.nb_wisps(2) for component in component_meshes]

    for c in component_cells:
        topomesh.remove_wisp(3,c)

    return [component_meshes[i] for i in np.argsort(component_size)][::-1]









