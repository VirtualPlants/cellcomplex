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

from openalea.cellcomplex.property_topomesh import PropertyTopomesh
from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property


def spatial_image_from_topomesh(topomesh,shape,offset=np.array([0.,0.,0.]),scale=1.0):
    """"todo"""

    start_time = time()
    print "--> Drawing topomesh image"

    if not topomesh.has_wisp_property('barycenter',degree=3,is_computed=True):
        compute_topomesh_property(topomesh,'barycenter',degree=3)
    if not topomesh.has_wisp_property('borders',degree=3,is_computed=True):
        compute_topomesh_property(topomesh,'borders',degree=3)
    if not topomesh.has_wisp_property('vertices',degree=2,is_computed=True):
        compute_topomesh_property(topomesh,'vertices',degree=2)

    if not topomesh.has_wisp_property('label',degree=3,is_computed=False):
        topomesh.add_wisp_property('label',degree=3)
    if not topomesh.has_wisp_property('label',degree=3,is_computed=True):
        topomesh.update_wisp_property('label',degree=3,values=np.array(list(topomesh.wisps(3)))+1)

    shape = tuple(np.array(shape)*scale)
    topomesh_img = np.ones(shape,np.uint16)

    # image_start_time = time()
    # print "  --> Generating coordinates ",shape
    # coords = np.mgrid[0:shape[0],0:shape[1],0:shape[2]]
    # coords = np.rollaxis(np.rollaxis(np.rollaxis(coords,3),3),3)/scale - offset
    # image_end_time = time()
    # print "  <-- Generating coordinates [",image_end_time-image_start_time,"s]"

    image_start_time = time()
    print "  --> Computing tetrahedra "
    cell_triangles = np.array(np.concatenate(topomesh.wisp_property('borders',degree=3).values()),int)
    cell_tetrahedra_cells = np.array(np.concatenate([[c for t in topomesh.wisp_property('borders',degree=3)[c]] for c in topomesh.wisps(3)]),int)
    cell_tetrahedra_labels = topomesh.wisp_property('label',degree=3).values(cell_tetrahedra_cells)
    cell_triangle_vertices = topomesh.wisp_property('vertices',degree=2).values(cell_triangles)
    cell_tetrahedra = np.concatenate([topomesh.wisp_property('barycenter',degree=3).values(cell_tetrahedra_cells)[:,np.newaxis], topomesh.wisp_property('barycenter',degree=0).values(cell_triangle_vertices)],axis=1)
    image_end_time = time()
    print "  <-- Computing tetrahedra (",cell_tetrahedra_cells.shape[0],") [",image_end_time-image_start_time,"s]"

    image_start_time = time()
    print "  --> Inverting matrices"
    cell_tetra_matrix = np.transpose(np.array([cell_tetrahedra[:,1],cell_tetrahedra[:,2],cell_tetrahedra[:,3]]) - cell_tetrahedra[:,0],axes=(1,2,0))
    cell_tetra_inv_matrix = np.linalg.inv(cell_tetra_matrix)
    image_end_time = time()
    print "  <-- Inverting matrices     [",image_end_time-image_start_time,"s]"

    # image_start_time = time()
    # print "  --> Detecting inner voxels"
    # xyz = (coords[:,:,:,np.newaxis] - cell_tetrahedra[:,0][np.newaxis,np.newaxis,np.newaxis,:])
    # lambdas = np.einsum('...ijk,...ik->...ij',cell_tetra_inv_matrix,xyz)
    # image_end_time = time()
    # print "  <-- Detecting inner voxels [",image_end_time-image_start_time,"s]"

    # image_start_time = time()
    # print "  --> Affecting image labels"
    # voxel_tetrahedra = ((lambdas[...,0]>=0)&(lambdas[...,1]>=0)&(lambdas[...,2]>=0)&(np.sum(lambdas,axis=4)<=1))
    # topomesh_img = np.array(np.max(voxel_tetrahedra*cell_tetrahedra_labels[np.newaxis,np.newaxis,np.newaxis],axis=3),np.uint16)
    # image_end_time = time()
    # print "  <-- Affecting image labels [",image_end_time-image_start_time,"s]"

    image_start_time = time()
    for t,tetra_label in enumerate(cell_tetrahedra_labels):
        tetra_vertices = cell_tetrahedra[t]

        coords = (np.mgrid[2.*scale*np.min(tetra_vertices[:,0]-1):2.*scale*np.max(tetra_vertices[:,0]+1),
                           2.*scale*np.min(tetra_vertices[:,1]-1):2.*scale*np.max(tetra_vertices[:,1]+1),
                           2.*scale*np.min(tetra_vertices[:,2]-1):2.*scale*np.max(tetra_vertices[:,2]+1)])/(2.*scale)
        coords = np.transpose(coords,(1,2,3,0))

        xyz = coords - tetra_vertices[0]
        lambdas = np.einsum('...jk,...k->...j',cell_tetra_inv_matrix[t],xyz)

        # voxel_tetrahedra = ((lambdas[...,0.]>=0)&(lambdas[...,1]>=0.)&(lambdas[...,2]>=0.)&(np.sum(lambdas,axis=3)<=1.))
        # tetra_where = np.maximum(np.minimum(np.array(scale*(coords+offset),int),np.array(shape,int)-1),0)
        # topomesh_img[tuple(np.transpose(tetra_where,axes=(3,0,1,2)))] = np.maximum(topomesh_img[tuple(np.transpose(tetra_where,axes=(3,0,1,2)))],np.array(voxel_tetrahedra*cell_tetrahedra_labels[t],np.uint16))

        tetra_where = np.array(scale*(coords[np.where((lambdas[...,0]>=0) & (lambdas[...,1]>=0) & (lambdas[...,2]>=0) & (np.sum(lambdas,axis=3)<=1))]+offset),int)
        tetra_where = tetra_where[np.where((tetra_where[...,0]>=0)&(tetra_where[...,1]>=0)&(tetra_where[...,2]>=0))]
        tetra_where = tetra_where[np.where((tetra_where[...,0]<shape[0])&(tetra_where[...,1]<shape[1])&(tetra_where[...,2]<shape[2]))]
        topomesh_img[tuple(np.transpose(tetra_where))] = tetra_label

        if (t%10000 == 0)&(t>0):
            image_end_time = time()
            print "  --> Drawing Tetrahedron ",t,"[",(image_end_time-image_start_time)/10000.,"s]"
            image_start_time = time()

    end_time = time()
    print "<-- Drawing topomesh image   [",end_time-start_time,"s]"

    return topomesh_img


def topomesh_line_rasterization(topomesh,shape):
    """ """
    start_time = time()
    print "--> Rasterizing topomesh"

    topomesh_img = np.ones(shape,np.uint16)

    if not topomesh.has_wisp_property('barycenter',degree=3,is_computed=True):
        compute_topomesh_property(topomesh,'barycenter',degree=3)
    if not topomesh.has_wisp_property('borders',degree=3,is_computed=True):
        compute_topomesh_property(topomesh,'borders',degree=3)
    if not topomesh.has_wisp_property('vertices',degree=2,is_computed=True):
        compute_topomesh_property(topomesh,'vertices',degree=2)

    if not topomesh.has_wisp_property('label',degree=3,is_computed=False):
        topomesh.add_wisp_property('label',degree=3)
    if not topomesh.has_wisp_property('label',degree=3,is_computed=True):
        topomesh.update_wisp_property('label',degree=3,values=np.array(list(topomesh.wisps(3)))+1)

    for c in topomesh.wisps(3):
        cell_start_time = time()
        print "  --> Rasterizing cell",c

        cell_img = np.zeros(shape,np.uint16)

        cell_triangles = topomesh.wisp_property('borders',degree=3)[c]
        cell_triangle_vertices = topomesh.wisp_property('vertices',degree=2).values(cell_triangles)
        cell_faces = topomesh.wisp_property('barycenter',degree=0).values(cell_triangle_vertices)

        coords = np.mgrid[np.min(cell_faces[:,0]-1):np.max(cell_faces[:,0]+1),
                          np.min(cell_faces[:,1]-1):np.max(cell_faces[:,1]+1)]
        coords = np.transpose(coords,(1,2,0))
        coords = np.concatenate([coords,np.zeros_like(coords[...,0])[...,np.newaxis]-10],axis=2)
        #coords = np.concatenate([coords,(np.min(cell_faces[:,2]-1))*np.ones_like(coords[...,0])[...,np.newaxis]],axis=2)

        cell_face_edge1 = cell_faces[:,1] - cell_faces[:,0]
        cell_face_edge2 = cell_faces[:,2] - cell_faces[:,0]
        cell_rays_t     = coords[:,:,np.newaxis] - cell_faces[:,0][np.newaxis,np.newaxis,:]
        cell_rays_d     = np.zeros_like(cell_rays_t)
        cell_rays_d[...,2] = 1.

        cell_face_p = np.cross(cell_rays_d,cell_face_edge2[np.newaxis,np.newaxis,:])
        cell_face_q = np.cross(cell_rays_t,cell_face_edge1[np.newaxis,np.newaxis,:])

        cell_face_norm = np.einsum('...ij,...ij->...i',cell_face_p,cell_face_edge1[np.newaxis,np.newaxis,:])
        cell_face_distance = np.einsum('...ij,...ij->...i',cell_face_q,cell_face_edge2[np.newaxis,np.newaxis,:])
        cell_face_projection_u = np.einsum('...ij,...ij->...i',cell_face_p,cell_rays_t)
        cell_face_projection_v = np.einsum('...ij,...ij->...i',cell_face_q,cell_rays_d)
        
        cell_face_ray_intersection = np.concatenate([cell_face_distance[...,np.newaxis],cell_face_projection_u[...,np.newaxis],cell_face_projection_v[...,np.newaxis]],axis=3)/cell_face_norm[...,np.newaxis]        

        for z in xrange(np.min(cell_faces[:,2]-1),np.max(cell_faces[:,2]+2)):
            rays_intersections = np.where((cell_face_ray_intersection[...,0]>z)&(cell_face_ray_intersection[...,1]>0)&(cell_face_ray_intersection[...,1]<1)&(cell_face_ray_intersection[...,2]>0)&(cell_face_ray_intersection[...,2]<1))
            layer = np.zeros_like(coords[:,:,0])

        cell_end_time = time()
        print "  <-- Rasterizing cell",c,"     [",end_time-start_time,"s]"


    end_time = time()
    print "<-- Rasterizing topomesh     [",end_time-start_time,"s]"

    