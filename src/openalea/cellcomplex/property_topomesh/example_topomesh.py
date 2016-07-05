# -*- coding: utf-8 -*-
# -*- python -*-
#
#       PropertyTopomesh
#
#       Copyright 2016 INRIA - CIRAD - INRA
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

from openalea.mesh.property_topomesh_creation import triangle_topomesh

def square_topomesh(side_length = 1):
    points = {}
    points[0] = [0,0,0]
    points[1] = [side_length,0,0]
    points[2] = [0,side_length,0]
    points[3] = [side_length,side_length,0]

    triangles = [[0,1,2],[3,2,1]]

    return triangle_topomesh(triangles, points)


def hexagon_topomesh(side_length = 1):
    points = {}
    points[0] = [0,0,0]
    for p in xrange(6):
        points[p+1] = [side_length*np.cos(p*np.pi/3.),side_length*np.sin(p*np.pi/3.),0]

    triangles = []
    for p in xrange(6):
        triangles += [[0,p + 1,(p+1)%6 + 1]]

    return triangle_topomesh(triangles, points)

def icosahedron_topomesh(size=1.0):
    ico_points = {}
    ico_points[0] = size*np.array([ 0.      , -1.      ,  0.      ])
    ico_points[1] = size*np.array([ 0.7236  , -0.447215,  0.52572 ])
    ico_points[2] = size*np.array([-0.276385, -0.447215,  0.85064 ])
    ico_points[3] = size*np.array([-0.894425, -0.447215,  0.      ])
    ico_points[4] = size*np.array([-0.276385, -0.447215, -0.85064 ])
    ico_points[5] = size*np.array([ 0.7236  , -0.447215, -0.52572 ])
    ico_points[6] = size*np.array([ 0.276385,  0.447215,  0.85064 ])
    ico_points[7] = size*np.array([-0.7236  ,  0.447215,  0.52572 ])
    ico_points[8] = size*np.array([-0.7236  ,  0.447215, -0.52572 ])
    ico_points[9] = size*np.array([ 0.276385,  0.447215, -0.85064 ])
    ico_points[10] = size*np.array([ 0.894425,  0.447215,  0.      ])
    ico_points[11] = size*np.array([ 0.      ,  1.      ,  0.      ])

    ico_triangles = [[ 0,  1,  2],
        [ 0,  1,  5],
        [ 0,  2,  3],
        [ 0,  3,  4],
        [ 0,  4,  5],
        [ 1,  5, 10],
        [ 1,  2,  6],
        [ 2,  3,  7],
        [ 3,  4,  8],
        [ 4,  5,  9],
        [ 1,  6, 10],
        [ 2,  6,  7],
        [ 3,  7,  8],
        [ 4,  8,  9],
        [ 5,  9, 10],
        [ 6, 10, 11],
        [ 6,  7, 11],
        [ 7,  8, 11],
        [ 8,  9, 11],
        [ 9, 10, 11]]

    return triangle_topomesh(ico_triangles, ico_points)


def sphere_topomesh(radius=1.0,center=np.zeros(3)):
    from openalea.container import array_dict
    from openalea.cellcomplex.property_topomesh.property_topomesh_optimization import topomesh_triangle_split

    ico = icosahedron_topomesh()
    topomesh = topomesh_triangle_split(topomesh_triangle_split(ico))

    positions = topomesh.wisp_property('barycenter',0)
    new_positions = array_dict(center + radius*positions.values()/np.linalg.norm(positions.values(),axis=1)[:,np.newaxis],positions.keys())
    topomesh.update_wisp_property('barycenter',0,new_positions)

    return topomesh


def vtk_ellipsoid_topomesh(ellipsoid_radius=50.0, ellipsoid_scales=[1,1,1], ellipsoid_axes=np.diag(np.ones(3)), ellipsoid_center=np.zeros(3)):
    """
    """      
    import vtk
    from openalea.cellcomplex.property_topomesh.utils.image_tools import vtk_polydata_to_triangular_mesh

    ico = vtk.vtkPlatonicSolidSource()
    ico.SetSolidTypeToIcosahedron()
    ico.Update()

    subdivide = vtk.vtkLoopSubdivisionFilter()
    subdivide.SetNumberOfSubdivisions(3)
    subdivide.SetInputConnection(ico.GetOutputPort())
    subdivide.Update()

    scale_transform = vtk.vtkTransform()
    scale_factor = ellipsoid_radius/(np.sqrt(2)/2.)
    scale_transform.Scale(scale_factor,scale_factor,scale_factor)

    ellipsoid_sphere = vtk.vtkTransformPolyDataFilter()
    ellipsoid_sphere.SetInput(subdivide.GetOutput())
    ellipsoid_sphere.SetTransform(scale_transform)
    ellipsoid_sphere.Update()

    ellipsoid_transform = vtk.vtkTransform()
    axes_transform = vtk.vtkLandmarkTransform()
    source_points = vtk.vtkPoints()
    source_points.InsertNextPoint([1,0,0])
    source_points.InsertNextPoint([0,1,0])
    source_points.InsertNextPoint([0,0,1])
    target_points = vtk.vtkPoints()
    target_points.InsertNextPoint(ellipsoid_axes[0])
    target_points.InsertNextPoint(ellipsoid_axes[1])
    target_points.InsertNextPoint(ellipsoid_axes[2])
    axes_transform.SetSourceLandmarks(source_points)
    axes_transform.SetTargetLandmarks(target_points)
    axes_transform.SetModeToRigidBody()
    axes_transform.Update()
    ellipsoid_transform.SetMatrix(axes_transform.GetMatrix())
    ellipsoid_transform.Scale(ellipsoid_scales[0],ellipsoid_scales[1],ellipsoid_scales[2])
    center_transform = vtk.vtkTransform()
    center_transform.Translate(ellipsoid_center[0],ellipsoid_center[1],ellipsoid_center[2])
    center_transform.Concatenate(ellipsoid_transform)
    ellipsoid_ellipsoid = vtk.vtkTransformPolyDataFilter()
    ellipsoid_ellipsoid.SetInput(ellipsoid_sphere.GetOutput())
    ellipsoid_ellipsoid.SetTransform(center_transform)
    ellipsoid_ellipsoid.Update()

    ellipsoid_mesh = vtk_polydata_to_triangular_mesh(ellipsoid_ellipsoid.GetOutput())

    topomesh = triangle_topomesh(ellipsoid_mesh.triangles.values(),ellipsoid_mesh.points)

    return topomesh


