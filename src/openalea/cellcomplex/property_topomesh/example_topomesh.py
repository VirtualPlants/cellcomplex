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
from openalea.container import array_dict

from openalea.cellcomplex.property_topomesh.property_topomesh_creation import triangle_topomesh, quad_topomesh, poly_topomesh
from openalea.cellcomplex.property_topomesh.utils.delaunay_tools import delaunay_triangulation
from openalea.cellcomplex.property_topomesh.property_topomesh_creation import dual_topomesh

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


def triangular_grid_topomesh(size=1, resolution=1., center=np.zeros(3)):
    """
    """

    return topomesh


def square_grid_topomesh(size=1, resolution=1., center=np.zeros(3)):
    """
    """

    x = y = np.linspace(-size,size,2*size+1)*resolution
    x += center[0]
    y += center[1]

    xx,yy = np.meshgrid(x,y)

    zz = np.zeros_like(xx) + center[2]

    square_bl = np.concatenate([np.arange(2*size+1)[:-1]+i*(2*size+1) for i in np.arange(2*size+1)[:-1]])
    square_br = np.concatenate([np.arange(2*size+1)[1:]+i*(2*size+1) for i in np.arange(2*size+1)[:-1]])
    square_tl = np.concatenate([np.arange(2*size+1)[:-1]+i*(2*size+1) for i in np.arange(2*size+1)[1:]])
    square_tr = np.concatenate([np.arange(2*size+1)[1:]+i*(2*size+1) for i in np.arange(2*size+1)[1:]])

    squares = np.transpose([square_bl,square_tl,square_tr,square_br])

    positions = dict(zip(range(len(np.ravel(xx))),np.transpose([np.ravel(xx),np.ravel(yy),np.ravel(zz)])))

    topomesh = quad_topomesh(squares,positions,faces_as_cells=True)

    return topomesh


def hexagonal_grid_topomesh(size=1, resolution=1., center=np.zeros(3)):
    """
    """

    x = []
    y = []
    p = []
    for n in range(size+1,2*size+2):
        for r in [0,1]:
            p += list(len(p)+np.arange(n+r))
            x += list((np.arange(n+r)-((n+r-1)/2.))*resolution*np.sqrt(3.))
            row_y = (((2*size+1 - n)+1)/2)*3
            row_side = (((2*size+1 - n)+1)%2)
            y += list((row_y+(2*row_side-1)-(row_side==r)*(2*row_side-1)/2.)*np.ones(n+r)*resolution)


    for n in range(size+1,2*size+2)[::-1]:
        for r in [1,0]:
            p += list(len(p)+np.arange(n+r))
            x += list((np.arange(n+r)-((n+r-1)/2.))*resolution*np.sqrt(3.))
            row_y = (((2*size+1 - n)+1)/2)*3
            row_side = (((2*size+1 - n)+1)%2)
            y += list(-(row_y+(2*row_side-1)-(row_side==r)*(2*row_side-1)/2.)*np.ones(n+r)*resolution)

    x = np.array(x) + center[0]
    y = np.array(y) + center[1]
    z = np.zeros_like(x) + center[2]

    positions = dict(zip(range(len(x)),np.transpose([x,y,z])))

    index_list = []
    row_length_list = []
    row_gaps = []
    start = 0
    for n in range(size+1,2*size+1):
        index_list += list(np.arange(n)+start)
        start += n+n+1
        row_length_list += list(n*np.ones(n,int))
        for h in xrange(n):
            row_gaps += [[n+1,n+1,n+1,-(n+2),-(n+1)]]
    for n in [2*size+1]:
        index_list += list(np.arange(n)+start)
        start += n+n+2
        row_length_list += list(n*np.ones(n,int))
        for h in xrange(n):
            row_gaps += [[n+1,n+1,n,-(n+1),-(n+1)]]
    for n in range(2*size,size,-1):
        index_list += list(np.arange(n)+start)
        start += n+n+3
        row_length_list += list(n*np.ones(n,int))
        for h in xrange(n):
            row_gaps += [[n+2,n+1,n,-(n+1),-(n+1)]]

    print index_list
    print row_length_list

    print "Hexagons"

    hexagons = []
    for i,n,gaps in zip(index_list,row_length_list,row_gaps):
        hexagon = [i]
        for p, gap in zip(xrange(5),gaps):
            hexagon += [hexagon[-1]+gap]
        print hexagon
        hexagons += [hexagon]

    topomesh = poly_topomesh(hexagons,positions,faces_as_cells=True)

    return topomesh

def circle_voronoi_topomesh(size = 1,resolution = 1.,circle_size = 100.,z_coef = 0.):
    n_cells = 3*size*(size-1)+1
    radius = size*resolution

    circle_thetas = np.linspace(-np.pi,np.pi-2*np.pi/float(circle_size),circle_size)
    circle_points = np.transpose([radius*np.cos(circle_thetas),radius*np.sin(circle_thetas)])

    cell_thetas = np.array([np.pi*np.random.randint(-180,180)/180. for c in xrange(n_cells)])
    cell_distances = 0.5*radius*np.sqrt([np.random.rand() for c in xrange(n_cells)])

    cell_points = np.transpose([cell_distances*np.cos(cell_thetas),cell_distances*np.sin(cell_thetas)])

    omega_forces = dict(repulsion=0.5)
    sigma_deformation = 2.*radius/float(n_cells)

    for iteration in xrange(n_cells/2):
        cell_to_cell_vectors = np.array([[p-q for q in cell_points] for p in cell_points])
        cell_to_cell_distances = np.linalg.norm(cell_to_cell_vectors,axis=2)/radius
        cell_to_circle_vectors = np.array([[p-q for q in circle_points] for p in cell_points])
        cell_to_circle_distances = np.linalg.norm(cell_to_circle_vectors,axis=2)/radius
        
        deformation_force = np.zeros_like(cell_points)
        
        cell_repulsion_force = np.nansum(cell_to_cell_vectors/np.power(cell_to_cell_distances,3)[:,:,np.newaxis],axis=1)
        circle_repulsion_force = np.nansum(cell_to_circle_vectors/np.power(cell_to_circle_distances,3)[:,:,np.newaxis],axis=1)

        deformation_force += omega_forces['repulsion']*cell_repulsion_force
        deformation_force += 1.5*(n_cells/float(circle_size))*omega_forces['repulsion']*circle_repulsion_force

        deformation_force_amplitude = np.linalg.norm(deformation_force,axis=1)
        deformation_force = np.minimum(1.0,sigma_deformation/deformation_force_amplitude)[:,np.newaxis] * deformation_force

        cell_points += deformation_force
        cell_points = np.minimum(1.0,radius/(np.linalg.norm(cell_points,axis=1)))[:,np.newaxis] * cell_points

    all_positions = array_dict(np.transpose([np.concatenate([cell_points[:,0],circle_points[:,0]]),np.concatenate([cell_points[:,1],circle_points[:,1]]),np.zeros(n_cells+circle_size)]),keys=np.arange(n_cells+circle_size).astype(int))
    triangles = all_positions.keys()[delaunay_triangulation(all_positions.values())]

    radial_distances = np.linalg.norm(all_positions.values(),axis=1)
    radial_z = z_coef*np.power(radial_distances/radius,2)
    all_positions = array_dict(np.transpose([all_positions.values()[:,0],all_positions.values()[:,1],radial_z]),keys=all_positions.keys())

    triangulation_topomesh = triangle_topomesh(triangles,all_positions)
    cell_topomesh = dual_topomesh(triangulation_topomesh,2,vertex_positions='voronoi')

    return cell_topomesh


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


