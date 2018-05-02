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
from scipy.cluster.vq import vq

from openalea.container import array_dict

try:
    import vtk

    def SetInput(obj, _input):
        if vtk.VTK_MAJOR_VERSION <= 5:
            obj.SetInput(_input)
        else:
            obj.SetInputData(_input)

except ImportError:
    print "VTK needs to be installed to use these functionalities"
    raise

from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import is_triangular, compute_topomesh_property
from openalea.cellcomplex.property_topomesh.property_topomesh_creation import triangle_topomesh

from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh
from openalea.cellcomplex.triangular_mesh import TriangularMesh

def property_topomesh_vtk_smoothing_decimation(input_topomesh, smoothing=2, decimation=2):
    if not is_triangular(input_topomesh):
        raise "The input mesh should be triangular!"

    # input_polydata = topomesh_to_triangular_mesh(input_topomesh,2)[0]._repr_vtk_()
    compute_topomesh_property(input_topomesh,'vertices',2)

    input_mesh = TriangularMesh()
    input_mesh.points = input_topomesh.wisp_property('barycenter',0).to_dict()
    input_mesh.triangles = input_topomesh.wisp_property('vertices',2).to_dict()

    input_polydata = input_mesh._repr_vtk_()    

    if smoothing>0:
        smoother = vtk.vtkWindowedSincPolyDataFilter()
        SetInput(smoother,input_polydata)
        smoother.BoundarySmoothingOn()
        smoother.FeatureEdgeSmoothingOn()
        smoother.SetFeatureAngle(120.0)
        smoother.SetPassBand(1)
        smoother.SetNumberOfIterations(smoothing)
        # smoother.NonManifoldSmoothingOn()
        # smoother.NormalizeCoordinatesOn()
        smoother.Update()

        print "Smoothing :",smoother.GetOutput().GetNumberOfPoints()

    if decimation>0:
        # decimate = vtk.vtkQuadricClustering()
        decimate = vtk.vtkQuadricDecimation()
        # decimate = vtk.vtkDecimatePro()

        if smoothing>0:
            SetInput(decimate,smoother.GetOutput())
        else:
            SetInput(decimate,input_polydata)
        decimate.SetTargetReduction(1-1./float(decimation))
        # decimate.SetNumberOfDivisions(int(decimation),int(decimation),int(decimation))
        # decimate.SetFeaturePointsAngle(60.0)
        decimate.Update()

        print "Decimation :",decimate.GetOutput().GetNumberOfPoints()

    if decimation>0:
        polydata = decimate.GetOutput()
    elif smoothing>0:
        polydata = smoother.GetOutput()
    else:
        polydata = input_polydata

    # input_polydata_points = np.array([input_polydata.GetPoints().GetPoint(p) for p in xrange(input_polydata.GetPoints().GetNumberOfPoints())])
    polydata_points = np.array([polydata.GetPoints().GetPoint(p) for p in xrange(polydata.GetPoints().GetNumberOfPoints())])

    polydata_triangles =  np.array([[polydata.GetCell(t).GetPointIds().GetId(i) for i in xrange(3)] for t in xrange(polydata.GetNumberOfCells())])

    if len(polydata_triangles)>0:
        return triangle_topomesh(polydata_triangles,array_dict(polydata_points,np.arange(len(polydata_points))))
    else:
        return None




