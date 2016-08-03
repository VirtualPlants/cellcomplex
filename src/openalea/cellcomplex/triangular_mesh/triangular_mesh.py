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

def isiterable(obj):
    try:
        iter(obj)
        return True
    except:
        return False

class TriangularMesh(object):
    def __init__(self):
        self.points = {}
        self.point_data = {}
        self.edges = {}
        self.edge_data = {}
        self.triangles = {}
        self.triangle_data = {}

        self.point_radius = 2.0
        self.char_dimension = None

    def _repr_geom_(self):
        import openalea.plantgl.all as pgl
        import vplants.plantgl.ext.color as color
        from openalea.container import array_dict
        
        scene = pgl.Scene()
        colormap = color.GlasbeyMap(0,256)

        if len(self.triangles) > 0:
            triangle_points = array_dict(self.points).values(self.triangles.values())
            triangle_normals = np.cross(triangle_points[:,1]-triangle_points[:,0],triangle_points[:,2]-triangle_points[:,0])
            mesh_center = np.mean(self.points.values(),axis=0)
            reversed_normals = np.array(self.triangles.keys())[np.where(np.einsum('ij,ij->i',triangle_normals,triangle_points[:,0]-mesh_center) < 0)[0]]
            for t in reversed_normals:
                self.triangles[t] = list(reversed(self.triangles[t]))
                                         

            points_index = array_dict(np.arange(len(self.points)),self.points.keys())

            if isiterable(self.triangle_data.values()[0]):
                colors = [pgl.Color4(colormap(self.triangle_data[t][0]%256).i3tuple()) if self.triangle_data.has_key(t) else pgl.Color4(colormap(0).i3tuple()) for t in self.triangles.keys()]
            else:
                colors = [pgl.Color4(colormap(self.triangle_data[t]%256).i3tuple()) if self.triangle_data.has_key(t) else pgl.Color4(colormap(0).i3tuple()) for t in self.triangles.keys()]

            # scene += pgl.Shape(pgl.FaceSet(self.points.values(),list(points_index.values(self.triangles.values()))),pgl.Material((255,255,255)))
            scene += pgl.Shape(pgl.FaceSet(self.points.values(),list(points_index.values(self.triangles.values())),colorList=colors,colorPerVertex=False))
            #for t in self.triangles.keys():
            #    scene += pgl.Shape(pgl.FaceSet([self.points[p] for p in self.triangles[t]],[list(range(3))]),pgl.Material(colormap(self.triangle_data[t]%256).i3tuple()),id=t)
        else:
            for p in self.points.keys():
                mat = pgl.Material(colormap(p%256).i3tuple(),transparency=0.0,name=p)
                scene += pgl.Shape(pgl.Translated(self.points[p],pgl.Sphere(self.point_radius,slices=16,stacks=16)),mat,id=p)
        return scene
    
    def _repr_vtk_(self):
        import vtk
        from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        from openalea.container import array_dict

        if len(self.triangles) > 0:
            vtk_mesh = vtk.vtkPolyData()
            
            vtk_points = vtk.vtkPoints()
            vtk_point_data = vtk.vtkDoubleArray()
            
            vtk_triangles = vtk.vtkCellArray()
            vtk_triangle_data = vtk.vtkDoubleArray()
            
            if len(self.triangle_data)>0 and np.array(self.triangle_data.values()).ndim > 1:
                if np.array(self.triangle_data.values()).ndim==2:
                    vtk_point_data.SetNumberOfComponents(np.array(self.triangle_data.values()).shape[1])
                elif np.array(self.triangle_data.values()).ndim==3:
                    vtk_point_data.SetNumberOfComponents(np.array(self.triangle_data.values()).shape[1]*np.array(self.triangle_data.values()).shape[2])
                
                mesh_points = []
                positions = array_dict(self.points)
                for t in self.triangles.keys():
                    triangle_center = positions.values(self.triangles[t]).mean(axis=0)
                    tid = vtk_points.InsertNextPoint(triangle_center)
                    mesh_points.append(tid)

                    if self.triangle_data.has_key(t):
                        if isiterable(self.triangle_data[t]):
                            if np.array(self.triangle_data[t]).ndim==1:
                                vtk_point_data.InsertTuple(tid,self.triangle_data[t])
                            else:
                                vtk_point_data.InsertTuple(tid,np.concatenate(self.triangle_data[t]))
                mesh_points = array_dict(mesh_points,self.triangles.keys())
                
                vtk_mesh.SetPoints(vtk_points)

                if np.array(self.triangle_data.values()).ndim==2:
                    vtk_mesh.GetPointData().SetVectors(vtk_point_data)
                elif np.array(self.triangle_data.values()).ndim==3:
                    vtk_mesh.GetPointData().SetTensors(vtk_point_data)

            else:
                from time import time
                start_time = time()

                double_array = vtk.vtkDoubleArray()
                double_array = numpy_to_vtk(np.array(self.points.values()), deep=True, array_type=vtk.VTK_DOUBLE)
                vtk_points.SetData(double_array)

                if self.point_data.has_key(self.points.keys()[0]):
                    if isiterable(self.point_data[self.points.keys()[0]]):
                        vtk_point_data = numpy_to_vtk(np.array([v[0] for v in self.point_data.values()]), deep=True, array_type=vtk.VTK_DOUBLE)
                    else:
                        vtk_point_data = numpy_to_vtk(np.array(self.point_data.values()), deep=True, array_type=vtk.VTK_DOUBLE)
                vtk_point_data.SetNumberOfComponents(1)
                
                # mesh_points = []

                # for p in self.points.keys():
                    # pid = vtk_points.InsertNextPoint(self.points[p])
                    # mesh_points.append(pid)
                    # if self.point_data.has_key(p):
                        # if isiterable(self.point_data[p]):
                            # vtk_point_data.InsertValue(pid,self.point_data[p][0])
                        # else:
                            # vtk_point_data.InsertValue(pid,self.point_data[p])

                mesh_points = array_dict(np.arange(len(self.points)),self.points.keys())
                # mesh_points = array_dict(mesh_points,self.points.keys())
                
                if len(self.point_data) > 0:
                    vtk_mesh.GetPointData().SetScalars(vtk_point_data)

                triangles_vtk = np.concatenate([3*np.ones(len(self.triangles),int)[:,np.newaxis],mesh_points.values(np.array(self.triangles.values()))],axis=1)
                vtk_ids = numpy_to_vtkIdTypeArray(triangles_vtk, deep=True)
                vtk_triangles.SetCells(len(self.triangles), vtk_ids)

                if self.triangle_data.has_key(self.triangles.keys()[0]):
                    if isiterable(self.triangle_data[self.triangles.keys()[0]]):
                        if np.array(self.triangle_data[self.triangles.keys()[0]]).ndim==1:
                            vtk_triangle_data = numpy_to_vtk(np.array(self.triangle_data.values()), deep=True, array_type=vtk.VTK_DOUBLE)
                        else:
                            vtk_triangle_data = numpy_to_vtk(np.array([np.concatenate(v) for v in self.triangle_data.values()]), deep=True, array_type=vtk.VTK_DOUBLE)
                    else:
                        vtk_triangle_data = numpy_to_vtk(np.array(self.triangle_data.values()), deep=True, array_type=vtk.VTK_DOUBLE)

                # for t in self.triangles.keys():
                #     poly = vtk_triangles.InsertNextCell(3)
                #     for i in xrange(3):
                #         vtk_triangles.InsertCellPoint(mesh_points[self.triangles[t][i]])
                    # if self.triangle_data.has_key(t):
                    #     if isiterable(self.triangle_data[t]):
                    #         if np.array(self.triangle_data[t]).ndim==1:
                    #             vtk_triangle_data.InsertTuple(poly,self.triangle_data[t])
                    #         else:
                    #             vtk_triangle_data.InsertTuple(poly,np.concatenate(self.triangle_data[t]))
                    #         # vtk_triangle_data.InsertValue(poly,self.triangle_data[t][0])
                    #     else:
                    #         vtk_triangle_data.InsertValue(poly,self.triangle_data[t])
                
                vtk_mesh.SetPoints(vtk_points)
                vtk_mesh.SetPolys(vtk_triangles)

                if len(self.triangle_data) > 0:
                    vtk_mesh.GetCellData().SetScalars(vtk_triangle_data)

                end_time = time()
                print "--> Converting to VTK PolyData     [",end_time-start_time,"s]"

            return vtk_mesh

        elif len(self.edges) > 0:
            vtk_mesh = vtk.vtkPolyData()
            vtk_points = vtk.vtkPoints()
            vtk_point_data = vtk.vtkDoubleArray()
            vtk_lines = vtk.vtkCellArray()
            vtk_line_data = vtk.vtkDoubleArray()

            mesh_points = []
            for p in self.points.keys():
                pid = vtk_points.InsertNextPoint(self.points[p])
                mesh_points.append(pid)
                if self.point_data.has_key(p):
                    if isiterable(self.point_data[p]):
                        vtk_point_data.InsertValue(pid,self.point_data[p][0])
                    else:
                        vtk_point_data.InsertValue(pid,self.point_data[p])
            mesh_points = array_dict(mesh_points,self.points.keys())
            if len(self.point_data) > 0:
                vtk_mesh.GetPointData().SetScalars(vtk_point_data)

            for e in self.edges.keys():
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0,mesh_points[self.edges[e][0]])
                line.GetPointIds().SetId(1,mesh_points[self.edges[e][1]])
                edge = vtk_lines.InsertNextCell(line)
                if self.edge_data.has_key(e):
                    if isiterable(self.edge_data[e]):
                        vtk_line_data.InsertValue(edge,self.edge_data[e][0])
                    else:
                        vtk_line_data.InsertValue(edge,self.edge_data[e])
                else:
                    vtk_line_data.InsertValue(edge,0)

            vtk_mesh.SetPoints(vtk_points)
            vtk_mesh.SetLines(vtk_lines)
            vtk_mesh.GetCellData().SetScalars(vtk_line_data)

            return vtk_mesh

        else:
            vtk_mesh = vtk.vtkPolyData()
            vtk_points = vtk.vtkPoints()
            
            vtk_cells = vtk.vtkDoubleArray()
            if len(self.point_data)>0 and np.array(self.point_data.values()).ndim==2:
                vtk_cells.SetNumberOfComponents(np.array(self.point_data.values()).shape[1])
            elif len(self.point_data)>0 and np.array(self.point_data.values()).ndim==3:
                vtk_cells.SetNumberOfComponents(np.array(self.point_data.values()).shape[1]*np.array(self.point_data.values()).shape[2])
            for p in self.points.keys():
                pid = vtk_points.InsertNextPoint(self.points[p])
                if self.point_data.has_key(p):
                    if isiterable(self.point_data[p]):
                        if np.array(self.point_data[p]).ndim==1:
                            cell = vtk_cells.InsertNextTuple(self.point_data[p])
                        else:
                            cell = vtk_cells.InsertNextTuple(np.concatenate(self.point_data[p]))
                    else:
                        vtk_cells.InsertValue(pid,self.point_data[p])
                else:
                    vtk_cells.InsertValue(pid,p)
            vtk_mesh.SetPoints(vtk_points)
            if len(self.point_data)>0 and np.array(self.point_data.values()).ndim==2:
                vtk_mesh.GetPointData().SetVectors(vtk_cells)
            elif len(self.point_data)>0 and np.array(self.point_data.values()).ndim==3:
                vtk_mesh.GetPointData().SetTensors(vtk_cells)
            else:
                vtk_mesh.GetPointData().SetScalars(vtk_cells)


            return vtk_mesh

    def data(self):
        if len(self.triangle_data) > 0:
            data = np.array(self.triangle_data.values())
        elif len(self.point_data) > 0:
            data = np.array(self.point_data.values())
        elif len(self.triangles) > 0:
            data = np.array(self.triangles.keys())
        else:
            data = np.array(self.points.keys())
        return data

    def min(self):
        if self.data().ndim == 1:
            return np.nanmin(self.data())
        elif self.data().ndim == 2:
            return np.nanmin(np.linalg.norm(self.data(),axis=1))
        elif self.data().ndim == 3:
            return np.nanmin(np.sqrt(np.trace(np.power(self.data(),2),axis1=1,axis2=2)))

    def max(self):
        if self.data().ndim == 1:
            return np.nanmax(self.data())
        elif self.data().ndim == 2:
            return np.nanmax(np.linalg.norm(self.data(),axis=1))
        elif self.data().ndim == 3:
            return np.nanmax(np.sqrt(np.trace(np.power(self.data(),2),axis1=1,axis2=2)))

    def mean(self):
        if self.data().ndim == 1:
            return np.nanmean(self.data())
        elif self.data().ndim == 2:
            return np.nanmean(np.linalg.norm(self.data(),axis=1))
        elif self.data().ndim == 3:
            return np.nanmean(np.sqrt(np.trace(np.power(self.data(),2),axis1=1,axis2=2)))

    def bounding_box(self):
        if len(self.points)>0:
            extent_min = (np.min(self.points.values(),axis=0))
            extent_max = (np.max(self.points.values(),axis=0))
            return zip(extent_min,extent_max)
        else:
            return zip([0,0,0],[0,0,0])

    def characteristic_dimension(self):
        if self.char_dimension is None:
            if len(self.points)>1:
                if len(self.triangles)>0:
                    triangle_edge_list = [[1,2],[0,2],[0,1]]
                    triangle_edges = np.concatenate(np.array(self.triangles.values())[:,triangle_edge_list])
                    triangle_edge_points = array_dict(self.points).values(triangle_edges)
                    triangle_edge_vectors = triangle_edge_points[:,1] - triangle_edge_points[:,0]
                    triangle_edge_lengths = np.linalg.norm(triangle_edge_vectors,axis=1)
                    self.char_dimension = triangle_edge_lengths.mean()
                elif len(self.edges)>0:
                    edges = np.array(self.edges.values())
                    edge_points = array_dict(self.points).values(edges)
                    edge_vectors = edge_points[:,1] - edge_points[:,0]
                    edge_lengths = np.linalg.norm(edge_vectors,axis=1)
                    self.char_dimension = edge_lengths.mean()
                else:
                    #from scipy.cluster.vq import vq
                    #point_distances = np.sort([vq(np.array(self.points.values()),np.array([self.points[p]]))[1] for p in self.points.keys()])
                    # self.char_dimension = point_distances[:,1].mean()
                    bbox = np.array(self.bounding_box())
                    bbox_volume = np.prod(bbox[:,1] - bbox[:,0])
                    point_volume = bbox_volume/float(2.*len(self.points))
                    self.char_dimension = np.power(3.*point_volume/(4.*np.pi),1/3.)

                return self.char_dimension
            else:
                return 1.
        else:
            return self.char_dimension



def point_triangular_mesh(point_positions, point_data=None):
    points_mesh = TriangularMesh()
    points_mesh.points = point_positions
    if point_data is not None:
        points_mesh.point_data = point_data
    return points_mesh


def save_ply_triangular_mesh(mesh,ply_filename,intensity_range=None):
    from time import time

    start_time =time()
    print "--> Saving .ply"

    ply_file = open(ply_filename,'w+')

    ply_file.write("ply\n")
    ply_file.write("format ascii 1.0\n")
    ply_file.write("element vertex "+str(len(mesh.points))+"\n")
    ply_file.write("property float x\n")
    ply_file.write("property float y\n")
    ply_file.write("property float z\n")
    if len(mesh.point_data)>0 and np.ndim(mesh.point_data.values()[0])==0:
        ply_file.write("property uchar red\n")
        ply_file.write("property uchar green\n")
        ply_file.write("property uchar blue\n")
        if intensity_range is None:
            point_data = array_dict(np.array(mesh.point_data.values(),int)%256,mesh.point_data.keys())
        else:
            point_data = array_dict((np.array(mesh.point_data.values())-intensity_range[0])/(intensity_range[1]-intensity_range[0]),mesh.point_data.keys())

    ply_file.write("element face "+str(len(mesh.triangles))+"\n")
    ply_file.write("property list uchar int vertex_indices\n")
    if len(mesh.triangle_data)>0 and np.ndim(mesh.triangle_data.values()[0])==0:
        ply_file.write("property uchar red\n")
        ply_file.write("property uchar green\n")
        ply_file.write("property uchar blue\n")
        if intensity_range is None:
            triangle_data = array_dict(np.array(mesh.triangle_data.values(),int)%256,mesh.triangle_data.keys())
        else:
            triangle_data = array_dict((np.array(mesh.triangle_data.values())-intensity_range[0])/(intensity_range[1]-intensity_range[0]),mesh.triangle_data.keys())
    ply_file.write("end_header\n")

    vertex_index = {}
    for v,p in enumerate(mesh.points.keys()):
        ply_file.write(str(mesh.points[p][0])+" ")
        ply_file.write(str(mesh.points[p][1])+" ")
        ply_file.write(str(mesh.points[p][2]))
        if len(mesh.point_data)>0 and np.ndim(mesh.point_data.values()[0])==0:
            ply_file.write(" "+str(point_data[p]))
            ply_file.write(" "+str(point_data[p]))
            ply_file.write(" "+str(point_data[p]))
        ply_file.write("\n")
        vertex_index[p] = v

    for t in mesh.triangles.keys():
        ply_file.write("3 ")
        ply_file.write(str(vertex_index[mesh.triangles[t][0]])+" ")
        ply_file.write(str(vertex_index[mesh.triangles[t][1]])+" ")
        ply_file.write(str(vertex_index[mesh.triangles[t][2]]))
        if len(mesh.triangle_data)>0 and np.ndim(mesh.triangle_data.values()[0])==0:
            ply_file.write(" "+str(triangle_data[t]))
            ply_file.write(" "+str(triangle_data[t]))
            ply_file.write(" "+str(triangle_data[t]))
        ply_file.write("\n")

    ply_file.flush()
    ply_file.close()

    end_time = time()
    print "<-- Saving .ply        [",end_time-start_time,"s]"



def save_ply_triangle_mesh(ply_filename, positions, triangles={}, edges={}, vertex_properties={}, triangle_properties={}, edge_properties={}):
    """
    """

    from time import time

    if isinstance(positions,list) or isinstance(positions,np.ndarray):
        positions = dict(zip(range(len(positions)),positions))
    if isinstance(triangles,list) or isinstance(triangles,np.ndarray):
        triangles = dict(zip(range(len(triangles)),triangles))
    if isinstance(edges,list) or isinstance(edges,np.ndarray):
        edges = dict(zip(range(len(edges)),edges))


    property_types = {}
    property_types['bool'] = "int"
    property_types['int'] = "int"
    property_types['int32'] = "int"
    property_types['int64'] = "int"
    property_types['float'] = "float"
    property_types['float32'] = "float"
    property_types['float64'] = "float"
    # property_types['object'] = "list"

    start_time =time()
    print "--> Saving .ply"

    ply_file = open(ply_filename,'w+')

    ply_file.write("ply\n")
    ply_file.write("format ascii 1.0\n")

    # Declaring vertices and vertex properties
    ply_file.write("element vertex "+str(len(positions))+"\n")
    ply_file.write("property float x\n")
    ply_file.write("property float y\n")
    ply_file.write("property float z\n")

    vertex_property_data = {}
    for property_name in vertex_properties.keys():
        property_data = vertex_properties[property_name]
        if isinstance(property_data,dict):
            property_data = np.array([property_data[p] for p in positions.keys()])
        elif isinstance(property_data,list) or isinstance(property_data,tuple):
            property_data = np.array(property_data)
        vertex_property_data[property_name] = property_data

        if property_data.ndim == 1:
            property_type = property_types[str(property_data.dtype)]
        elif property_data.ndim == 2:
            property_type = "list int "+property_types[str(property_data.dtype)]
        else:
            property_type = "tensor "
            for i in xrange(property_data.ndim-1):
                property_type += "int "
            property_type += property_types[str(property_data.dtype)]

        ply_file.write("property "+property_type+" "+property_name+"\n")

    ply_file.write("element face "+str(len(triangles))+"\n")
    ply_file.write("property list int int vertex_index\n")

    triangle_property_data = {}
    for property_name in triangle_properties.keys():
        property_data = triangle_properties[property_name]
        if isinstance(property_data,dict):
            property_data = np.array([property_data[t] for t in triangles.keys()])
        elif isinstance(property_data,list) or isinstance(property_data,tuple):
            property_data = np.array(property_data)
        triangle_property_data[property_name] = property_data

        if property_data.ndim == 1:
            property_type = property_types[str(property_data.dtype)]
        elif property_data.ndim == 2:       
            property_type = "list int "+property_types[str(property_data.dtype)]
        else:
            property_type = "tensor "
            for i in xrange(property_data.ndim-1):
                property_type += "int "
            property_type += property_types[str(property_data.dtype)]
        ply_file.write("property "+property_type+" "+property_name+"\n")
    ply_file.write("end_header\n")

    # Writing property data
    vertex_index = {}
    for pid, p in enumerate(positions.keys()):
        ply_file.write(str(positions[p][0])+" ")
        ply_file.write(str(positions[p][1])+" ")
        ply_file.write(str(positions[p][2]))
        for property_name in vertex_properties.keys():
            data = np.array(vertex_property_data[property_name][pid])
            if data.ndim == 0:
                ply_file.write(" "+str(data))
            else:
                ply_file.write(multidim_data_ply_string(data))
        ply_file.write("\n")
        vertex_index[p] = pid

    for tid, t in enumerate(triangles.keys()):
        ply_file.write("3 ")
        ply_file.write(str(vertex_index[triangles[t][0]])+" ")
        ply_file.write(str(vertex_index[triangles[t][1]])+" ")
        ply_file.write(str(vertex_index[triangles[t][2]]))
        for property_name in triangle_properties.keys():
            data = np.array(triangle_property_data[property_name][tid])
            if data.ndim == 0:
                ply_file.write(" "+str(data))
            else:
                ply_file.write(multidim_data_ply_string(data))
        ply_file.write("\n")
    ply_file.flush()
    ply_file.close()

    end_time = time()
    print "<-- Saving .ply        [",end_time-start_time,"s]"


def multidim_data_ply_string(data):
    data_string = ""
    for dim in xrange(data.ndim):
        data_string += " "+str(data.shape[dim])
    # print data," (",len(data),")"
    # data_string = " "+str(len(data))
    for d in data.ravel():  
        data_string += " "+str(d)
        # data_string += (" "+str(d)) if data.ndim==1 else multidim_data_ply_string(d)
    return data_string


