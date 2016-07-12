
Computing Geometrical Properties on a PropertyTopomesh
======================================================

In CellComplex, cell tissues are represented as meshes (more exactly
cellular complexes implemented as incidence graphs) in an object called
**PropertyTopomesh**. This data structure contains the geometry and
topology of the cells, and can also bear other properties defined on
each of its elements (vertices, edges, faces or cells).

The CellComplex.PropertyTopomesh library comes with tools to
automatically compute a whole range of geometrical properties on the
tissue mesh (*compute\_topomesh\_property*) and to propagate properties
between elements of different dimensions
(*compute\_topomesh\_vertex\_property\_from\_faces* for instance).

.. code:: python

    import numpy as np
    import os
    
    from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property
    from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_vertex_property_from_faces
    from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_cell_property_from_faces

It is possible to construct a **PropertyTopomesh** from standard mesh
file formats, such as PLY. Here an example of a meristematic L1 tissue
saved as a PLY file, provided as an example in the *share/data*
directory of the package, is loaded.

.. code:: python

    import openalea.cellcomplex.property_topomesh
    from openalea.deploy.shared_data import shared_data
    dirname = shared_data(openalea.cellcomplex.property_topomesh)
    filename = os.path.join(dirname,"p194-t4_L1_topomesh.ply")
    
    from openalea.cellcomplex.property_topomesh.property_topomesh_io import read_ply_property_topomesh
    topomesh = read_ply_property_topomesh(filename,verbose=False)


.. parsed-literal::

    property float32 x
    
    property float32 y
    
    property float32 z
    
    property list uchar int32 vertex_indices
    
    property uchar red
    
    property uchar green
    
    property uchar blue
    
    property float32 area
    
    property int epidermis
    
    property int source
    
    property int target
    
    property list int int face_index
    
    property float length
    
    property list int int face_index
    
    property int label
    
    property barycenter is undefined on elements of degree 0
    Creating property  barycenter  for degree  0
    property length is undefined on elements of degree 1
    Creating property  length  for degree  1
    property area is undefined on elements of degree 2
    Creating property  area  for degree  2
    property epidermis is undefined on elements of degree 2
    Creating property  epidermis  for degree  2


The tissue is a surfacic triangular mesh where triangles are linked to a
labelled cell. Each label can be mapped to an independent color for
visualization.

.. code:: python

    # Specific function for 3D visualization in notebooks
    from vplants.meshing.notebook_tools import vtk_show_polydata
    from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh
    
    mesh,_ = topomesh_to_triangular_mesh(topomesh,3,coef=0.95)
    vtk_show_polydata(mesh._repr_vtk_(), width=800, height=300, position=(0,-50,-190))


.. parsed-literal::

    --> Creating triangular mesh
    <-- Creating triangular mesh [ 0.267352104187 s]




.. image:: output_5_1.png



Based on this initial mesh, we would like to compute the curvature of
the surface. To do so, we will use the function
*compute\_topomesh\_property* to compute the geometrical properties
allowing to obtain this information.

The first property to compute is the **normal** of the faces (elements
of degree = 2). Given the implementation as an incidence graph, the
order of vertices for each face is not unambiguous, hence the
orientation of faces is not contained in the surface. It is necessary to
have a way of re-orienting consistently the normals. Here we choose the
*orientation* mode that propagates the orientation infomration in a
topologically accurate way on a connected surfacic mesh.

Next, we compute the **area** of all faces, and use it to compute the
**normal** property on the vertices (elements of degree = 0). This is
done by averaging the normals of faces neighboring each vertex using the
function *compute\_topomesh\_vertex\_property\_from\_faces*.

Once each vertex bears a consistent normal, it is possible to compute
the **mean curvature** of the each face (actually the **3D curvature
tensor** containing all the information on the principal curvatures, in
particular mean, gaussian, min and max curvatures). This property can be
propagated to the vertices using the same method as for the normals.

.. code:: python

    compute_topomesh_property(topomesh,'normal',2,normal_method='orientation')
    compute_topomesh_property(topomesh,'area',2)
    compute_topomesh_vertex_property_from_faces(topomesh,'normal',weighting='area',adjacency_sigma=1.2,neighborhood=3)
    compute_topomesh_property(topomesh,'mean_curvature',2)
    compute_topomesh_vertex_property_from_faces(topomesh,'mean_curvature',weighting='area',adjacency_sigma=1.2,neighborhood=3)


.. parsed-literal::

    --> Computing vertex property from faces
    <-- Computing vertex property from faces [ 1.45693802834 s]
    --> Computing vertex property from faces
    <-- Computing vertex property from faces [ 1.46020102501 s]


This numerical property can be visualized on the triangles of the mesh
using the property on the vertices, creating a smooth interpolation of
curvature values.

.. code:: python

    mesh,_ = topomesh_to_triangular_mesh(topomesh,2,property_name='mean_curvature',property_degree=0)
    vtk_show_polydata(mesh._repr_vtk_(), width=800, height=300, colormap_name='curvature', 
                      intensity_range=(-0.05,0.05), position=(0,-50,-190))


.. parsed-literal::

    --> Creating triangular mesh
    <-- Creating triangular mesh [ 0.50506901741 s]




.. image:: output_9_1.png



Alternatively, the mean curvature could be propagated to the cells
(elements of degree = 3) to have unique value for each cell of the
tissue. This value can then be conveniently visualized on a
representation of the mesh that specifically isolates the cells.

.. code:: python

    compute_topomesh_cell_property_from_faces(topomesh,'mean_curvature')
    
    mesh,_ = topomesh_to_triangular_mesh(topomesh,3,coef=0.95,property_name='mean_curvature')
    vtk_show_polydata(mesh._repr_vtk_(), width=800, height=300, colormap_name='curvature', 
                      intensity_range=(-0.05,0.05), position=(0,-50,-190))


.. parsed-literal::

    --> Creating triangular mesh
    <-- Creating triangular mesh [ 0.282058000565 s]




.. image:: output_11_1.png



All the properties are stored inside the topomesh that can then be saved
with only the desired properties in a standard PLY file.

.. code:: python

    save_filename = os.path.join(dirname,"p194-t4_L1_topomesh_curvature.ply")
    
    from openalea.cellcomplex.property_topomesh.property_topomesh_io import save_ply_property_topomesh
    save_ply_property_topomesh(topomesh,save_filename,
                               properties_to_save=dict([(0,[]),(1,[]),(2,[]),(3,['mean_curvature'])]))



.. parsed-literal::

    --> Saving .ply
    <-- Saving .ply        [ 0.177224874496 s]


