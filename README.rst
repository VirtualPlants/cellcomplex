============================
OpenAlea CellComplex Library
============================

.. {# pkglts, doc

.. #}

The Topological Plant Tissue datastructure

Authors:
--------
* Frederic Boudon (frederic.boudon@cirad.fr)
* Guillaume Cerutti (guillaume.cerutti@inria.fr)


Institutes:
-----------

* Inria (http://www.inria.fr)
* Cirad (http://www.cirad.fr)


License: 
--------

* `Cecill-C`


Description
-----------

OpenAlea.CellComplex is a library providing data structures and algorithms to represent and analyze plant cell tissues in 3D. It offers an implementation of the topological structure of cellular complex as an Incidence Graph in a class named **PropertyTopomesh**.

.. image:: ../tissue.png

The structure comes with algorithms allowing to
	* Create a structure from more basic representations
	* Compute geometrical and topological properties and store them within the structure
	* Edit the structure by local topological operations
	* Read and export the structure from/to a standard PLY format (http://sainsburyworkshop2015.wikispaces.com/file/view/PlyFormat.pdf)


Mesh-OAlab
----------

A set of plugins and components for the best mesh experience in OpenAleaLab (and TissueLab)


The TopomeshControls service (currently an applet)
==================================================


Add the applet to your OALab environment as any regular applet :
	* In a workspace right click and select "Edit Layout"
	* Add a new tab (right click + "Add Tab") if necessary
	* Select the Topomesh Control applet in the scrolling list
	* Finalize your layout by right click and "Lock Layout"

Mesh objects stored as PropertyTopomesh structures can now be visualized simply by the command

.. code-block:: python

	world.add(topomesh,"topomesh")


VisuAlea components for mesh processing
=======================================

Add functionalities handling PropertyTopomesh objects directly as visual progamming bricks


Requirements
------------

* SconsX (https://github.com/openalea/sconsx)
* OpenAlea.Deploy (https://github.com/openalea/deploy)
* OpenAlea (https://github.com/openalea/openalea)
* OpenAlea.Container (https://github.com/openalea/openalea-components)
* NumPy / SciPy



Installation
------------

.. code-block:: python

	python setup.py develop


