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

import openalea.cellcomplex
from openalea.cellcomplex.property_topomesh.property_topomesh_io import read_ply_property_topomesh
from openalea.cellcomplex.property_topomesh.property_topomesh_io import read_obj_property_topomesh

from openalea.deploy.shared_data import shared_data

import os


def test_ply_reading():
    dirname = shared_data(openalea.cellcomplex)
    filename = os.path.join(dirname,"p194-t4_L1_topomesh.ply")
    topomesh = read_ply_property_topomesh(filename)

def test_obj_reading():
    dirname = shared_data(openalea.cellcomplex)
    filename = os.path.join(dirname,"icosahedron.obj")
    topomesh = read_obj_property_topomesh(filename)

    print topomesh.nb_wisps(0)
    assert topomesh.nb_wisps(0) == 12
    assert topomesh.nb_wisps(2) == 20
    


