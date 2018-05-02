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

from openalea.cellcomplex.property_topomesh.property_topomesh_optimization import property_topomesh_isotropic_remeshing
from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property
from openalea.cellcomplex.property_topomesh.property_topomesh_extraction import star_interface_topomesh
from openalea.cellcomplex.property_topomesh.example_topomesh import hexagonal_grid_topomesh

def test_isotropic_remeshing():
    topomesh = star_interface_topomesh(hexagonal_grid_topomesh(size=2))
    compute_topomesh_property(topomesh,'length',1)

    optimized_topomesh = property_topomesh_isotropic_remeshing(topomesh,maximal_length=1.,iterations=4)
    compute_topomesh_property(optimized_topomesh,'vertices',1)
    compute_topomesh_property(optimized_topomesh,'length',1)

    assert topomesh.wisp_property('length',1).values().mean() == 1.
    assert optimized_topomesh.wisp_property('length',1).values().max() < 1.
