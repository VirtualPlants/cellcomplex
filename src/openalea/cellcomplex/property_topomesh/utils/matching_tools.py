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
from scipy.spatial import cKDTree

def kd_tree_match(obs, codebook):
    data = cKDTree(obs)
    res1 = data.query_ball_tree(cKDTree(codebook), 1e-5, 1, 1e-5)
    # print res1
    if np.array(res1).ndim == 2:
        res1 = np.array(res1)[:,0]
    else:
        res1 = np.array([r[0] if len(r)>0 else None for r in res1])
    #print res1
    #res2 = vq(obs, codebook)[0]
    #print res2
    #assert np.array_equal(res1,res2)
    return res1