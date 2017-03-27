# -*- python -*-
#
#       Openalea.CellComplex
#
#       Copyright 2006-2009 INRIA - CIRAD - INRA
#
#       File author(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################

import numpy as np
from scipy import ndimage as nd
from array import array

from openalea.container import PropertyTopomesh, array_dict
from openalea.container.utils import IdDict


class TemporalPropertyTopomesh(PropertyTopomesh):

    def __init__(self, degree=3, topomesh=None, **kwds):
        self._successors = [IdDict(idgenerator = "set") for i in xrange(degree+1)]
        self._predecessors = [IdDict(idgenerator = "set") for i in xrange(degree+1)]

        PropertyTopomesh.__init__(self, degree, topomesh, **kwds)
        if topomesh is not None:
            for d in xrange(degree+1):
                for w in self.wisps(d):
                    self._successors[d][w] = array("L")
                    self._predecessors[d][w] = array("L")

    def add_wisp(self, degree, wid=None):
        wid = super(PropertyTopomesh, self).add_wisp(degree,wid)
        self._successors[degree][wid] = array("L")
        self._predecessors[degree][wid] = array("L")
        return wid


    def relate(self, degree, wid, successor_wid):
        self._successors[degree][wid].append(successor_wid)
        self._predecessors[degree][successor_wid].append(wid)

    def unrelate(self, degree, wid, successor_wid):
        self._successors[degree][wid].remove(successor_wid)
        self._predecessors[degree][successor_wid].remove(wid)

    def successors(self, degree, wid):
        return iter(self._successors[degree][wid])

    def predecessors(self, degree, wid):
        return iter(self._predecessors[degree][wid])

    def has_predecessor(self, degree, wid):
        return len(self._predecessors[degree][wid])>0

    def _ancestors(self, degree, wids):
        ret = set()
        for wid in wids :
            if self.has_predecessor(degree,wid):
                ret |= set(self._ancestors(degree,self.predecessors(degree,wid)))
            else:
                ret |= {wid}
        return iter(ret)

    def ancestors(self, degree, wid):
        return self._ancestors(degree,[wid])


def compute_temporal_topomesh_property(topomesh, property_name, degree=0, positions=None, verbose=False):
    """
    """
    from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property

    if positions is None:
        positions = topomesh.wisp_property('barycenter',degree=0)
    
    if property_name == 'time':
        assert degree>0
        if not 'time' in topomesh.wisp_property_names(degree):
            topomesh.add_wisp_property('time',degree=degree)
        compute_topomesh_property(topomesh,'vertices',degree)
        topomesh.update_wisp_property('time',degree=degree,values=array_dict([np.unique(topomesh.wisp_property('time',0).values(topomesh.wisp_property('vertices',degree)[w]))[0] for w in topomesh.wisps(degree)],keys=list(topomesh.wisps(degree))))
       
    if property_name == 'mother_cell':
        assert degree==3
        if not 'mother_cell' in topomesh.wisp_property_names(degree):
            topomesh.add_wisp_property('mother_cell',degree=degree)
        topomesh.update_wisp_property('mother_cell',degree=degree,values=np.array([int(list(topomesh.ancestors(degree,w))[0]) for w in topomesh.wisps(degree)]),keys=np.array(list(topomesh.wisps(degree))))
       
       

