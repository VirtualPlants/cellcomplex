# -*- coding: utf-8 -*-
# -*- python -*-
#
#       PropertyTopomesh
#
#       Copyright 2015 INRIA - CIRAD - INRA
#
#       File author(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       File contributor(s): Guillaume Baty <guillaume.baty@inria.fr>, 
#                            Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       TissueLab Website : http://virtualplants.github.io/
#
###############################################################################

__revision__ = ""

from openalea.core.plugin import PluginDef
from openalea.core.authors import gcerutti


class WorldPlugin(object):
    implement = 'IWorldHandler'



class TopomeshControl(WorldPlugin):
    label = 'Topomesh Handler'
    icon = 'topomesh_control.png'
    authors = [gcerutti]
    __plugin__ = True

    def __call__(self):
        from openalea.cellcomplex.mesh_oalab.world.property_topomesh_handler import TopomeshHandler
        return TopomeshHandler

@PluginDef
class DataframeHandlerPlugin(WorldPlugin):
    label = 'Dataframe Handler'
    icon = 'dataframe_control.png'
    authors = [gcerutti]
    __plugin__ = True

    def __call__(self):
        from openalea.cellcomplex.mesh_oalab.world.dataframe_handler import DataframeHandler
        return DataframeHandler
