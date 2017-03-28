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
from openalea.oalab.mimedata import QMimeCodecPlugin
from openalea.core.authors import gcerutti


@PluginDef
class MeshFileCodecPlugin(QMimeCodecPlugin):
    authors = [gcerutti]

    qtdecode = [
        ('text/uri-list', 'openalea/interface.ITopomesh'),
        ('text/uri-list', 'text/plain'),
    ]
    qtencode = [
        ('cellcomplex/property_topomesh', 'openalea/interface.ITopomesh')
    ]

    mimetype_desc = {
        'openalea/interface.ITopomesh': dict(title='PropertyTopomesh Object'),
        'text/plain': dict(title='Plain Text'),
    }

    def __call__(self):
        from openalea.cellcomplex.mesh_oalab.mimedata.codec import MeshFileCodec
        return MeshFileCodec


@PluginDef
class DataFrameFileCodecPlugin(QMimeCodecPlugin):
    authors = [gcerutti]

    qtdecode = [
        ('text/uri-list', 'pandas/dataframe'),
        ('text/uri-list', 'text/plain'),
    ]

    mimetype_desc = {
        'pandas/dataframe': dict(title='Pandas DataFrame Object'),
        'text/plain': dict(title='Plain Text'),
    }

    def __call__(self):
        from openalea.cellcomplex.mesh_oalab.mimedata.codec import DataFrameFileCodec
        return DataFrameFileCodec
