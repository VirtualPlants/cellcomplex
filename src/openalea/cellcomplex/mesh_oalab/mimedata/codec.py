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

import mimetypes
from openalea.oalab.mimedata.qcodec import QMimeCodec

from openalea.core.path import path
from openalea.cellcomplex.property_topomesh import meshread
from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh

import sys

def is_mesh_path(path):
    mime, encoding = mimetypes.guess_type(path)
    if mime and mime.startswith('model/mesh'):
        return True
    if path.ext in ('.ply', '.obj', '.msh'):
        return True
    else:
        return False

def encode_mesh(category, data, mimetype_in, mimetype_out):
    # openalea://user@localhost:/project/projectname/data/dataname
    # uri = 'openalea://user@localhost:/openalea/topomesh/%s' % (data.__str__())
    print "Encoding Mesh!"
    return ('openalea/interface.ITopomesh', data)

def decode_mesh_file(filename, mimetype_in, mimetype_out):
    mesh = None

    if filename is not None:
        if filename.exists():
            try:
                mesh = meshread(filename)
            except:
                print "Error importing mesh file ",filename
                print sys.exc_info()
                mesh = None
    return mesh


class MeshFileCodec(QMimeCodec):

    def _raw_data(self, mimedata, mimetype_in, mimetype_out):
        """
        'text/uri-list' -> list of paths
        """

        if mimetype_in == 'text/uri-list':
            return [path(url.toLocalFile()) for url in mimedata.urls()]

    def quick_check(self, mimedata, mimetype_in, mimetype_out):
        raw_data = self._raw_data(mimedata, mimetype_in, mimetype_out)
        if not raw_data:
            return False
        elif mimetype_in == 'text/uri-list':
            url = raw_data[0]
        else:
            return False
        return is_mesh_path(url)


    def encode(self, data, mimetype_in, mimetype_out):
        return encode_mesh("property_topomesh", data, mimetype_in, mimetype_out)


    def qtdecode(self, mimedata, mimetype_in, mimetype_out):
        raw_data = self._raw_data(mimedata, mimetype_in, mimetype_out)
        if raw_data is None:
            return None, {}
        else:
            return self.decode(raw_data, mimetype_in, mimetype_out)

    def decode(self, raw_data, mimetype_in, mimetype_out):
        kwds = {}
        if raw_data[:4] == "file":
            local_file = str(raw_data)[7:-2]
        else:
            local_file = raw_data[0]
        local_file = path(local_file)
        kwds['name'] = local_file.namebase
        
        if mimetype_in == 'text/uri-list':
            data = decode_mesh_file(local_file, mimetype_in, mimetype_out)

            if mimetype_out == "openalea/interface.ITopomesh":
                return data, kwds
            elif mimetype_out == "text/plain":
                text = raw_data
                return text, kwds
            else:
                return data, kwds
        else:
            return None, {}


class TopomeshCodec(QMimeCodec):

    def quick_check(self, mimedata, mimetype_in, mimetype_out):
        if mimetype_in == "openalea/interface.ITopomesh":
            mesh = mimedata.data
        else:
            return False
        return mesh is not None

    def qtdecode(self, mimedata, mimetype_in, mimetype_out):
        return self.decode(mimedata, mimetype_in, mimetype_out)

    def decode(self, mimedata, mimetype_in, mimetype_out):
        if mimetype_in == 'openalea/interface.ITopomesh':
            mesh = mimedata.data
            kwds = {}

            if mimetype_out == "openalea/interface.ITopomesh":
                return mesh, kwds
            elif mimetype_out == "openalea/interface.ITriangularMesh":
                return topomesh_to_triangular_mesh(mesh)





def is_csv_path(path):
    mime, encoding = mimetypes.guess_type(path)
    if mime and mime.startswith('text/csv'):
        return True
    if path.ext in ('.csv'):
        return True
    else:
        return False

def decode_csv_file(filename, mimetype_in, mimetype_out):
    try:
        import pandas as pd
    except ImportError:
        print "Pandas should be installed to visualize CSV files"
        pass
    else:
        df = None
        if filename is not None:
            print filename[-4:]
            print "\""+filename+"\"",filename.exists()
            if filename.exists():
                if filename[-4:] == ".csv":
                    try:
                        df = pd.read_csv(filename)
                        if ";" in df.columns[0]:
                            df = pd.read_csv(filename,delimiter=";")
                    except:
                        df = None

    return df


class DataFrameFileCodec(QMimeCodec):

    def _raw_data(self, mimedata, mimetype_in, mimetype_out):
        """
        'text/uri-list' -> list of paths
        """

        if mimetype_in == 'text/uri-list':
            return [path(url.toLocalFile()) for url in mimedata.urls()]

    def quick_check(self, mimedata, mimetype_in, mimetype_out):
        raw_data = self._raw_data(mimedata, mimetype_in, mimetype_out)
        if not raw_data:
            return False
        elif mimetype_in == 'text/uri-list':
            url = raw_data[0]
        else:
            return False
        return is_csv_path(url)

    def qtdecode(self, mimedata, mimetype_in, mimetype_out):
        raw_data = self._raw_data(mimedata, mimetype_in, mimetype_out)
        if raw_data is None:
            return None, {}
        else:
            return self.decode(raw_data, mimetype_in, mimetype_out)

    def decode(self, raw_data, mimetype_in, mimetype_out):
        kwds = {}
        if raw_data[:4] == "file":
            local_file = str(raw_data)[7:-2]
        else:
            local_file = raw_data[0]
        local_file = path(local_file)
        kwds['name'] = local_file.namebase
        
        if mimetype_in == 'text/uri-list':
            data = decode_csv_file(local_file, mimetype_in, mimetype_out)

            if mimetype_out == "pandas/dataframe":
                return data, kwds
            elif mimetype_out == "text/plain":
                text = raw_data
                return text, kwds
            else:
                return data, kwds
        else:
            return None, {}

