
# Redirect path
import os

cdir = os.path.dirname(__file__)
pdir = os.path.join(cdir, "../cellcomplex/property_topomesh")
pdir = os.path.abspath(pdir)

__path__ = [pdir] + __path__[:]

from openalea.cellcomplex.property_topomesh.__init__ import *
