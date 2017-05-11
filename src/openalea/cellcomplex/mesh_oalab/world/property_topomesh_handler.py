# -*- coding: utf-8 -*-
# -*- python -*-
#
#       PropertyTopomesh
#
#       Copyright 2015-2016 INRIA - CIRAD - INRA
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

# import weakref
from openalea.vpltk.qt import QtGui, QtCore
from openalea.core.observer import AbstractListener
from openalea.core.control import Control
from openalea.oalab.control.manager import ControlManagerWidget
from openalea.core.service.ipython import interpreter as get_interpreter

from openalea.oalab.service.drag_and_drop import add_drop_callback
from openalea.oalab.widget.world import WorldModel

from openalea.container import PropertyTopomesh, array_dict

try:
    from openalea.cellcomplex.property_topomesh.temporal_property_topomesh import TemporalPropertyTopomesh
    from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh
    from openalea.cellcomplex.property_topomesh.property_topomesh_io import save_ply_property_topomesh
except:
    print "Openalea.Cellcomplex must be installed to use TopomeshHandler!"
    raise


import numpy as np
from scipy.cluster.vq import vq

from tissuelab.gui.vtkviewer.vtkworldviewer import setdefault, world_kwargs


element_names = dict(zip(range(4),['vertices','edges','faces','cells']))

cst_proba = dict(step=0.01, min=0, max=1)
cst_degree = dict(step=1,min=0,max=3)

attribute_definition = {}
attribute_definition['topomesh'] = {}
for degree in xrange(4):
    attribute_definition['topomesh']["display_"+str(degree)] = dict(value=False,interface="IBool",constraints={},label="Display "+element_names[degree])
    attribute_definition['topomesh']["property_degree_"+str(degree)] = dict(value=degree,interface="IInt",constraints=cst_degree,label="Degree") 
    attribute_definition['topomesh']["property_name_"+str(degree)] = dict(value="",interface="IEnumStr",constraints=dict(enum=[""]),label="Property")     
    attribute_definition['topomesh']["coef_"+str(degree)] = dict(value=1,interface="IFloat",constraints=cst_proba,label="Coef") 
attribute_definition['topomesh']["cell_edges"] = dict(value=False,interface="IBool",constraints={},label="Cell edges")

attribute_definition['temporal_topomesh'] = {}
attribute_definition['temporal_topomesh']['time_point'] = dict(value=0,interface="IFloat",constraints=dict(min=0,max=0,step=1),label="Time Point") 


def _property_names(world_object, attr_name, property_name, **kwargs):
    degree = int(attr_name[-1])
    property_degree = world_object["property_degree_"+str(degree)]
    print "New property_names ! ",property_degree
    topomesh = world_object.data
    constraints = dict(enum=[""]+list(topomesh.wisp_property_names(property_degree)))
    print constraints
    if property_name in constraints['enum']:
        return dict(value=property_name, constraints=constraints)
    else:
        return dict(value="", constraints=constraints)

def _time_points(world_object, attr_name, time_point, **kwargs):
    topomesh = world_object.data
    if topomesh.has_wisp_property('time',0):
        time_points = np.sort(np.unique(topomesh.wisp_property('time',0).values())).astype(float)
    else:
        time_points = np.zeros(1).astype(float)
    constraints = dict(min=time_points.min(),max=time_points.max(),step=(time_points.max()-time_points.min())/100.)
    print constraints
    value = np.maximum(np.minimum(time_point,constraints['max']),constraints['min'])
    return dict(value=value, constraints=constraints)


class TopomeshHandler(AbstractListener):

    element_names = dict(zip(range(4),['vertices','edges','faces','cells']))

    property_colormaps = {}
    property_colormaps['cells'] = 'glasbey'
    property_colormaps['volume'] = 'morocco'
    property_colormaps['eccentricity'] = 'jet'
    
    def __init__(self, parent=None, style=None):
        AbstractListener.__init__(self)

        self.world = None
        self.model = WorldModel()

        self._mesh = {}
        self._mesh_matching = {}


        self._current = None

        self.interpreter = get_interpreter()
        self.interpreter.locals['topomesh_control'] = self


    def initialize(self):
        from openalea.core.world.world import World
        from openalea.core.service.ipython import interpreter
        world = World()
        world.update_namespace(interpreter())
        self.set_world(world)

    def set_world(self, world):
        self.clear()

        self.world = world
        self.world.register_listener(self)

        for object_name in world.keys():
            if isinstance(world[object_name].data,PropertyTopomesh):
                self.refresh_world_object(world[object_name])
    
    def notify(self, sender, event=None):
        signal, data = event
        if signal == 'world_changed':
            world, old_object, new_object = data
            if isinstance(new_object.data,PropertyTopomesh):
                self.refresh()
        elif signal == 'world_object_removed':
            world, old_object = data
            if isinstance(old_object.data,PropertyTopomesh):
                for degree in xrange(4):
                    if world.has_key(old_object.name+"_"+self.element_names[degree]):
                        world.remove(old_object.name+"_"+self.element_names[degree])
                self.refresh()
        elif signal == 'world_object_changed':
            world, old_object, world_object = data
            if isinstance(world_object.data,PropertyTopomesh):

                print world_object.data,": ",isinstance(world_object.data,PropertyTopomesh), isinstance(world_object.data,TemporalPropertyTopomesh)
                # raw_input()

                self.refresh_world_object(world_object)
        elif signal == 'world_object_item_changed':
            world, world_object, item, old, new = data
            if isinstance(world_object.data,PropertyTopomesh):
                # self.refresh_manager(world_object)
                if item == 'attribute':
                    self.update_topomesh_display(world_object, new)
        elif signal == 'world_sync':
            self.refresh()


    def clear(self):
        if self.world:
            self.world.unregister_listener(self)
            self.world = None

    def refresh_world_object(self, world_object):
        if world_object:
            dtype = 'topomesh'

            temporal = isinstance(world_object.data,TemporalPropertyTopomesh)

            self._topomesh = world_object.data

            kwargs = world_kwargs(world_object)

            print "Set default attributes : ",world_object.name

            if temporal:
                setdefault(world_object, 'temporal_topomesh', 'time_point', conv=_time_points, attribute_definition=attribute_definition, **kwargs)
            
            world_object.silent = True
            for degree in np.arange(4)[::-1]:
                setdefault(world_object, dtype, 'display_'+str(degree), attribute_definition=attribute_definition, **kwargs)

                setdefault(world_object, dtype, 'property_degree_'+str(degree), attribute_definition=attribute_definition, **kwargs)
                setdefault(world_object, dtype, 'property_name_'+str(degree), conv=_property_names, attribute_definition=attribute_definition, **kwargs)
                if degree>1:
                    setdefault(world_object, dtype, 'coef_'+str(degree), attribute_definition=attribute_definition, **kwargs)
                elif degree == 1:
                    setdefault(world_object, dtype, 'cell_edges', attribute_definition=attribute_definition, **kwargs)
            world_object.silent = False
            
            if not self._mesh.has_key(world_object.name):
                self._mesh[world_object.name] = dict([(0,None),(1,None),(2,None),(3,None)])
                self._mesh_matching[world_object.name] = dict([(0,None),(1,None),(2,None),(3,None)])

            world_object.set_attribute("display_"+str(max([degree for degree in xrange(4) if world_object.data.nb_wisps(degree)>0])),True)
    
    def select_world_object(self, object_name):
        if object_name != self._current:
            self._current = object_name

    def refresh_item(self, world_object, item, old, new):
        object_name = world_object.name

    def refresh(self):
        if self.world is not None:
            self.set_world(self.world)

    def update_topomesh_display(self, world_object, attribute):
        if world_object:

            temporal = isinstance(world_object.data,TemporalPropertyTopomesh)

            if 'time' in attribute['name']:
                topomesh = world_object.data

                for display_degree in xrange(4):
                    if world_object['display_'+str(display_degree)]:
                        property_name = world_object['property_name_'+str(display_degree)]
                        property_degree = world_object['property_degree_'+str(display_degree)]
                        cell_edges = world_object['cell_edges']

                        if temporal:
                            time = world_object['time_point']
                            elements_times = topomesh.wisp_property('time',display_degree).values()
                            topomesh_times = np.unique(elements_times)
                            display_time = topomesh_times[vq(np.array([time]),topomesh_times)[0][0]]
                            print "Displaying time ",display_time
                            wids = np.array(list(topomesh.wisps(display_degree)))[elements_times==display_time]
                        else:
                            wids = None

                        if display_degree > 1:
                            coef = world_object['coef_'+str(display_degree)]
                        else:
                            coef = 1
                        print "Property : ",property_name," (",attribute['name'],")"
                        mesh, matching = topomesh_to_triangular_mesh(topomesh,degree=display_degree,coef=coef,wids=wids,mesh_center=[0,0,0],cell_edges=cell_edges,property_name=property_name,property_degree=property_degree)
                        
                        self._mesh[world_object.name][display_degree] = mesh
                        self._mesh_matching[world_object.name][display_degree] = matching

                        if self.world.has_key(world_object.name+"_"+self.element_names[display_degree]):
                            kwargs = world_kwargs(self.world[world_object.name+"_"+self.element_names[display_degree]])
                            if not 'coef_' in attribute['name']:
                                if kwargs.has_key('intensity_range'):
                                    kwargs.pop('intensity_range')
                        else:
                            kwargs = {}
                            kwargs['colormap'] = 'glasbey' if ((property_name == '')and(display_degree>0)) else self.property_colormaps.get(property_name,'grey')
                            # kwargs['position'] = world_object['position']


                        self.world.add(mesh,world_object.name+"_"+self.element_names[display_degree],**kwargs)
                    else:
                        self.world.remove(world_object.name+"_"+self.element_names[display_degree])


            if 'display_' in attribute['name'] or 'coef_' in attribute['name']:
                display_degree = int(attribute['name'][-1])
                if world_object['display_'+str(display_degree)]:
                    topomesh = world_object.data
                    property_name = world_object['property_name_'+str(display_degree)]
                    property_degree = world_object['property_degree_'+str(display_degree)]
                    cell_edges = world_object['cell_edges']

                    if temporal:
                        time = world_object['time_point']
                        
                        elements_times = topomesh.wisp_property('time',display_degree).values()
                        topomesh_times = np.unique(elements_times)
                        display_time = topomesh_times[vq(np.array([time]),topomesh_times)[0][0]]
                        print "Displaying time ",display_time
                        wids = np.array(list(topomesh.wisps(display_degree)))[elements_times==display_time]
                    else:
                        wids = None

                    if display_degree > 1:
                        coef = world_object['coef_'+str(display_degree)]
                    else:
                        coef = 1
                    print "Property : ",property_name," (",attribute['name'],")"
                    mesh, matching = topomesh_to_triangular_mesh(topomesh,degree=display_degree,coef=coef,wids=wids,mesh_center=[0,0,0],cell_edges=cell_edges,property_name=property_name,property_degree=property_degree)
                    
                    self._mesh[world_object.name][display_degree] = mesh
                    self._mesh_matching[world_object.name][display_degree] = matching

                    if self.world.has_key(world_object.name+"_"+self.element_names[display_degree]):
                        kwargs = world_kwargs(self.world[world_object.name+"_"+self.element_names[display_degree]])
                        if not 'coef_' in attribute['name']:
                            if kwargs.has_key('intensity_range'):
                                kwargs.pop('intensity_range')
                    else:
                        kwargs = {}
                        kwargs['colormap'] = 'glasbey' if ((property_name == '')and(display_degree>0)) else self.property_colormaps.get(property_name,'grey')
                        # kwargs['position'] = world_object['position']

                    self.world.add(mesh,world_object.name+"_"+self.element_names[display_degree],**kwargs)
                else:
                    self.world.remove(world_object.name+"_"+self.element_names[display_degree])

            elif 'property_name_' in attribute['name']:
                display_degree = int(attribute['name'][-1])
                if world_object['display_'+str(display_degree)]:
                    topomesh = world_object.data
                    property_name = world_object['property_name_'+str(display_degree)]
                    property_degree = world_object['property_degree_'+str(display_degree)]
                    mesh_element_matching = self._mesh_matching[world_object.name][display_degree][property_degree]
                    if topomesh.has_wisp_property(property_name,property_degree):
                        property_data = array_dict(topomesh.wisp_property(property_name,property_degree).values(mesh_element_matching.values()),mesh_element_matching.keys())
                    else:
                        property_data = array_dict()
                    if property_degree == 0:
                        self._mesh[world_object.name][display_degree].point_data = property_data.to_dict()
                    else:
                        self._mesh[world_object.name][display_degree].triangle_data = property_data.to_dict()

                    self.world[world_object.name+"_"+self.element_names[display_degree]].data = self._mesh[world_object.name][display_degree]
                    if len(property_data)>1:
                        self.world[world_object.name+"_"+self.element_names[display_degree]].set_attribute('intensity_range',(property_data.values().min(),property_data.values().max()))

            elif 'property_degree_' in attribute['name']:
                dtype = 'topomesh'
                kwargs = world_kwargs(world_object)
                display_degree = int(attribute['name'][-1])
                print world_object['property_degree_'+str(display_degree)]
                world_object.silent = True
                setdefault(world_object, dtype, 'property_name_'+str(display_degree), conv=_property_names, attribute_definition=attribute_definition, **kwargs)
                world_object.silent = False





