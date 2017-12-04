# -*- coding: utf-8 -*-
# -*- python -*-
#
#       TissueLab - DataFrames
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

from openalea.vpltk.qt import QtGui, QtCore

from openalea.core.observer import AbstractListener
from openalea.core.control import Control
from openalea.oalab.control.manager import ControlManagerWidget
from openalea.core.service.ipython import interpreter as get_interpreter

from openalea.core.interface import IBool, IInt, IFloat, ITuple, IEnumStr
from openalea.oalab.interface import IIntRange, IColormap

from openalea.oalab.utils import ModalDialog

from openalea.oalab.service.drag_and_drop import add_drop_callback
from openalea.oalab.widget.world import WorldModel

from openalea.oalab.colormap.colormap_utils import Colormap

from copy import copy

try:
    import pandas as pd
except:
    print "Pandas must be installed to use DataframeControls! (conda install pandas)"
    raise

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.patches as patch
except:
    print "Matplotlib must be installed to use DataframeControls! (conda install matplotlib)"
    raise

try:
    import sklearn
except:
    print "SciKit-Learn should be installed to use PCA analysis tools!"
    pass

import numpy as np
import time

from cute_plot import simple_plot, smooth_plot, average_plot, histo_plot, bar_histo_plot, violin_plot, box_plot, density_contour_plot, map_plot
from pca_tools import pca_analysis

from tissuelab.gui.vtkviewer.vtkworldviewer import setdefault, world_kwargs, _colormap


cst_figure = dict(step=1,min=0,max=9)
cst_plots = dict(enum=['average curve','scatter','line','histogram','distribution','cumulative','boxplot','violin','density','map','PCA'])
cst_regression = dict(enum=['','linear','quadratic','cubic','invert','invert_quadratic','exponential','logarithmic','density','fireworks'])
cst_legend = dict(enum=['','top_right','bottom_right','top_left','bottom_left'])
cst_size = dict(step=5, min=0, max=200)
cst_width = dict(step=1, min=0, max=10)
cst_proba = dict(step=0.01, min=0, max=1)
cst_extent_range = dict(step=1, min=-50, max=150)
cst_points = dict(step=1, min=1, max=100)


dataframe_attributes = {}
dataframe_attributes['dataframe'] = {}
dataframe_attributes['dataframe']['figure'] = dict(value=0,interface="IInt",constraints=cst_figure,label=u"Figure Number")
dataframe_attributes['dataframe']['title'] = dict(value="",interface="IStr",constraints={},label=u"Figure Title")
for variable_type in ['X','Y','class','label','size','subplot','filter']:
    dataframe_attributes['dataframe'][variable_type+"_variable"] = dict(value="",interface="IEnumStr",constraints=dict(enum=[""]),label=variable_type.capitalize()+" Variable")  
    dataframe_attributes['dataframe'][variable_type+'_colormap'] = dict(value=dict(name='grey', color_points=dict([(0, (0, 0, 0)), (1, (1, 1, 1))])), interface=IColormap, label="Colormap")
    dataframe_attributes['dataframe'][variable_type+'_range'] = dict(value=(-1,101),interface=IIntRange,constraints=cst_extent_range,label=variable_type.capitalize()+" Range")
dataframe_attributes['dataframe']['plot'] = dict(value=0,interface="IEnumStr",constraints=cst_plots,label=u"Plot Type")

dataframe_attributes['dataframe']['markersize'] = dict(value=50,interface="IInt",constraints=cst_size,label=u"Marker Size")
dataframe_attributes['dataframe']['linewidth'] = dict(value=1,interface="IInt",constraints=cst_width,label=u"Linewidth")
dataframe_attributes['dataframe']['alpha'] = dict(value=1.0, interface=IFloat, constraints=cst_proba,label=u"Alpha")
dataframe_attributes['dataframe']['smooth_factor'] = dict(value=0.0, interface=IFloat, constraints=cst_proba,label=u"Smoothing")
dataframe_attributes['dataframe']['regression'] = dict(value='', interface=IEnumStr, constraints=cst_regression,label=u"Regression")
dataframe_attributes['dataframe']['n_points'] = dict(value=10, interface=IInt, constraints=cst_points,label=u"Number of Points")
dataframe_attributes['dataframe']['legend'] = dict(value='',interface="IEnumStr",constraints=cst_legend,label="Legend")

def _dataframe_columns(world_object, attr_name, variable_name, **kwargs):
    df = world_object.data
    constraints = dict(enum=[""]+[k for k in df.keys()])
    # print constraints
    if variable_name in constraints['enum']:
        return dict(value=variable_name, constraints=constraints)
    else:
        return dict(value="", constraints=constraints)


class DataframeHandler(AbstractListener):

    def __init__(self, parent=None, style=None):
        AbstractListener.__init__(self)

        self.world = None
        self.model = WorldModel()

        self.variables = {}

        self.interpreter = get_interpreter()
        self.interpreter.locals['dataframe_control'] = self

        self._figures = {}
        self._subplots = {}
        self._colorbar_axes = {}
        self._plot_legend = {}

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

        self.model.set_world(world)

        for object_name in world.keys():
            if isinstance(world[object_name].data,pd.DataFrame):
                self.refresh_world_object(world[object_name])

    def notify(self, sender, event=None):
        signal, data = event

        if signal == 'world_changed':
            world, old_object, new_object = data
            if isinstance(new_object.data,pd.DataFrame):
                self.refresh()
        elif signal == 'world_object_removed':
            world, old_object = data
            if isinstance(old_object.data,pd.DataFrame):
                figure = plt.figure(self._figures[old_object.name])
                figure.clf()
                del self._figures[old_object.name]
                self.refresh()
        elif signal == 'world_object_changed':
            world, old_object, world_object = data
            if isinstance(world_object.data,pd.DataFrame):
                print world_object.name," : ",world_object.kwargs
                self.refresh_world_object(world_object)
        elif signal == 'world_object_item_changed':
            world, world_object, item, old, new = data
            if isinstance(world_object.data,pd.DataFrame):
                # self.refresh_manager(world_object)
                if item == 'attribute':
                    if new['name'] in ['X_range','Y_range']:
                        self.update_plot_limits(world_object)
                    elif new['name'] in ['title']:
                        self.update_figure_title(world_object)
                    else:
                        self.update_dataframe_figure(world_object)
        elif signal == 'world_sync':
            self.refresh()

        elif signal == 'data_point_picked':
            fig, ax, point = data
            self.get_matching_data(fig,ax,point)

    def clear(self):
        if self.world:
            self.world.unregister_listener(self)
            self.world = None

        for o in self._figures.keys():
            plt.figure(self._figures[o]).clf()
            plt.draw()

    def refresh_world_object(self, world_object):
        if world_object:
            dtype = 'dataframe'

            dataframe = world_object.data
            kwargs = world_kwargs(world_object)
            for name,value in world_object.kwargs.items():
                kwargs[name] = value

            world_object.silent = True

            setdefault(world_object, dtype, 'figure', attribute_definition=dataframe_attributes, **kwargs)
            
            for variable_type in ['subplot']:
                setdefault(world_object, dtype, variable_type+'_variable', conv=_dataframe_columns, attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'title', attribute_definition=dataframe_attributes, **kwargs)
            
            for variable_type in ['filter']:
                setdefault(world_object, dtype, variable_type+'_variable', conv=_dataframe_columns, attribute_definition=dataframe_attributes, **kwargs)
                setdefault(world_object, dtype, variable_type+'_range', attribute_definition=dataframe_attributes, **kwargs)
             
            for variable_type in ['class']:
                setdefault(world_object, dtype, variable_type+'_variable', conv=_dataframe_columns, attribute_definition=dataframe_attributes, **kwargs)

            for variable_type in ['label']:
                setdefault(world_object, dtype, variable_type+'_variable', conv=_dataframe_columns, attribute_definition=dataframe_attributes, **kwargs)
                setdefault(world_object, dtype,  variable_type+'_colormap', conv=_colormap, attribute_definition=dataframe_attributes, **kwargs)
                setdefault(world_object, dtype, variable_type+'_range', attribute_definition=dataframe_attributes, **kwargs)

            for variable_type in ['X','Y']:
                setdefault(world_object, dtype, variable_type+'_variable', conv=_dataframe_columns, attribute_definition=dataframe_attributes, **kwargs)
                setdefault(world_object, dtype, variable_type+'_range', attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'plot', attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'markersize', attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'size_variable', conv=_dataframe_columns, attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'linewidth', attribute_definition=dataframe_attributes, **kwargs)
            setdefault(world_object, dtype, 'alpha', attribute_definition=dataframe_attributes, **kwargs)
            setdefault(world_object, dtype, 'smooth_factor', attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'regression', attribute_definition=dataframe_attributes, **kwargs)

            setdefault(world_object, dtype, 'n_points', attribute_definition=dataframe_attributes, **kwargs)

            world_object.silent = False

            setdefault(world_object, dtype, 'legend', attribute_definition=dataframe_attributes, **kwargs)

    def select_world_object(self, object_name):
        if object_name != self._current:
            self._current = object_name

    def refresh_item(self, world_object, item, old, new):
        object_name = world_object.name

    def refresh(self):
        if self.world is not None:
            self.set_world(self.world)

    def update_dataframe_figure(self, world_object):
        if world_object:
            # print "Display figure!"
            if self._figures.has_key(world_object.name):
                old_figure = plt.figure(self._figures[world_object.name])
                # old_figure.unregister_listener(self)
                old_figure.clf()
                plt.draw()

            self._figures[world_object.name] = world_object['figure']
            figure = plt.figure(self._figures[world_object.name])
            figure.clf()
            figure.register_listener(self)

            variables = [v for v in world_object.data.keys() if not 'Unnamed:' in v]
            variables = [v for v in variables if not np.array(world_object.data[v]).dtype == np.dtype('O')]
            all_data = np.transpose([world_object.data[variable] for variable in variables])

            subplot_variable = world_object['subplot_variable']
            if subplot_variable != "":
                subplots = np.array(world_object.data[subplot_variable])
            else:
                subplots = np.ones_like(world_object.data.index)
            if subplots.dtype != np.dtype('O') and (len(np.unique(subplots)) > 9):
                subplots = np.ones_like(world_object.data.index)

            X_variable = world_object['X_variable']
            if X_variable != "":
                all_X = np.array(world_object.data[X_variable])
            else:
                all_X = np.array(world_object.data.index)

            Y_variable = world_object['Y_variable']
            if Y_variable != "":
                all_Y = np.array(world_object.data[Y_variable])
            else:
                all_Y = np.array(world_object.data.index)

            size_variable = world_object['size_variable']
            if size_variable != "":
                all_sizes = np.array(world_object.data[size_variable])
                if all_sizes.dtype == np.dtype('O'):
                    all_sizes = np.ones_like(world_object.data.index)
            else:
                all_sizes = np.ones_like(world_object.data.index)

            class_variable = world_object['class_variable']
            if class_variable != "":
                all_classes = np.array(world_object.data[class_variable])
            else:
                all_classes = np.ones_like(world_object.data.index)

            label_variable = world_object['label_variable']
            if label_variable != "":
                all_labels = np.array(world_object.data[label_variable])
            else:
                all_labels = all_classes

            filter_variable = world_object['filter_variable']
            if filter_variable != "":
                all_filters = np.array(world_object.data[filter_variable])
            else:
                all_filters = all_classes


            if world_object['plot'] != 'PCA':
                point_invalidity = np.zeros(len(all_classes)).astype(bool)
                if subplots.dtype != np.dtype('O'):
                    point_invalidity = point_invalidity | np.isnan(subplots)
                else:
                    point_invalidity = point_invalidity | pd.isnull(subplots)

                point_invalidity = point_invalidity | np.isnan(all_X)
                point_invalidity = point_invalidity | np.isnan(all_Y)
                point_invalidity = point_invalidity | np.isnan(all_filters)
                point_invalidity = point_invalidity | np.isnan(all_sizes)
                if all_classes.dtype != np.dtype('O'):
                    point_invalidity = point_invalidity | np.isnan(all_classes)
                else:
                    point_invalidity = point_invalidity | pd.isnull(all_classes)

                if all_labels.dtype != np.dtype('O'):
                    point_invalidity = point_invalidity | np.isnan(all_labels)
                else:
                    point_invalidity = point_invalidity | pd.isnull(all_labels)
                # print "Invalid points : ",point_invalidity.sum()
                valid_points = np.where(True - point_invalidity)
            else:
                point_invalidity = np.zeros(len(all_classes)).astype(bool)
                for variable in variables:
                    if np.array(world_object.data[variable]).dtype != np.dtype('O'):
                        point_invalidity = point_invalidity | np.isnan(np.array(world_object.data[variable]))
                valid_points = np.where(True - point_invalidity)

            subplots = subplots[valid_points]
            all_data = all_data[valid_points]
            all_X = all_X[valid_points]
            all_Y = all_Y[valid_points]
            all_sizes = all_sizes[valid_points]
            all_classes = all_classes[valid_points]
            all_labels = all_labels[valid_points]
            all_filters = all_filters[valid_points]

            if all_filters.dtype != np.dtype('O'):
                filter_range = (all_filters.min() + world_object['filter_range'][0]*(all_filters.max()-all_filters.min())/100.,all_filters.min() + world_object['filter_range'][1]*(all_filters.max()-all_filters.min())/100.) 
                subplots = subplots[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
                all_X = all_X[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
                all_Y = all_Y[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
                all_sizes = all_sizes[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
                all_classes = all_classes[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
                all_labels = all_labels[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
                all_data = all_data[(all_filters>=filter_range[0])&(all_filters<=filter_range[1])]
            
            subplot_index = np.sort(np.unique(subplots))
            n_subplots = len(subplot_index)
            n_rows = 1 if n_subplots < 4 else 2 if n_subplots<8 else 3
            n_cols = np.ceil(n_subplots/float(n_rows))
            # subplot_label = "".join([w.capitalize()+" " for w in subplot_variable.split('_')])

            for sub, subplot in enumerate(subplot_index):
                figure.add_subplot(n_rows,n_cols,sub+1)
                if n_subplots>1:
                    # figure.gca().set_title(subplot_label + " = " + str(subplot))
                    # figure.gca().set_title("".join([w.capitalize()+" " for w in str(subplot).split('_')]))
                    figure.gca().set_title(str(subplot))

                X = all_X[subplots == subplot]
                Y = all_Y[subplots == subplot]
                sizes = all_sizes[subplots == subplot]
                classes = all_classes[subplots == subplot]
                labels = all_labels[subplots == subplot]
                data = all_data[subplots == subplot]

                n_slices = world_object['n_points']
                if classes.dtype != np.dtype('O') and (len(np.unique(classes)) > n_slices):
                    magnitude = np.power(10,np.around(4*np.log10(np.nanmean(classes)+np.nanstd(classes)+1e-7))/4+0.5)
                    magnitude = np.around(magnitude,int(-np.log10(magnitude))+1)
                    class_values = np.array(np.maximum(np.minimum(np.around(10*classes/magnitude),10),0),int)
                    classes = (np.arange(11)*magnitude)[class_values]

                label_variable = world_object['label_variable']
                if label_variable == "":
                    labels = classes

                # class_list = np.sort(np.unique(classes))
                class_list = np.sort(np.unique(all_classes))
                n_classes = len(class_list)

                # label_list = np.sort(np.unique(labels))
                label_list = np.sort(np.unique(all_labels))
                n_labels = len(label_list)

                cmap = Colormap(world_object['label_colormap']['name'])
                cmap._color_points = world_object['label_colormap']['color_points']
                cmap._compute()

                markersize = world_object['markersize']
                linewidth = world_object['linewidth']
                alpha = world_object['alpha']
                smooth = world_object['smooth_factor']

                xlabel = u"".join([w.capitalize()+" " for w in X_variable.split('_')])
                ylabel = u"".join([w.capitalize()+" " for w in Y_variable.split('_')])
                label_label = u"".join([w.capitalize()+" " for w in label_variable.split('_')])

                # if (n_classes > 1) or (n_labels == 1) or np.all(labels == classes) or (labels.dtype == np.dtype('O')):
                if (n_classes > 1) or (n_labels == 1) or np.all(labels == classes):
                    plot_legend = ['labels']


                    are_classes_sublabels = np.all([len(np.unique(labels[classes==c]))==1 for c in class_list])

                    if not are_classes_sublabels:
                        plot_legend += ['colorbar']
                        if labels.dtype == np.dtype('O'):
                            labels = np.sum([i*(labels==l).astype(int) for i,l in enumerate(label_list)],axis=0)+1
                            label_range = (0,labels.max()+1.)
                            cbar_ticks = np.arange(labels.max())+1
                            cbar_ticklabels = label_list
                            label_list = np.sort(np.unique(labels))
                        else:
                            label_range = (all_labels.min() + world_object['label_range'][0]*(all_labels.max()-all_labels.min())/100.,all_labels.min() + world_object['label_range'][1]*(all_labels.max()-all_labels.min())/100.)
                            #label_range = (labels.min() + world_object['label_range'][0]*(labels.max()-labels.min())/100.,labels.min() + world_object['label_range'][1]*(labels.max()-labels.min())/100.)
                            cbar_ticks = np.linspace(label_range[0],label_range[1],5)
                            cbar_ticklabels = np.round(cbar_ticks,decimals=3)
                        mpl_norm = mpl.colors.Normalize(vmin=label_range[0], vmax=label_range[1])
                        

                        color_dict = dict(red=[],green=[],blue=[])
                        for k,c in enumerate(['red','green','blue']):
                            color_dict[c] += [(0.0,0.0,0.0)]
                            color_dict[c] += [(1.0,1.0,1.0)]
                        for c in ['red','green','blue']:
                            color_dict[c] = tuple(color_dict[c])
                        mpl_cmap = mpl.colors.LinearSegmentedColormap('intensity', color_dict)


                    print class_list
                    print label_list

                    for i_class, c in enumerate(class_list):
                        if (classes == c).sum()>0:

                            label = np.unique(labels[classes==c])[0]
                            print label
                            i_label = np.where(label_list==label)[0][0]

                            if world_object['label_colormap']['name'] != 'glasbey':
                                if are_classes_sublabels:
                                    class_color = np.round(cmap.get_color((i_label+1)/float(n_labels+1)),decimals=3)
                                else:
                                    class_color = np.round(cmap.get_color((i_class+1)/float(n_classes+1)),decimals=3)
                            else:
                                if are_classes_sublabels:
                                    class_color = np.round(cmap.get_color(((i_label+1) % 256)/255.),decimals=3)
                                else:
                                    class_color = np.round(cmap.get_color(((i_class+1) % 256)/255.),decimals=3)


                            _,figure_labels = figure.gca().get_legend_handles_labels()
                            if are_classes_sublabels:
                                plot_label = str(label) if str(label) not in figure_labels else None
                            else:
                                plot_label = u""+str(c) if u""+str(c) not in figure_labels else None


                            if world_object['plot'] in ['scatter','line','density']:
                                if world_object['regression'] == 'fireworks':
                                    class_center = np.array([X[classes==c].mean(),Y[classes==c].mean()])
                                    [figure.gca().plot([x,class_center[0]],[y,class_center[1]],color=class_color,alpha=alpha/5.,linewidth=linewidth) for x,y in zip(X[classes==c],Y[classes==c])]
                            

                            if world_object['plot'] == 'scatter':
                                if are_classes_sublabels:
                                    simple_plot(figure,X[classes==c],Y[classes==c],class_color,xlabel=xlabel,ylabel=ylabel,linked=False,marker_size=markersize*sizes[classes==c]/sizes.mean(),linewidth=linewidth,alpha=alpha,label=plot_label)
                                else:
                                    class_labels = labels[classes==c]
                                    normalized_class_labels = np.maximum(0,np.minimum(1,(class_labels-label_range[0])/float(label_range[1]-label_range[0])))
                                    #label_colors = np.array([(1-0.8*np.abs(0.5-l)/0.5)*class_color + np.maximum(0,0.8*(0.5-l)/0.5)*np.zeros(3) + np.maximum(0,0.8*(l-0.5)/0.5)*np.ones(3) for l in normalized_class_labels])
                                    label_colors = np.array([l*class_color + (1-l)*(0.9*np.ones(3) + 0.1*class_color) for l in normalized_class_labels])
                                    

                                    simple_plot(figure,X[classes==c],Y[classes==c],label_colors,xlabel=xlabel,ylabel=ylabel,linked=False,marker_size=markersize*sizes[classes==c]/sizes.mean(),linewidth=linewidth,alpha=alpha,label=plot_label)
                                

                            elif world_object['plot'] == 'line':
                                simple_plot(figure,np.sort(X[classes==c]),Y[classes==c][np.argsort(X[classes==c])],class_color,xlabel=xlabel,ylabel=ylabel,linked=True,marker_size=0,linewidth=(linewidth+1)/2.,alpha=alpha,label=plot_label)
                                # smooth_plot(figure,np.sort(X[classes==c]),Y[classes==c][np.argsort(X[classes==c])],class_color,class_color,xlabel=xlabel,ylabel=ylabel,smooth_factor=smooth*Y.mean(),linewidth=linewidth,alpha=alpha,label=plot_label)
                            elif world_object['plot'] == 'average curve':
                                average_plot(figure,X[classes==c],Y[classes==c],class_color,xlabel=xlabel,ylabel=ylabel,linewidth=(linewidth+1)/2.,n_points=n_slices,smooth_factor=smooth,alpha=alpha,label=plot_label)
                            elif world_object['plot'] == 'histogram':
                                #x_range = (X.min() + world_object['X_range'][0]*(X.max()-X.min())/100.,X.min() + world_object['X_range'][1]*(X.max()-X.min())/100.)
                                x_range = (all_X.min() + world_object['X_range'][0]*(all_X.max()-all_X.min())/100.,all_X.min() + world_object['X_range'][1]*(all_X.max()-all_X.min())/100.)
                                bar_histo_plot(figure,X[classes==c],class_color,xlabel=xlabel,ylabel="Number of Elements",n_bins=n_slices,X_range=x_range,linewidth=linewidth,alpha=alpha,label=plot_label)
                            elif world_object['plot'] == 'distribution':
                                histo_plot(figure,X[classes==c],class_color,xlabel=xlabel,ylabel="Number of Elements (%)",cumul=False,bar=False,smooth_factor=smooth*10,linewidth=linewidth,alpha=alpha,label=plot_label)
                            elif world_object['plot'] == 'cumulative':
                                histo_plot(figure,X[classes==c],class_color,xlabel=xlabel,ylabel="Number of Elements (%)",cumul=True,bar=False,smooth_factor=smooth*10,linewidth=linewidth,alpha=alpha,label=plot_label)
                            elif world_object['plot'] == 'density':
                                #data_range = [[X.min()-0.5*(X.max()-X.min()),X.max()+0.5*(X.max()-X.min())],[Y.min()-0.5*(Y.max()-Y.min()),Y.max()+0.5*(Y.max()-Y.min())]]
                                data_range = [[all_X.min()-0.5*(all_X.max()-all_X.min()),all_X.max()+0.5*(all_X.max()-all_X.min())],[all_Y.min()-0.5*(all_Y.max()-all_Y.min()),all_Y.max()+0.5*(all_Y.max()-all_Y.min())]]
                                density_contour_plot(figure,X[classes==c],Y[classes==c],class_color,data_range,xlabel=xlabel,ylabel=ylabel,smooth_factor=smooth*10,marker_size=markersize*sizes[classes==c]/sizes.mean(),linewidth=linewidth,alpha=alpha,label=plot_label)

                            if world_object['plot'] in ['average curve','scatter','line','density']:
                                if world_object['regression'] in ['linear','quadratic','cubic','invert','invert_quadratic','exponential','logarithmic','density']:
                                    

                                    reg_X = np.linspace(X.min() -0.5*(X.max()-X.min()),X.min() +1.5*(X.max()-X.min()),400)
                                    if world_object['regression'] in ['linear','quadratic','cubic']:
                                        degree = 1 if world_object['regression']=='linear' else 2 if world_object['regression']=='quadratic' else 3
                                        p = np.polyfit(X[classes==c],Y[classes==c],deg=degree)
                                        reg_Y = np.polyval(p,reg_X)
                                    elif world_object['regression'] in ['invert']:
                                        p = np.polyfit(1./X[classes==c],Y[classes==c],deg=1)
                                        reg_Y = np.polyval(p,1./reg_X)
                                    elif world_object['regression'] in ['invert_quadratic']:
                                        p = np.polyfit(1./np.power(X[classes==c],0.5),Y[classes==c],deg=1)
                                        reg_Y = np.polyval(p,1./np.power(reg_X,0.5))
                                    elif world_object['regression'] in ['exponential']:
                                        p = np.polyfit(X[classes==c],np.log(Y[classes==c]),deg=1)
                                        reg_Y = np.exp(np.polyval(p,reg_X))
                                    elif world_object['regression'] in ['logarithmic']:
                                        p = np.polyfit(np.log(X[classes==c]),Y[classes==c],deg=1)
                                        reg_Y = np.polyval(p,np.log(reg_X))
                                    elif world_object['regression'] in ['density']:
                                        radius = 10.0*X.std()/float(len(reg_X))
                                        density_k = 100.*np.exp(-10.*smooth/2.)
                                        x_distances = np.power(np.power(X[classes==c][np.newaxis] - reg_X[:,np.newaxis],2),0.5)
                                        x_potential = 1./2. * (1. - np.tanh(density_k*(x_distances - radius)))
                                        # x_potential = np.transpose(x_potential)
                                        x_density = x_potential.sum(axis=1)
                                        x_membership = x_potential/x_density[:,np.newaxis]
                                        reg_Y = (x_membership*Y[classes==c]).sum(axis=1)

                                    # regression_linewidth = 2*(linewidth+1)+np.sqrt(markersize)/np.pi
                                    regression_linewidth = (linewidth+1)
                                    simple_plot(figure,reg_X[True-np.isnan(reg_Y)],reg_Y[True-np.isnan(reg_Y)],class_color,xlabel=xlabel,ylabel=ylabel,linked=True,marker_size=0,linewidth=regression_linewidth,alpha=(alpha+1)/2.)
                                elif world_object['regression'] in ['fireworks']:
                                    class_covariance = np.cov(np.transpose([X[classes==c],Y[classes==c]]),rowvar=False)
                                    vals, vecs = np.linalg.eigh(class_covariance)
                                    order = vals.argsort()[::-1]
                                    vals = vals[order]
                                    vecs = vecs[:,order]
                                    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
                                    width, height = 2.*np.sqrt(2)*np.sqrt(vals)
                                    class_ellipse = patch.Ellipse(xy=class_center, width=width, height=height, angle=theta,color=class_color,alpha=alpha/10., linewidth=2.*linewidth)
                                    figure.gca().add_patch(class_ellipse)
                                    figure.gca().scatter([class_center[0]],[class_center[1]],marker=u'o',c=class_color,s=2.*markersize,linewidth=1.5*linewidth,alpha=(alpha+1.)/2.)
                else:
                    plot_legend = ['colorbar']

                    if labels.dtype == np.dtype('O'):
                        labels = np.sum([i*(labels==l).astype(int) for i,l in enumerate(label_list)],axis=0)+1
                        # print labels
                        label_colors = np.array([np.round(cmap.get_color(l),decimals=3) for l in labels/float(labels.max()+1.)])
                        mpl_norm = mpl.colors.Normalize(vmin=0, vmax=labels.max()+1.)
                        cbar_ticks = np.arange(labels.max())+1
                        cbar_ticklabels = label_list
                    else:
                        label_range = (labels.min() + world_object['label_range'][0]*(labels.max()-labels.min())/100.,labels.min() + world_object['label_range'][1]*(labels.max()-labels.min())/100.)
                        label_colors = np.array([np.round(cmap.get_color(l),decimals=3) for l in (labels-label_range[0])/float(label_range[1]-label_range[0])])
                        mpl_norm = mpl.colors.Normalize(vmin=label_range[0], vmax=label_range[1])
                        cbar_ticks = np.linspace(label_range[0],label_range[1],5)
                        cbar_ticklabels = np.round(cbar_ticks,decimals=3)

                        # print (labels.min(),labels.max()), world_object['label_range']
                        # print label_range, cbar_ticks, cbar_ticklabels

                    color_dict = dict(red=[],green=[],blue=[])
                    for p in np.sort(cmap._color_points.keys()):
                        for k,c in enumerate(['red','green','blue']):
                            color_dict[c] += [(p,cmap._color_points[p][k],cmap._color_points[p][k])]
                    for c in ['red','green','blue']:
                        color_dict[c] = tuple(color_dict[c])
                    # print color_dict

                    mpl_cmap = mpl.colors.LinearSegmentedColormap(cmap.name, color_dict)

                    if world_object['plot'] == 'scatter':
                        simple_plot(figure,X,Y,label_colors,xlabel=xlabel,ylabel=ylabel,linked=False,marker_size=markersize*sizes/sizes.mean(),linewidth=linewidth,alpha=alpha)
                    elif world_object['plot'] == 'line':
                        simple_plot(figure,np.sort(X),Y[np.argsort(X)],label_colors,xlabel=xlabel,ylabel=ylabel,linked=True,marker_size=0,linewidth=(linewidth+1)/2.,alpha=alpha)
                    elif world_object['plot'] == 'map':
                        #data_range = [[X.min()-0.5*(X.max()-X.min()),X.max()+0.5*(X.max()-X.min())],[Y.min()-0.5*(Y.max()-Y.min()),Y.max()+0.5*(Y.max()-Y.min())]]
                        data_range = [[all_X.min()-0.5*(all_X.max()-all_X.min()),all_X.max()+0.5*(all_X.max()-all_X.min())],[all_Y.min()-0.5*(all_Y.max()-all_Y.min()),all_Y.max()+0.5*(all_Y.max()-all_Y.min())]]
                        # print mpl_norm(labels)
                        map_plot(figure,X,Y,mpl_norm(labels),mpl_cmap,data_range,xlabel=xlabel,ylabel=ylabel,smooth_factor=smooth*10,n_points=10.*world_object['n_points'],linewidth=linewidth,alpha=alpha)
                    elif world_object['plot'] == 'histogram' and label_variable == X_variable:
                        x_range = (all_X.min() + world_object['X_range'][0]*(all_X.max()-all_X.min())/100.,all_X.min() + world_object['X_range'][1]*(all_X.max()-all_X.min())/100.)
                        x_bins = np.linspace(x_range[0],x_range[1],n_slices+1)
                        bin_colors = np.array([np.round(cmap.get_color(b),decimals=3) for b in (x_bins-label_range[0])/float(label_range[1]-label_range[0])])
                        
                        bar_histo_plot(figure,X,bin_colors,xlabel=xlabel,ylabel="Number of Elements",n_bins=n_slices,X_range=x_range,linewidth=linewidth,alpha=alpha)
                            


                if world_object['plot'] in ['boxplot','violin','PCA']:
                    class_labels = [np.unique(labels[classes==c])[0] for c in class_list]
                    class_i_labels = [np.where(label_list==l)[0][0] for l in class_labels]
                    if world_object['label_colormap']['name'] != 'glasbey':
                        # class_colors = [np.array(cmap.get_color((i_class+1)/float(n_classes+1))) for i_class in xrange(n_classes)]
                        class_colors = [np.round(cmap.get_color((i_label+1)/float(n_labels+1)),decimals=3) for i_label in class_i_labels]
                    else:
                        # class_colors = [np.array(cmap.get_color((i_class+1)/255.)) for i_class in xrange(n_classes)]
                        class_colors = [np.round(cmap.get_color((i_label+1)/255.),decimals=3) for i_label in class_i_labels]

                    if world_object['plot'] in ['boxplot','violin']:

                        if Y_variable == "":
                            if world_object['plot'] == 'violin':
                                violin_plot(figure,np.arange(n_classes),np.array([X[classes==c] for c in class_list]),class_colors,xlabel="",ylabel=xlabel,linewidth=linewidth,marker_size=(markersize+50)/10)
                            elif world_object['plot'] == 'boxplot':
                                box_plot(figure,np.arange(n_classes),np.array([X[classes==c] for c in class_list]),class_colors,xlabel="",ylabel=xlabel,linewidth=linewidth,marker_size=(markersize+50)/10,alpha=alpha)
                        else:

                            plot_legend = ['labels']

                            n_slices = world_object['n_points']
                            if X.dtype != np.dtype('O') and len(np.unique(X)) > len(X)/5. or (len(np.unique(X)) > 50):
                                magnitude = np.power(10,np.around(4*np.log10(np.nanmean(X)+np.nanstd(X)+1e-7))/4+0.5)
                                magnitude = np.around(magnitude,int(-np.log10(magnitude))+1)

                                print magnitude,"(",n_slices,")  -> ",np.arange(-n_slices-1,n_slices+1)*magnitude

                                print X
                                class_values = np.array(np.maximum(np.minimum(np.around(n_slices*X/magnitude),n_slices),-n_slices),int)
                                
                                print class_values
                                X = (np.arange(-n_slices-1,n_slices+1)*magnitude/float(n_slices))[class_values+n_slices+1]


                            X_values = np.sort(np.unique(X))

                            class_median = {}
                            for c in class_list:
                                class_median[c] = []

                            for i,x in enumerate(X_values):
                                if n_classes > 1:
                                    x_classes = i + np.arange(n_classes)/(3.*float(n_classes-1)) - 1/(2*3.)
                                    width = 1./(9.*float(n_classes-1))
                                else:
                                    x_classes = np.array([i])
                                    width = None

                                y_data = np.array([Y[X==x][classes[X==x]==c] for c in class_list if (classes[X==x]==c).sum()>1])
                                x_classes = np.array([x_classes[i] for i,c in enumerate(class_list) if (classes[X==x]==c).sum()>1])
                                x_class_colors = np.array([class_colors[i] for i,c in enumerate(class_list) if (classes[X==x]==c).sum()>1])
                                x_represented_classes = np.array([c for c in class_list if (classes[X==x]==c).sum()>1])

                                # print y_data

                                if len(x_classes)>0:

                                    if world_object['plot'] == 'violin':
                                        violin_plot(figure,x_classes,y_data,x_class_colors,xlabel=xlabel,ylabel=ylabel,violin_width=width,linewidth=linewidth,marker_size=(markersize+50)/10)
                                    else:
                                        box_plot(figure,x_classes,y_data,x_class_colors,xlabel=xlabel,ylabel=ylabel,box_width=width,linewidth=linewidth,marker_size=markersize,alpha=alpha)
                       
                        
                                for c,x,y in zip(x_represented_classes,x_classes,y_data):
                                    class_median[c] += [[x,np.percentile(y,50)]]

                            for c in class_list:
                                label = np.unique(labels[classes==c])[0]
                                i_label = np.where(label_list==label)[0][0]

                                if world_object['label_colormap']['name'] != 'glasbey':
                                    # class_color = np.array(cmap.get_color((i_class+1)/float(n_classes+1)))
                                    class_color = np.round(cmap.get_color((i_label+1)/float(n_labels+1)),decimals=3)
                                else:
                                    #class_color = np.array(cmap.get_color((i_class+1)/255.))
                                    class_color = np.round(cmap.get_color((i_label+1)/255.),decimals=3)

                                _,figure_labels = figure.gca().get_legend_handles_labels()
                                plot_label = str(label) if str(label) not in figure_labels else None

                                if len(class_median[c])>0:
                                    simple_plot(figure,np.array(class_median[c])[:,0],np.array(class_median[c])[:,1],class_color,xlabel=xlabel,ylabel=ylabel,linked=True,marker_size=0,linewidth=(linewidth+1)/2.,alpha=alpha/2,label=plot_label)
                        


                    elif world_object['plot'] == 'PCA':
                        plot_legend = ['labels']
                        pca, projected_data = pca_analysis(data,classes,class_colors,class_labels,variables,pca_figure=figure,linewidth=linewidth,marker_size=markersize,alpha=alpha,draw_classes=world_object['regression'] in ['fireworks'])  
                        x_min,x_max = np.percentile(projected_data[:,0],1)-4,np.percentile(projected_data[:,0],99)+4
                        y_min,y_max = np.percentile(projected_data[:,1],1)-4,np.percentile(projected_data[:,1],99)+4    
                        figure.gca().set_xlim(x_min + world_object['X_range'][0]*(x_max-x_min)/100.,x_min + world_object['X_range'][1]*(x_max-x_min)/100.)
                        figure.gca().set_xticklabels(figure.gca().get_xticks())
                        figure.gca().set_ylim(y_min + world_object['Y_range'][0]*(y_max-y_min)/100.,y_min + world_object['Y_range'][1]*(y_max-y_min)/100.)
                        figure.gca().set_yticklabels(figure.gca().get_yticks())        


                if world_object['plot'] in ['boxplot','violin']:

                    ticks = figure.gca().get_xticks()

                    if Y_variable == "":
                        ticklabels = class_labels
                        ticks = range(n_classes)
                    else:
                        x_range = (-1 + world_object['X_range'][0]*(len(X_values))/100.,-1 + world_object['X_range'][1]*(len(X_values))/100.)
                        figure.gca().set_xlim(*x_range)
                        ticks = [i for i in range(len(X_values)) if (i>=x_range[0]) and (i<=x_range[1])]
                        ticklabels = [X_values[i] for i in ticks]
                        print ticks,' -> ',ticklabels
                        # print X_values
                    # ticklabels = []
                    # for t in xrange(len(ticks)):
                        # if np.any(np.isclose(ticks[t],range(n_classes),1e-3)):
                            # for c in xrange(n_classes):
                                # if np.isclose(ticks[t],c,1e-3):
                                    # ticklabels += [class_labels[c]]
                        # else:
                            # ticklabels += [""]
                    figure.gca().set_xticklabels(ticklabels)
                    figure.gca().set_xticks(ticks)

                # elif world_object['plot'] != 'PCA':
                #     x_range = (X.min() + world_object['X_range'][0]*(X.max()-X.min())/100.,X.min() + world_object['X_range'][1]*(X.max()-X.min())/100.)
                #     figure.gca().set_xlim(*x_range)
                #     figure.gca().set_xticklabels(figure.gca().get_xticks())

                # if world_object['plot'] in ['distribution','cumulative']: 
                #     figure.gca().set_ylim(*world_object['Y_range'])
                #     figure.gca().set_yticklabels(figure.gca().get_yticks())
                # elif world_object['plot'] in ['boxplot','violin'] and  Y_variable == "":
                #         x_range = (X.min() + world_object['X_range'][0]*(X.max()-X.min())/100.,X.min() + world_object['X_range'][1]*(X.max()-X.min())/100.)
                #         figure.gca().set_ylim(*x_range)
                #         figure.gca().set_yticklabels(figure.gca().get_yticks())
                # elif world_object['plot'] != 'PCA':
                #     y_range = (Y.min() + world_object['Y_range'][0]*(Y.max()-Y.min())/100.,Y.min() + world_object['Y_range'][1]*(Y.max()-Y.min())/100.)
                #     figure.gca().set_ylim(*y_range)
                #     figure.gca().set_yticklabels(figure.gca().get_yticks())


                legend_locations = dict(top_right=1,bottom_right=4,top_left=2,bottom_left=3)
                colorbar_locations = dict(top_right=[0.12, 0.98, 0.78, 0.02],bottom_right=[0.9, 0.1, 0.02, 0.8],top_left=[0.0, 0.1, 0.02, 0.8],bottom_left=[0.12, 0.05, 0.78, 0.02])
                colorbar_orientations = dict(top_right='horizontal',bottom_right='vertical',top_left='vertical',bottom_left='horizontal')
                self._plot_legend[world_object.name] = []
                self._colorbar_axes[world_object.name] = None
                if world_object['legend'] != "":
                    self._plot_legend[world_object.name] = plot_legend
                    if 'labels' in plot_legend:
                        figure.gca().legend(loc=legend_locations[world_object['legend']])
                    if 'colorbar' in plot_legend:
                        ax = figure.gca()
                        cb_ax = figure.add_axes(colorbar_locations[world_object['legend']])
                        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=mpl_cmap, norm=mpl_norm, orientation=colorbar_orientations[world_object['legend']])
                        cb.set_label(label_label)

                        cb.set_ticks(cbar_ticks)
                        cb.set_ticklabels(cbar_ticklabels)

                        self._colorbar_axes[world_object.name] = cb_ax
                        figure.sca(ax)
                
                self.update_plot_limits(world_object)
                self.update_figure_title(world_object)
                # time.sleep(0.05)

    def update_plot_limits(self, world_object):
        if world_object:
            # print "Display figure!"
            figure = plt.figure(self._figures[world_object.name])

            variables = [v for v in world_object.data.keys() if not 'Unnamed:' in v]
            variables = [v for v in variables if not np.array(world_object.data[v]).dtype == np.dtype('O')]
            data = np.transpose([world_object.data[variable] for variable in variables])

            X_variable = world_object['X_variable']
            if X_variable != "":
                X = np.array(world_object.data[X_variable])
            else:
                X = np.array(world_object.data.index)

            Y_variable = world_object['Y_variable']
            if Y_variable != "":
                Y = np.array(world_object.data[Y_variable])
            else:
                Y = np.array(world_object.data.index)

            class_variable = world_object['class_variable']
            if class_variable != "":
                classes = np.array(world_object.data[class_variable])
            else:
                classes = np.ones_like(world_object.data.index)

            label_variable = world_object['label_variable']
            if label_variable != "":
                labels = np.array(world_object.data[label_variable])
            else:
                labels = classes

            if world_object['plot'] != 'PCA':
                point_invalidity = np.zeros(len(classes)).astype(bool)
                point_invalidity = point_invalidity | np.isnan(X)
                point_invalidity = point_invalidity | np.isnan(Y)
                if classes.dtype != np.dtype('O'):
                    point_invalidity = point_invalidity | np.isnan(classes)
                if labels.dtype != np.dtype('O'):
                    point_invalidity = point_invalidity | np.isnan(labels)
                # print "Invalid points : ",point_invalidity.sum()
                valid_points = np.where(True - point_invalidity)
            else:
                point_invalidity = np.zeros(len(classes)).astype(bool)
                for variable in variables:
                    if np.array(world_object.data[variable]).dtype != np.dtype('O'):
                        point_invalidity = point_invalidity | np.isnan(np.array(world_object.data[variable]))
                valid_points = np.where(True - point_invalidity)

            data = data[valid_points]
            X = X[valid_points]
            Y = Y[valid_points]
            classes = classes[valid_points]
            labels = labels[valid_points]

            n_rows, n_cols = figure.get_axes()[0].get_subplotspec().get_gridspec().get_geometry()
            plot_legend = self._plot_legend[world_object.name]
            if 'colorbar' in plot_legend:
                n_cols += 1

            for i_ax, axes in enumerate(figure.get_axes()):
                if not axes == self._colorbar_axes[world_object.name]:
                    if world_object['plot'] != 'PCA':
                        if world_object['plot'] in ['boxplot','violin']:
                            if Y_variable == "":
                                pass
                            else:
                                n_slices = world_object['n_points']
                                if X.dtype != np.dtype('O') and len(np.unique(X)) > len(X)/5. or (len(np.unique(X)) > 50):
                                    magnitude = np.power(10,np.around(4*np.log10(np.nanmean(X)+np.nanstd(X)+1e-7))/4+0.5)
                                    magnitude = np.around(magnitude,int(-np.log10(magnitude))+1)
                                    class_values = np.array(np.maximum(np.minimum(np.around(n_slices*X/magnitude),n_slices),-n_slices),int)
                                    X = (np.arange(-n_slices-1,n_slices+1)*magnitude/float(n_slices))[class_values+n_slices+1]


                                X_values = np.sort(np.unique(X))

                                x_range = (-1 + world_object['X_range'][0]*(len(X_values))/100.,-1 + world_object['X_range'][1]*(len(X_values))/100.)
                                figure.gca().set_xlim(*x_range)
                                ticks = [i for i in range(len(X_values)) if (i>=x_range[0]) and (i<=x_range[1])]
                                ticklabels = [X_values[i] for i in ticks]

                        else:
                            x_range = (X.min() + world_object['X_range'][0]*(X.max()-X.min())/100.,X.min() + world_object['X_range'][1]*(X.max()-X.min())/100.)
                            axes.set_xlim(*x_range)
                            
                            if i_ax/n_cols == n_rows-1:
                                axes.set_xticklabels(axes.get_xticks())
                            else:
                                axes.set_xticklabels([])
                                axes.set_xlabel("")

                    if world_object['plot'] in ['distribution','cumulative']: 
                        axes.set_ylim(*world_object['Y_range'])
                        axes.set_yticklabels(axes.get_yticks())
                    elif world_object['plot'] in ['histgram']:
                        y_range = (world_object['Y_range'][0]*len(all_X)/100.,world_object['Y_range'][1]*len(all_X)/100.)
                        axes.set_ylim(*y_range)
                        axes.set_yticklabels(axes.get_yticks())
                    elif world_object['plot'] in ['boxplot','violin'] and  Y_variable == "":
                            x_range = (X.min() + world_object['X_range'][0]*(X.max()-X.min())/100.,X.min() + world_object['X_range'][1]*(X.max()-X.min())/100.)
                            axes.set_ylim(*x_range)
                            axes.set_yticklabels(axes.get_yticks())
                    elif world_object['plot'] != 'PCA':
                        y_range = (Y.min() + world_object['Y_range'][0]*(Y.max()-Y.min())/100.,Y.min() + world_object['Y_range'][1]*(Y.max()-Y.min())/100.)
                        axes.set_ylim(*y_range)

                        if i_ax%n_cols == 0:
                            axes.set_yticklabels(axes.get_yticks())
                        else:
                            axes.set_yticklabels([])
                            axes.set_ylabel("")

            plt.draw()

    def update_figure_title(self, world_object):
        print world_object['title']
        if world_object:
            # print "Display figure!"
            figure = plt.figure(self._figures[world_object.name])
            title = world_object['title']
            figure.set_label(title)
            figure.suptitle(title,fontsize=18)
            plt.draw()


    def get_matching_data(self, figure, axes, point):
        world_object_names = [w for w in self._figures.keys() if self._figures[w] == figure.number]
        if len(world_object_names)>0:
            world_object = self.world[world_object_names[0]]

            variables = [v for v in world_object.data.keys() if not 'Unnamed:' in v]
            # variables = [v for v in variables if not np.array(world_object.data[v]).dtype == np.dtype('O')]
            # data = np.transpose([world_object.data[variable] for variable in variables])

            front_variables = []


            subplot_variable = world_object['subplot_variable']
            if subplot_variable != "":
                subplots = np.array(world_object.data[subplot_variable])
                if not subplot_variable in front_variables:
                    front_variables += [subplot_variable]
            else:
                subplots = np.ones_like(world_object.data.index)

            class_variable = world_object['class_variable']
            if class_variable != "":
                if not class_variable in front_variables:
                    front_variables += [class_variable]

            label_variable = world_object['label_variable']
            if label_variable != "":
                if not label_variable in front_variables:
                    front_variables += [label_variable]

            size_variable = world_object['size_variable']
            if size_variable != "":
                if not size_variable in front_variables:
                    front_variables += [size_variable]

            X_variable = world_object['X_variable']
            if X_variable != "":
                X = np.array(world_object.data[X_variable])
                if not X_variable in front_variables:
                    front_variables += [X_variable]
            else:
                X = np.array(world_object.data.index)

            Y_variable = world_object['Y_variable']
            if Y_variable != "":
                Y = np.array(world_object.data[Y_variable])
                if not Y_variable in front_variables:
                    front_variables += [Y_variable]
            else:
                Y = np.array(world_object.data.index)

            for v in front_variables:
                if v in variables:
                    variables.remove(v)
            variables = front_variables + variables
            
            # if len(variables)>10:
                # variables = variables[:10]

            if subplot_variable != "":
                axes_subplot = axes.get_title()
                data_index = world_object.data.index[np.where((np.array([str(s) for s in subplots])==axes_subplot)&(X==point[0])&(Y==point[1]))]
            else:
                data_index = world_object.data.index[np.where((X==point[0])&(Y==point[1]))]

            if len(data_index)>0:
                # data_row = world_object.data.iloc[data_index[0]]
                data_row = world_object.data.loc[data_index[0]]
                # print ""
                # print data_row[variables]
                # print ""
                w = DataPointWidget(data_row[variables])
                dialog = ModalDialog(w)
                dialog.setWindowTitle(world_object.name)
                dialog.exec_()



class DataPointWidget(QtGui.QWidget):

    def __init__(self, data, parent=None):
        QtGui.QWidget.__init__(self, parent=parent)

        self._layout = QtGui.QGridLayout(self)

        self.max_variables = 6

        self.variables = data.axes[0]
        self.extra_variables = []
        self.data = data

        if len(self.variables)>self.max_variables:
            self.extra_variables = self.variables[self.max_variables:]
            self.variables = self.variables[:self.max_variables]

        for i,v in enumerate(self.variables):
            # self.data_labels[v] = QtGui.QLabel("<b>"+str(v)+"</b> : "+str(data[v]))
            self._layout.addWidget(QtGui.QLabel("<b>"+str(v)+"</b>"),i,0)
            self._layout.addWidget(QtGui.QLabel(str(self.data[v])),i,1) 

        self._variable_cbox = None
        self._variable_label = None

        if len(self.extra_variables)>0:
            self._variable_cbox = QtGui.QComboBox(self)
            for v in self.extra_variables:
                self._variable_cbox.addItem(v)
            self._variable_cbox.currentIndexChanged.connect(self._extra_variable_changed)
            self._variable_label = QtGui.QLabel(str(self.data[self._variable_cbox.currentText()]))
            self._layout.addWidget(self._variable_cbox,len(self.variables),0)
            self._layout.addWidget(self._variable_label,len(self.variables),1)

    def _extra_variable_changed(self, idx):
        self._variable_label.setText(str(self.data[self._variable_cbox.currentText()]))








