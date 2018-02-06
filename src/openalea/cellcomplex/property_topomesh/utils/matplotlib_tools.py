import numpy as np

from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property, is_triangular

from matplotlib import cm
from matplotlib import tri
from matplotlib.colors import Normalize
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt

def mpl_draw_topomesh(topomesh,figure,degree=2,coef=1,property_name="",property_degree=None,colormap='viridis',color='k',alpha=1.0,cell_edges=False,intensity_range=None,linewidth=1,size=20):

    if property_degree is None:
        property_degree = degree

    positions = topomesh.wisp_property('barycenter',0)

    compute_topomesh_property(topomesh,'vertices',2)
    triangles = topomesh.wisp_property('vertices',2).values(list(topomesh.wisps(2)))
    
    compute_topomesh_property(topomesh,'barycenter',2)
    triangle_positions = np.concatenate([b + coef*(p-b) for p,b in zip(positions.values(triangles),topomesh.wisp_property('barycenter',2).values())])
    triangle_triangles = np.arange(3*len(triangles)).reshape((len(triangles),3))
    
    triangulation = tri.Triangulation(triangle_positions[:,0],triangle_positions[:,1],triangle_triangles)
    if degree==2:
        if property_name == "":
            if topomesh.has_wisp_property('color',2):
                colors = topomesh.wisp_property('color',2).values()
            else:
                colors = np.zeros((len(triangles),3))
            figure.gca().add_collection(PolyCollection(triangle_positions[:,:2].reshape((len(triangles),3,2)),facecolors=colors,linewidth=0.5,alpha=0.8*alpha))
        else:
            if property_degree == 2:
                triangle_property = topomesh.wisp_property(property_name,property_degree).values()
            elif property_degree == 0:
                triangle_property = np.concatenate(topomesh.wisp_property(property_name,property_degree).values(triangles))
            if triangle_property.ndim == 1:
                if intensity_range is None:
                    intensity_range = (triangle_property.min(),triangle_property.max())
                # mpl_colormap = cm.ScalarMappable(norm=Normalize(vmin=intensity_range[0], vmax=intensity_range[1]),cmap=cm.cmap_d[colormap])
                # colors = mpl_colormap.to_rgba(triangle_property)[:,:3]
                # print colors
                # figure.gca().add_collection(PolyCollection(triangle_positions[:,:2].reshape((len(triangles),3,2)),facecolors=colors,cmap=colormap,linewidth=0.5,alpha=0.8*alpha))
                if property_degree == 2:
                    figure.gca().tripcolor(triangulation,triangle_property,cmap=colormap,alpha=alpha,vmin=intensity_range[0],vmax=intensity_range[1])
                elif property_degree == 0:
                    figure.gca().tripcolor(triangulation,triangle_property,cmap=colormap,alpha=alpha,vmin=intensity_range[0],vmax=intensity_range[1],shading='gouraud')


    elif degree==1:
        if is_triangular(topomesh) and not cell_edges:
            figure.gca().triplot(triangulation,color=color,linewidth=linewidth,alpha=alpha)
        else:
            compute_topomesh_property(topomesh,'vertices',1)
            if cell_edges:
                boundary_edges = [e for e in topomesh.wisps(1) if topomesh.nb_regions(1,e)==1]
                cell_boundary_edges = [e for e in topomesh.wisps(1) if len(list(topomesh.regions(1,e,2)))>1]
                considered_edges = np.unique(cell_boundary_edges+boundary_edges)
            else:
                considered_edges = list(topomesh.wisps(1))
            edge_points = topomesh.wisp_property('barycenter',0).values(topomesh.wisp_property('vertices',1).values(considered_edges))
            for p in edge_points:
                figure.gca().plot(p[:,0],p[:,1],color=color,linewidth=linewidth,alpha=alpha)

    elif degree==0:
        if property_name == "":
            if topomesh.has_wisp_property('color',0):
                colors = topomesh.wisp_property('color',0).values()
            else:
                colors = color
            figure.gca().scatter(positions.values()[:,0],positions.values()[:,1],s=20,edgecolor=color,color=colors,alpha=alpha)
        else:
            vertex_property = topomesh.wisp_property(property_name,0).values(list(topomesh.wisps(0)))
            if vertex_property.ndim == 1:
                if intensity_range is None:
                    intensity_range = (vertex_property.min(),vertex_property.max())
                figure.gca().scatter(positions.values(list(topomesh.wisps(0)))[:,0],positions.values(list(topomesh.wisps(0)))[:,1],c=vertex_property,s=size,linewidth=linewidth,cmap=colormap,alpha=alpha,vmin=intensity_range[0],vmax=intensity_range[1])

    figure.canvas.draw()
    plt.pause(1e-3)
