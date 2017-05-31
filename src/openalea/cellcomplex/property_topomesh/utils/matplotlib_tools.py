import numpy as np

from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property, is_triangular

from matplotlib import cm
from matplotlib import tri
from matplotlib.colors import Normalize
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt

def mpl_draw_topomesh(topomesh,figure,degree=2,coef=1,property_name="",colormap='viridis',color='k',alpha=1.0,cell_edges=False):

    positions = topomesh.wisp_property('barycenter',0)

    compute_topomesh_property(topomesh,'vertices',2)
    triangles = topomesh.wisp_property('vertices',2).values()
    
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
        else:
            triangle_property = topomesh.wisp_property(property_name,2).values()
            if triangle_property.ndim == 1:
                mpl_colormap = cm.ScalarMappable(norm=Normalize(vmin=triangle_property.min(), vmax=triangle_property.max()),cmap=cm.cmap_d[colormap])
                colors = mpl_colormap.to_rgba(triangle_property)[:,:3]
                print colors
        figure.gca().add_collection(PolyCollection(triangle_positions[:,:2].reshape((len(triangles),3,2)),facecolors=colors,cmap=colormap,linewidth=0.5,alpha=0.8*alpha))
        #figure.gca().tripcolor(triangulation,facecolors=colors,alpha=0.8)
    elif degree==1:
        if is_triangular(topomesh) and not cell_edges:
            figure.gca().triplot(triangulation,color=color,linewidth=2,alpha=alpha)
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
                figure.gca().plot(p[:,0],p[:,1],color=color,linewidth=2,alpha=alpha)

    elif degree==0:
        if topomesh.has_wisp_property('color',0):
            colors = topomesh.wisp_property('color',0).values()
        else:
            colors = color
        figure.gca().scatter(positions.values()[:,0],positions.values()[:,1],s=20,edgecolor=color,color=colors,alpha=alpha)

    figure.canvas.draw()
    plt.pause(1e-3)
