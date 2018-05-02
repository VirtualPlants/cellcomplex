import numpy as np

from openalea.cellcomplex.property_topomesh import PropertyTopomesh
from openalea.cellcomplex.property_topomesh.property_topomesh_creation import triangle_topomesh, poly_topomesh
# from openalea.cellcomplex.property_topomesh.property_topomesh_extraction import clean_topomesh_properties
from openalea.cellcomplex.property_topomesh.property_topomesh_analysis import compute_topomesh_property

# from openalea.cellcomplex.property_topomesh.utils.matching_tools import kd_tree_match as match
# from openalea.cellcomplex.property_topomesh.utils.matching_tools import vector_quantization_match as match

from openalea.container import array_dict

from copy import deepcopy
from time import time as current_time

def append_topomesh(input_topomesh, topomesh_to_append, copy=False, properties_to_append=None):
    if copy:
        start_time = current_time()
        topomesh = deepcopy(input_topomesh)
        print "    --> Copying input mesh [",current_time()-start_time,"s]"
    else:
        topomesh = input_topomesh

    if properties_to_append is None:
        properties_to_append = dict([(degree,list(topomesh.wisp_property_names(degree))) for degree in xrange(4)])
    if not 'barycenter' in properties_to_append[0]:
        properties_to_append[0] += ['barycenter']

    start_time = current_time()
    wisp_dict = {}
    for degree in xrange(4):
        wisp_dict[degree] = {}
        for w in topomesh_to_append.wisps(degree):
            new_w = topomesh.add_wisp(degree)
            wisp_dict[degree][w] = new_w

            if degree>0:
                for b in topomesh_to_append.borders(degree,w):
                    topomesh.link(degree,wisp_dict[degree][w],wisp_dict[degree-1][b])
    print "    --> Appending mesh elements [",current_time()-start_time,"s]"

    start_time = current_time()
    for degree in xrange(4):
        for p in properties_to_append[degree]:
            property_dict = topomesh.wisp_property(p,degree)
            if topomesh.has_wisp_property(p,degree) and topomesh_to_append.has_wisp_property(p,degree,is_computed=True):
                property_to_append = topomesh_to_append.wisp_property(p,degree)
                # new_property_to_append = dict(zip([wisp_dict[degree][w] for w in topomesh_to_append.wisps(degree)],[property_to_append[w] for w in topomesh_to_append.wisps(degree)]))
                # new_property = property_dict.to_dict()
                # new_property.update(new_property_to_append)
                new_property = dict(zip(list(property_dict.keys())+[wisp_dict[degree][w] for w in topomesh_to_append.wisps(degree)],list(property_dict.values())+[property_to_append[w] for w in topomesh_to_append.wisps(degree)]))
                topomesh.update_wisp_property(p,degree,new_property)
                # for w in topomesh_to_append.wisps(degree):
                    # print w," (",wisp_dict[degree][w]," ) :",property_to_append[w]
                    # property_dict[wisp_dict[degree][w]] = property_to_append[w]
            else:
                raise KeyError("Property"+p+" not computed on elements of degree "+str(degree)+"!")
    print "    --> Appending mesh properties [",current_time()-start_time,"s]"

    return topomesh, wisp_dict


def merge_topomesh_vertices(topomesh, vertex_to_keep, vertex_to_merge, verbose=False):
    if verbose:
        print "--> Merging vertex ",vertex_to_merge," ->",vertex_to_keep
    for e in topomesh.regions(0,vertex_to_merge):
        topomesh.unlink(1,e,vertex_to_merge)
        topomesh.link(1,e,vertex_to_keep)
    topomesh.remove_wisp(0,vertex_to_merge)


def merge_topomesh_edges(topomesh, edge_to_keep, edge_to_merge, verbose=False):
    if verbose:
        print "--> Merging edge ",edge_to_merge," ->",edge_to_keep
    vertices_to_keep = list(topomesh.borders(1,edge_to_keep))
    vertices_to_merge = list(topomesh.borders(1,edge_to_merge))

    if topomesh.has_wisp_property('barycenter',0):
        # TODO
        # Reorder edge vertices to match vertex positions
        if np.any(topomesh.wisp_property('barycenter',0)[vertices_to_keep[0]] != topomesh.wisp_property('barycenter',0)[vertices_to_merge[0]]):
            vertices_to_merge = vertices_to_merge[::-1]

    for vertex_to_keep, vertex_to_merge in zip(vertices_to_keep,vertices_to_merge):
        if vertex_to_keep != vertex_to_merge:
            merge_topomesh_vertices(topomesh,vertex_to_keep,vertex_to_merge)

    for f in topomesh.regions(1,edge_to_merge):
        topomesh.unlink(2,f,edge_to_merge)
        topomesh.link(2,f,edge_to_keep)

    topomesh.remove_wisp(1,edge_to_merge)


def merge_topomesh_faces(topomesh, face_to_keep, face_to_merge, verbose=False):
    if verbose:
        print "--> Merging face ",face_to_merge," ->",face_to_keep
    edges_to_keep = list(topomesh.borders(2,face_to_keep))
    vertices_to_keep = np.unique(list(topomesh.borders(2,face_to_keep,2)))

    edges_to_merge = list(topomesh.borders(2,face_to_merge))
    vertices_to_merge = np.unique(list(topomesh.borders(2,face_to_merge,2)))

    oriented_edges_to_keep = [edges_to_keep[0]]
    remaining_edges_to_keep = set(edges_to_keep).difference(set(oriented_edges_to_keep))
    oriented_vertices_to_keep = [list(topomesh.borders(1,edges_to_keep[0]))[0]]
    remaining_vertices_to_keep = set(vertices_to_keep)
    while len(remaining_edges_to_keep)>0:
        current_edge = oriented_edges_to_keep[-1]
        previous_vertex, current_vertex = topomesh.borders(1,current_edge)
        if current_vertex in oriented_vertices_to_keep:
            current_vertex, previous_vertex = topomesh.borders(1,current_edge)
        next_vertex = list(set(topomesh.region_neighbors(0,current_vertex)).difference(set([previous_vertex])).intersection(remaining_vertices_to_keep))[0]
        next_edge = list(set(remaining_edges_to_keep).intersection(set(topomesh.regions(0,next_vertex))).intersection(set(topomesh.regions(0,current_vertex))))[0]
        oriented_edges_to_keep += [next_edge]
        remaining_edges_to_keep -= set([next_edge])
        oriented_vertices_to_keep += [current_vertex]
        remaining_vertices_to_keep -= set([current_vertex])

    oriented_edges_to_merge = [edges_to_merge[0]]
    remaining_edges_to_merge = set(edges_to_merge).difference(set(oriented_edges_to_merge))
    oriented_vertices_to_merge = [list(topomesh.borders(1,edges_to_merge[0]))[0]]
    remaining_vertices_to_merge = set(vertices_to_merge)
    while len(remaining_edges_to_merge)>0:
        current_edge = oriented_edges_to_merge[-1]
        previous_vertex, current_vertex = topomesh.borders(1,current_edge)
        if current_vertex in oriented_vertices_to_merge:
            current_vertex, previous_vertex = topomesh.borders(1,current_edge)
        next_vertex = list(set(topomesh.region_neighbors(0,current_vertex)).difference(set([previous_vertex])).intersection(remaining_vertices_to_merge))[0]
        next_edge = list(set(remaining_edges_to_merge).intersection(set(topomesh.regions(0,next_vertex))).intersection(set(topomesh.regions(0,current_vertex))))[0]
        oriented_edges_to_merge += [next_edge]
        remaining_edges_to_merge -= set([next_edge])
        oriented_vertices_to_merge += [current_vertex]
        remaining_vertices_to_merge -= set([current_vertex])

    # print oriented_vertices_to_merge," -> ",oriented_vertices_to_keep
    # print oriented_edges_to_merge," -> ",oriented_edges_to_keep

    if topomesh.has_wisp_property('barycenter',0):
        # TODO
        # Reorder face edges to match vertex positions (need reordering edges first?)
        oriented_vertex_matching = match(topomesh.wisp_property('barycenter',0).values(oriented_vertices_to_keep),topomesh.wisp_property('barycenter',0).values(oriented_vertices_to_merge))
        # print oriented_vertex_matching

        oriented_edge_vertices_to_keep = np.sort([list(topomesh.borders(1,e)) for e in oriented_edges_to_keep])
        oriented_edge_vertices_to_merge = np.sort(np.transpose([oriented_vertices_to_merge,list(oriented_vertices_to_merge[1:])+[oriented_vertices_to_merge[0]]]))
        
        # print oriented_edge_vertices_to_merge," -> ",oriented_edge_vertices_to_keep

        oriented_edge_matching = match(oriented_edge_vertices_to_keep,oriented_edge_vertices_to_merge)
        # print oriented_edge_matching

        oriented_vertices_to_merge = np.array(oriented_vertices_to_merge)[oriented_vertex_matching]
        oriented_edges_to_merge = np.array(oriented_edges_to_merge)[oriented_edge_matching]

    # print oriented_vertices_to_merge," -> ",oriented_vertices_to_keep
    # print oriented_edges_to_merge," -> ",oriented_edges_to_keep

    for edge_to_keep, edge_to_merge in zip(oriented_edges_to_keep,oriented_edges_to_merge):
    # for edge_to_keep, edge_to_merge in zip(edges_to_keep,edges_to_merge):
        if edge_to_keep != edge_to_merge:
            merge_topomesh_edges(topomesh,edge_to_keep,edge_to_merge)

    for c in topomesh.regions(2,face_to_merge):
        topomesh.unlink(3,c,face_to_merge)
        topomesh.link(3,c,face_to_keep)

    topomesh.remove_wisp(2,face_to_merge)


def fuse_topomesh_identical_vertices(topomesh):
    positions = topomesh.wisp_property('barycenter',0)
    points = positions.values(list(topomesh.wisps(0))).astype(np.float16)

    vertex_matching = np.array(list(topomesh.wisps(0)))[match(points,points)]

    vertex_to_fuse = {}
    for v in np.unique(vertex_matching):
        vertex_to_fuse[v] = np.array(list(topomesh.wisps(0)))[vertex_matching==v]
        # print v," -> ",vertex_to_fuse[v]," (",positions.values(vertex_to_fuse[v]),")"

    for vertex_to_keep, vertices_to_fuse in vertex_to_fuse.items():
        for vertex_to_merge in [v for v in vertices_to_fuse if v != vertex_to_keep]:
            merge_topomesh_vertices(topomesh,vertex_to_keep,vertex_to_merge)

    edge_vertices = np.sort([list(topomesh.borders(1,e)) for e in topomesh.wisps(1)])

    edge_matching = np.array(list(topomesh.wisps(1)))[match(edge_vertices,edge_vertices)]

    edge_to_fuse = {}
    for e in np.unique(edge_matching):
        edge_to_fuse[e] = np.array(list(topomesh.wisps(1)))[edge_matching==e]

    for edge_to_keep, edges_to_fuse in edge_to_fuse.items():
        for edge_to_merge in [e for e in edges_to_fuse if e != edge_to_keep]:
            merge_topomesh_edges(topomesh,edge_to_keep,edge_to_merge)

    face_vertices = np.array([np.sort(list(topomesh.borders(2,f,2))) for f in topomesh.wisps(2)])

    if face_vertices.ndim == 2:

        face_matching = np.array(list(topomesh.wisps(2)))[match(face_vertices,face_vertices)]

        face_to_fuse = {}
        for f in np.unique(face_matching):
            face_to_fuse[f] = np.array(list(topomesh.wisps(2)))[face_matching==f]

        for face_to_keep, faces_to_fuse in face_to_fuse.items():
            for face_to_merge in [f for f in faces_to_fuse if f != face_to_keep]:
                merge_topomesh_faces(topomesh,face_to_keep,face_to_merge)

    else:
        face_lengths = np.unique(map(len,face_vertices))

        for l in face_lengths:
            length_face_vertices = np.array([f_v for f_v in face_vertices if len(f_v) == l])
            length_faces = np.array([f for (f,f_v) in zip(topomesh.wisps(2),face_vertices) if len(f_v) == l])

            # print length_face_vertices.shape, len(length_faces)
            face_matching = length_faces[match(length_face_vertices,length_face_vertices)]

            face_to_fuse = {}
            for f in np.unique(face_matching):
                face_to_fuse[f] = length_faces[face_matching==f]

            for face_to_keep, faces_to_fuse in face_to_fuse.items():
                for face_to_merge in [f for f in faces_to_fuse if f != face_to_keep]:
                    merge_topomesh_faces(topomesh,face_to_keep,face_to_merge)

            face_vertices = np.array([np.sort(list(topomesh.borders(2,f,2))) for f in topomesh.wisps(2)])   

    clean_topomesh_properties(topomesh)


def fuse_topomesh_cells(topomesh, cell_to_keep, cell_to_merge):
    faces_to_keep = list(topomesh.borders(3,cell_to_keep))
    faces_to_merge = list(topomesh.borders(3,cell_to_merge))

    for face_to_keep, face_to_merge in zip(faces_to_keep, faces_to_merge):
        topomesh.unlink(3,cell_to_merge,face_to_merge)
        topomesh.link(3,cell_to_keep,face_to_merge)

    topomesh.remove_wisp(3,cell_to_merge)
