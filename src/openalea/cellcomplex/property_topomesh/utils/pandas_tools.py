import numpy as np
import pandas as pd

def topomesh_to_dataframe(topomesh, degree=3, properties=None):
    
    if properties is None:
        properties = topomesh.wisp_properties(degree).keys()
    
    dataframe = pd.DataFrame()
    dataframe['id'] = np.array(list(topomesh.wisps(degree)))

    for property_name in properties:
        if np.array(topomesh.wisp_property(property_name,degree).values()[0]).ndim == 0:
            print "  --> Adding column ",property_name
            dataframe[property_name] = topomesh.wisp_property(property_name,degree).values(dataframe['id'].values)
    
    if topomesh.has_wisp_property('barycenter',degree):
        for i,k in enumerate(['x','y','z']):
             dataframe['center_'+k] = topomesh.wisp_property('barycenter',degree).values(dataframe['id'].values)[:,i]

    dataframe = dataframe.set_index('id')
    dataframe.index.name = None
        
    return dataframe

def normalized_data_range(data_range, dataframe, property_name):
    return tuple([100.*(s-dataframe[property_name].min())/(dataframe[property_name].max()-dataframe[property_name].min()) for s in data_range])
        


