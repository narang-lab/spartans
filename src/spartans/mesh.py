import numpy as np
from spartans.mesh_utils import tetrahedra_volumes, triangle_surface_areas

class Mesh:
    '''
    A Mesh is an object representing properties of a given mesh
    identifier. It queries the Database and assigns properties as methods.
    '''

    def __init__(self, database_dict, mesh_identifier):
        '''
        Each mesh has a unique identifier of the form mesh({mesh_identifier})
        '''
        self.mesh_identifier = mesh_identifier
        
        # Required mesh properties
        self.vertices=database_dict[self.mesh_identifier]['vertices']
        
        self.tetrahedra_indices=database_dict[self.mesh_identifier]['tetrahedra_indices']
        self.tetrahedra_volumes = tetrahedra_volumes(self.vertices,self.tetrahedra_indices)

        self.surface_triangle_indices = database_dict[self.mesh_identifier]['triangle_indices']
        self.triangle_surface_areas = triangle_surface_areas(self.vertices,self.surface_triangle_indices)
        
        # Optional mesh property
        self.surface_triangles_normals = database_dict[self.mesh_identifier].get('surface_normals',None)
        if self.surface_triangles_normals is None:
            self.surface_triangles_normals = construct_surface_normals(self.vertices,self.surface_triangle_indices)

    def __str__(self):
        '''
        Prints mesh identifier
        '''
        return self.mesh_identifier

