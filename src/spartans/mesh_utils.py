import numpy as np

from spartans.interface import Interface

'''
Mesh methods mainly used to construct Interface objects
'''

def tetrahedra_volumes(vertices, tetrahedra_indices):
    '''
    Computes tetrahedra volumes
    from https://mathworld.wolfram.com/Tetrahedron.html
    '''
    tetrahedra_vertices = vertices[tetrahedra_indices]
    num_tetrahedra,num_vertices_per_tet, dim = tetrahedra_vertices.shape

    # Jacobian matrix for tetrahedra volumes
    jac = np.zeros((num_tetrahedra, dim+1,dim+1))
    jac[:, :, 0] = np.ones(dim+1)
    jac[:, :, 1:] = tetrahedra_vertices
    jac = np.swapaxes(jac, 1, 2)
    vols = np.abs(np.linalg.det(jac))/6

    return vols

def triangle_surface_areas(vertices,triangle_indices):
    '''
    Computes triangle surface areas for 3D-embedded triangles
    from https://mathworld.wolfram.com/TriangleArea.html
    '''
    surface_triangle_coords   = vertices[triangle_indices]
    num_triangles,num_vertices_per_triangle, dim = surface_triangle_coords.shape

    # Jacobian matrices for surface areas
    jac1 = np.zeros((num_triangles, dim , dim ))
    jac1[:, :, 0 ]  = np.ones(dim)
    jac1[:, :, 1:]  = surface_triangle_coords[:,:,[0,1]]
    detJac1 = np.abs(np.linalg.det(jac1))**2
    
    jac2 = np.zeros((num_triangles, dim, dim ))
    jac2[:, :, 0 ]  = np.ones(dim)
    jac2[:, :, 1:]  = surface_triangle_coords[:,:,[1,2]]
    detJac2 = np.abs(np.linalg.det(jac2))**2

    jac3 = np.zeros((num_triangles, dim, dim ))
    jac3[:, :, 0 ]  = np.ones(dim)
    jac3[:, :, 1:]  = surface_triangle_coords[:,:,[2,0]]
    detJac3 = np.abs(np.linalg.det(jac3))**2

    surface_areas = np.sqrt(detJac1 + detJac2 + detJac3)/2

    return surface_areas

def body_volume_integral(body_array, tetrahedra_indices, tetrahedra_volumes):
    '''
    Computes volume integral of input body_array
    '''

    dist  = body_array.copy()

    dist  = dist[:,tetrahedra_indices]
    dist *= tetrahedra_volumes[:,np.newaxis]
    dist /= 4

    return np.sum(dist)

def surface_area_integral(surface_array, triangle_surface_areas):
    '''
    Computes surface area integral of input surface_array
    '''
    dist  = surface_array.copy()

    dist *= triangle_surface_areas[:,np.newaxis]
    dist /= 3

    return np.sum(dist)

def construct_inside_outside_mask(surface_normals,velocities):
    '''
    Constructs a mask of shape (n_states, n_triangles) which is 1 when the given
    state is pointing inward (into the geometry on the given triangle)
    and 0 when it is pointing outward (out of the geometry).

    Inputs:
        velocities - (n_states, 3)
        surface_normals - (n_triangles, 3), point outwards from the geometry

    '''
    cos_theta = np.tensordot(velocities,surface_normals,axes=[1,1]) # Surface normals are outward
    mask = np.zeros_like(cos_theta)
    mask[cos_theta < 0] = 1 # 1 when pointing inward, 0 otherwise
    return mask

def compute_surface_normals(tetrahedra, triangles, vertices):
    '''
    Calculates the signed surface normals (positive points outward) of each triangle.
    The signs are determined by finding the corresponding tetrahedraon the triangle
    is a part of and checking if the normal points into that tetrahedron or out of it.

    Inputs:
        tetrahedra - (n_tetrahedra, 4) Indexes into vertices
        triangles - (n_triangles, 3) Indexes into vertices
        vertices - (n_vertices, 3) Locations of the vertices
    '''
    # Calculate the surface normal for each triangle
    normals = np.cross(vertices[triangles[:,0]] - vertices[triangles[:,1]], vertices[triangles[:,0]] - vertices[triangles[:,2]])

    # The rest of this routine is making sure the sign is right, so that the normals point outwards.

    # Construct the four triangles on each tetrahedron, has shape (4, n_tets, 3)
    parts = [(0,1,2),(0,1,3),(1,2,3),(0,2,3)]
    tet_tris = {tuple(sorted(tetrahedra[i,parts[k]])):i for i in range(len(tetrahedra)) for k in range(len(parts))}

    # Sort the indices in triangles
    tri_sorted = np.sort(triangles, axis=1)

    # Find the tetrahedron corresponding to each triangle
    tets = list(tet_tris[tuple(tri)] for tri in tri_sorted)

    for i in range(len(triangles)):
        tet = tetrahedra[tets[i]]
        tri = triangles[i]
        apex = vertices[list(set(tet) - set(tri))[0]]

        vert = vertices[tri[0]] # Can be any vertex on the triangle

        sign = np.dot(vert - apex, normals[i])
        if sign < 0:
            normals[i] *= -1

    return normals


def map_triangles_from_cpp_to_user(surface_indexing_cpp, surface_indexing_usr, array_to_map):
    '''
    Inputs:
        surface_indexing_cpp - Integer array of shape (n_triangles, 3) specifying the vertices of the triangles in cpp's ordering.
        surface_indexing_usr - Integer array of shape (n_triangles, 3) specifying the vertices of the triangles in the user's ordering.
        array_to_map - Array of shape (n_states,n_triangles, 3) to be remapped from the cpp ordering to the user's ordering.

    Assumptions:
        All entries in surface_indexing_cpp appear in surface_indexing_user, and vice-versa.
    '''

    # First, we sort each triangle in cpp, and each triangle in usr, recording the indexing we used for that.
    vert_inds_usr = np.argsort(surface_indexing_usr, axis=1)
    vert_inds_cpp = np.argsort(surface_indexing_cpp, axis=1)

    usr_vert_sorted = np.take_along_axis(surface_indexing_usr, vert_inds_usr, axis=1)
    cpp_vert_sorted = np.take_along_axis(surface_indexing_cpp, vert_inds_cpp, axis=1)

    # Now we use the fun property that if inds(x) is an array of indices that
    # sorts x, then argsort(inds(x)) returns the indexing that puts the sorted
    # array back into the original order.
    # We know that cpp and usr are the same after lex-sorting, so
    # to get the mapping from cpp -> usr, we compute lexsort(cpp),
    # which sorts cpp, then argsort(lexsort(usr)), which runs from the
    # sorted ordering to the usr ordering, and then apply those in sequence.

    inds_cpp = np.lexsort(cpp_vert_sorted.T) # We have to transpose because lexsort wants the tuples it sorts on to span the first axis.
    inds_usr = np.argsort(np.lexsort(usr_vert_sorted.T)) # Same as above.

    # If we do cpp[inds_cpp,:][inds_usr,:], this sorts cpp to have
    # the same triangle ordering as usr. And we can similarly reorder
    # the input array.

    # All we need to do now is find out how to reorder the vertices
    # of each triangle to match between cpp and usr.
    # We do this with the same trick as before:

    # First we reorder the 'triangle' axis of the vertex orderings from the beginning:
    vert_inds_cpp = vert_inds_cpp[inds_cpp,:][inds_usr,:]

    # Now vert_inds_cpp has the same triangle ordering as vert_inds_usr, so we can pull the 
    # same trick as before with argsort(inds):
    vert_inds_usr_sort = np.argsort(vert_inds_usr, axis=1)

    # So now we can sort cpp to match usr with
    # cpp_tri_sorted = cpp[inds_cpp,:][inds_usr,:]
    # cpp_sorted = np.take_along_axis(np.take_along_axis(cpp_tri_sorted, vert_inds_cpp, axis=1), vert_inds_usr_sort, axis=1)
    #
    # and we can do the same for the input array, taking care to incremenent the axes being sorted
    # because it has an additional 'state' axis:
    array_to_map = array_to_map[:,inds_cpp,:][:,inds_usr,:]
    array_to_map = np.take_along_axis(np.take_along_axis(array_to_map, vert_inds_cpp[np.newaxis,...], axis=2), vert_inds_usr_sort[np.newaxis,...], axis=2)

    return array_to_map

def map_triangles_from_user_to_cpp(surface_indexing_cpp, surface_indexing_usr, array_to_map):
    '''
    Inputs:
        surface_indexing_cpp - Integer array of shape (n_triangles, 3) specifying the vertices of the triangles in cpp's ordering.
        surface_indexing_usr - Integer array of shape (n_triangles, 3) specifying the vertices of the triangles in the user's ordering.
        array_to_map - Array of shape (n_states,n_triangles, 3) to be remapped from the user's ordering to the cpp ordering.

    Assumptions:
        All entries in surface_indexing_cpp appear in surface_indexing_user, and vice-versa.

    '''
    # We use the fact that this is the inverse of the mapping from cpp -> usr, so we just call that mapping function and swap
    # the cpp and usr indices.
    return map_triangles_from_cpp_to_user(surface_indexing_usr, surface_indexing_cpp, array_to_map)


def set_interfaces(config_dict, database_dict):
        '''
        This method reads the user-provided connectivity matrix
        and constructs Interface objects.
        '''

        # Read connectivity ordering list and make dictionary
        connectivity_ordering = database_dict['connectivity_info']['connectivity.ordering']
        connectivity_ordering = {'mesh({})'.format(key): i for i,key in enumerate(connectivity_ordering)}
        
        # Reverse dictionary to serve as lookup table
        ordering_connectivity = {v: k for k, v in connectivity_ordering.items()}

        # Read connectivity matrix
        connectivity_matrix = database_dict['connectivity_info']['connectivity']
        connectivity_rows, connectivity_cols = connectivity_matrix.shape

        # First pass to construct an interface object for each connected pair of meshes
        for row_index in range(connectivity_rows):
                for col_index in range(row_index,connectivity_cols):
                        if connectivity_matrix[row_index,col_index] == 1:
                                # Sort interface meshes since interface(s1,s2) == interface(s2,s1)
                                interface_meshes_identifiers = sorted([ordering_connectivity.get(index) for index in [row_index,col_index]])
                                # Extract structures from meshes
                                interface_structure_identifiers = [database_dict['mesh_to_structure'].get(mesh_id) for mesh_id in interface_meshes_identifiers]
                                interface_structures = [database_dict[structure_id]['structure_object'] for structure_id in interface_structure_identifiers]

                                # Make interface object
                                interface_identifier = 'interface({}--{})'.format(*interface_structure_identifiers)
                                interface = Interface(config_dict, database_dict, interface_structures)
                                
                                # Assign interface object to both structures
                                for structure in interface_structures:
                                        structure.interfaces[interface_identifier] = interface
