import numpy as np
from itertools import product
from collections import defaultdict
import logging

from spartans.logger_builder import setup_logger

class Interface:
    '''
    An Interface is an object representing properties of a given interface
    identifier. It queries the Database. The reason to make separate classes
    is to keep things tidy in the case of multiple interfaces.
    '''

    def __init__(self, config_dict, database_dict, sorted_int_structures):
        '''
        Each interface is defined using two structure objects.
        Note interface(S1,S2) and interface(S2,S1) define the same interface object, 
        and thus we pass Interface a sorted list of structures.
        '''

        self.config_dict = config_dict
        self.database_dict = database_dict

        # Set up interface_identifier
        self.first_structure, self.second_structure = sorted_int_structures
        first_structure_identifier = self.first_structure.structure_identifier
        second_structure_identifier = self.second_structure.structure_identifier
        self.interface_identifier = 'interface({}--{})'.format(first_structure_identifier,second_structure_identifier)
        self.interface_identifier_iosafe = self.first_structure.structure_identifier_iosafe + '__' + self.second_structure.structure_identifier_iosafe
        
        # Set up interface log file
        log_directory = config_dict['directories']['log_directory']
        setup_logger('{}_log'.format(self.interface_identifier_iosafe), '{}/{}_log.txt'.format(log_directory,self.interface_identifier_iosafe))
        self.interface_log = logging.getLogger('{}_log'.format(self.interface_identifier_iosafe))
        self.interface_log.info('Building interface object for {}'.format(self.interface_identifier))
        

        # Set up bounce tensors
        self.bounce_tensors = self.database_dict['connectivity_info']['bounce_tensors']

        # Set up surface mapping dictionaries
        self._set_surface_out_dictionary()
        self._set_triangle_mapping_dictionary()

        # add interface object to database dictionary
        self.database_dict[self.interface_identifier]={'interface_object':self}

    def __str__(self):
        return self.interface_identifier


    def _set_surface_out_dictionary(self):
        '''
        Initialize empty surface_out->surface_in.

        If Interface(S1--S1), i.e. with vacuum, returns only one dictionary:
        - S1_out -> S1_in, i.e. reflection of S1 into S1.

        If Interface(S1--S2), returns all four stacks, i.e.:
        - S1_out -> S1_in, i.e. reflection of S1 into S1.
        - S1_out -> S2_in, i.e. transmission of S1 into S2.
        - S2_out -> S2_in, i.e. reflection of S2 into S2.
        - S2_out -> S1_in, i.e. transmission of S2 into S1.
        '''
        tree = lambda: defaultdict(tree)
        self.surface_out_dictionary = tree()
        structures = (self.first_structure, self.second_structure)
        structure_permutations= set(product(structures,structures))

        for surface_out_structure, surface_in_structure in structure_permutations:
            surface_out_id = surface_out_structure.structure_identifier
            surface_in_id = surface_in_structure.structure_identifier
            self.surface_out_dictionary[surface_out_id][surface_in_id] = {}
    
    def _set_triangle_mapping_dictionary(self):
        '''
        Initialize surface_out->surface_in triangle mappings from user input.
        One mapping dictionary is defined for each bounce tensor, and has the form:
         [
         surface_out_triangle_index, 
         surface_out_triangle_vertex_index,
         surface_in_triangle_index,
         surface_in_triangle_vertex_index
         ]

        Here we group the mappings according to surface_out / surface_in for use later.
        '''
        tree = lambda: defaultdict(tree)
        self.triangle_mapping_dictionary = tree()
        structures = (self.first_structure, self.second_structure)
        structure_permutations= set(product(structures,structures))
        
        for surface_out_structure, surface_in_structure in structure_permutations:
            surface_out_id = surface_out_structure.structure_identifier
            surface_in_id = surface_in_structure.structure_identifier

            connectivity_identifier = 'connectivity({}--{})'.format(surface_out_id,surface_in_id)
            connectivity_dictionary = self.database_dict['connectivity_info'][connectivity_identifier]
            for bounce_id, mapping_dict in connectivity_dictionary.items():
                    self.triangle_mapping_dictionary[surface_out_id][surface_in_id][bounce_id]['extract_from_surface_out'] = mapping_dict[:,[0,1]]
                    self.triangle_mapping_dictionary[surface_out_id][surface_in_id][bounce_id]['deposit_into_surface_in'] =  mapping_dict[:,[2,3]]


    def _interface_needs_to_pull(self,structure_to_pull_from):
        '''
        Checks whether interface needs to pull from specified structure.
        Returns True/False
        '''
        structure_out_id = structure_to_pull_from.structure_identifier

        for surface_in_id, triangle_mapping_dict in self.triangle_mapping_dictionary[structure_out_id].items():
            for bounce_id, mapping_dict in triangle_mapping_dict.items():
                mapping_out = mapping_dict['extract_from_surface_out']
                triangles_out = mapping_out[:,0]
                indices_out   = mapping_out[:,1]
                minimal_surface_out = structure_to_pull_from.surface_out[:,triangles_out,indices_out]

                # check if any surface_out_dictionary is non-zero
                if not minimal_surface_out.any():
                    return True

        return False

    def pull_from(self,structure_to_pull_from):
        '''
        Pulls carriers from specified structure to interface (self).
        Sets carriers into surface_out_dictionary created above.
        '''
        self.interface_log.info('Pulling carriers to {interfaceId} from {structureId}'.format(interfaceId=self.interface_identifier,structureId=structure_to_pull_from.structure_identifier))
        structure_out_id = structure_to_pull_from.structure_identifier

        for surface_in_id, triangle_mapping_dict in self.triangle_mapping_dictionary[structure_out_id].items():
            for bounce_id, mapping_dict in triangle_mapping_dict.items():
                mapping_out = mapping_dict['extract_from_surface_out']
                triangles_out = mapping_out[:,0]
                indices_out   = mapping_out[:,1]
                self.surface_out_dictionary[structure_out_id][surface_in_id][bounce_id] = structure_to_pull_from.surface_out[:,triangles_out,indices_out]

    def _interface_needs_to_push(self,structure_to_push_to):
        '''
        Checks whether interface needs to push to specified structure.
        Returns True/False
        '''
        structure_in_id = structure_to_push_to.structure_identifier
        interface_surface_out_dictionary = self.surface_out_dictionary
        for surface_out_id, surface_in_dict in interface_surface_out_dictionary.items():
            surface_out_dict = surface_in_dict[structure_in_id]

            # If any of them are non-empty, exit loop and return True
            if surface_out_dict:
                return True

        return False

    def push_to(self,structure_to_push_to,weight=1.0):
        '''
        Pushes carriers from interface (self) to specified structure.
        Loops over all surface out dictionaries and bounce tensors and does the following:
        - applies bounce tensor to surface_out
        - ensures carrier-sum before and after bounce is maintained, by applying surface-normals mask
        - swaps bounced and scaled surface_out to appropriate surface_in triangles
        '''
        self.interface_log.info('Pushing carriers from {interfaceId} to {structureId}'.format(interfaceId=self.interface_identifier,structureId=structure_to_push_to.structure_identifier))
        
        structure_in_id = structure_to_push_to.structure_identifier
        material_in_id = structure_to_push_to.material.material_identifier

        for surface_out_id, surface_in_dict in self.surface_out_dictionary.items():
            surface_out_dict = surface_in_dict[structure_in_id]
            
            # check if surface out dictionary is empty
            if not surface_out_dict:
                continue

            triangle_mapping_dict = self.triangle_mapping_dictionary[surface_out_id][structure_in_id]

            for bounce_id, mapping_dict in triangle_mapping_dict.items():
                if not np.all(surface_out_dict[bounce_id] == 0.):

                    mapping_in = mapping_dict['deposit_into_surface_in']
                    triangles_in = mapping_in[:,0]
                    indices_in   = mapping_in[:,1]

                    bounce_tensor = self.bounce_tensors[bounce_id]
                    surface_out         = surface_out_dict[bounce_id]
                    surface_bounced_tmp = np.tensordot(bounce_tensor,surface_out,axes=1)

                    # Obtain binary mask for geometry-outgoing vs geomotry-ingoing carriers
                    surface_mask = structure_to_push_to.inside_outside_mask
                    minimal_surface_mask = surface_mask[:,triangles_in]
                
                    # Mask bounced carriers to only pick geometry-ingoing carriers
                    surface_in_masked = minimal_surface_mask*surface_bounced_tmp

                    # Ensure the geometry-outgoing carriers sum up to the same as the bounced geometry-ingoing carriers
                    scaling_factor= (np.sum(surface_out,axis=0)/np.sum(surface_in_masked,axis=0))[np.newaxis,...]
                    scaling_factor[np.logical_not(np.isfinite(scaling_factor))]=0.
                    surface_in_masked *= scaling_factor
                
                    structure_to_push_to.surface_injection[:,triangles_in,indices_in] += surface_in_masked

            # Reset surface out to an empty dictionary
            self.surface_out_dictionary[surface_out_id][structure_in_id] = {}

