import numpy as np

class Material:
    '''
    A Material is an object representing properties of a given material
    identifier. It queries the Database and assigns properties as methods.
    '''

    def __init__(self, database_dict, material_identifier):
        '''
        Each material has a unique identifier of the form material({material_identifier})
        '''
        self.material_identifier = material_identifier

        # Required material properties
        self.velocities=database_dict[self.material_identifier]['velocities']
        self.diagonal=database_dict[self.material_identifier]['diagonal']
        self.mixing_matrix=database_dict[self.material_identifier]['mixing_matrix']

        # Optional material properties
        self.weight=database_dict[self.material_identifier].get('mixing_matrix.weight',1.0)
        self.frequencies=database_dict[self.material_identifier].get('frequencies',None)

        # Pre-compute common properties
        self.num_carriers = self.frequencies.shape[0]
        self._set_inverse_diagonal()
        self._set_wavevectors()

    def __str__(self):
        '''
        Prints material identifier
        '''
        return self.material_identifier

    def _set_inverse_diagonal(self):
        '''
        Pre-computes lifetimes
        '''
        diagonal=self.diagonal.copy()
        assert not np.any(diagonal == 0.)
        self.inverse_diagonal = 1/diagonal

    def _set_wavevectors(self):
        '''
        Pre-computes wavevectors
        '''
        velocity_mag = np.linalg.norm(self.velocities,axis=1)
        self.wavevectors = self.velocities/((self.inverse_diagonal*velocity_mag**2)[:,np.newaxis])
