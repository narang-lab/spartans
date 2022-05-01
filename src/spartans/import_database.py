import numpy as np
import functools
import yaml
import logging
from glob import glob

from spartans.logger_builder import setup_logger
from spartans.import_utils import h5_to_dict, merge_nested_dicts

class Database:
    '''
    A Database is an object storing a dictionary with all the
    material, mesh, structure, and interface datasets.
    '''

    def __init__(self, config_dict):
        '''
        Database just needs to know where to look for property file(s).
        Property file(s) need to be *h5 or *hdf5 file(s) and 
        store dataset/attributes using the hierarchical group structure.
        '''

        # Set up directories
        log_directory = config_dict['directories']['log_directory']
        database_directory = config_dict['directories']['database_directory']

        # Set up database log file
        setup_logger('database_log', '{}/database_log.txt'.format(log_directory))
        self.database_log = logging.getLogger('database_log')
        self.database_log.info('Building Database object for files in {}'.format(database_directory))

        # Find all *.h*5 property files
        self._dataset_file_names=glob("{}/*.h*5".format(database_directory))

    def dry_run(self):
        '''
        Makes a dictionary with files to be read, 
        printing out attributes and the dimensions for each dataset.
        '''

        # Read each hdf5 file recursively and make a dummy property dictionary for each
        self.database_log.info('Preparing Database dry run.')
        dicts=[]
        for file_name in self._dataset_file_names:
            self.database_log.info('Parsing file {}'.format(file_name))
            dicts.append(h5_to_dict(file_name,dry_run=True))
            
        # Merge dummy nested dictionaries and rename ambiguous group names
        dictionary=functools.reduce(merge_nested_dicts,dicts)
        self._rename_dictionary_keys(dictionary)

        # Dump to log file
        yaml_str = 'SpaRTaNS will import the following datasets:\n'
        yaml_str+= '-'*80+'\n'
        yaml_str+= yaml.safe_dump(dictionary,default_flow_style=False)
        yaml_str+= '-'*80+'\n'
        yaml_str+= '\nABORT now if this is not what you expected.\n'
        self.database_log.warning(yaml_str)
        self.dry_run_dict = dictionary

    def import_datasets(self):
        '''
        Actually load all datasets, making a nested dictionary recursively.
        '''

        # Read each hdf5 file recursively and make a property dictionary for each
        self.database_log.info('Importing Database datasets.')
        dicts=[]
        for file_name in self._dataset_file_names:
            self.database_log.info('Parsing file {}'.format(file_name))
            dicts.append(h5_to_dict(file_name,dry_run=False,log_file=self.database_log))

        # Merge nested dictionaries and rename ambiguous group names
        self.database_log.info('Merging Database dictionaries.')
        self.database_dict=functools.reduce(merge_nested_dicts,dicts)
        self._rename_dictionary_keys(self.database_dict)

    def _rename_dictionary_keys(self,dictionary):
        '''
        Wrap material({}) and mesh({}) identifiers in explicit tags.
        
        This extends available tags to:
    
        - material({material_id})
        - mesh({mesh_id})
        - structure(mesh({mesh_id})--material({material_id}))
        - connectivity(structure(mesh({mesh_id_A})--material({material_id_A}))--structure(mesh({mesh_id_B})--material({material_id_B})))
        
        After the initialization runs, we finalize available tags by creating

        - interface(structure(mesh({mesh_id_A})--material({material_id_A}))--structure(mesh({mesh_id_B})--material({material_id_B})))
        '''
        material_keys = [k for k,v in dictionary.items() if 'velocities' in v]
        for old_key in material_keys:
            dictionary['material({})'.format(old_key)] = dictionary.pop(old_key)

        mesh_keys = [k for k,v in dictionary.items() if 'vertices' in v]
        for old_key in mesh_keys:
            dictionary['mesh({})'.format(old_key)] = dictionary.pop(old_key)
