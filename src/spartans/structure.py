import numpy as np
import h5py
import logging
from parse import *

from spartans.material import Material
from spartans.mesh import Mesh
from spartans.mesh_utils import construct_inside_outside_mask,map_triangles_from_user_to_cpp, map_triangles_from_cpp_to_user
from spartans.interface import Interface
from spartans.logger_builder import setup_logger
from spartans.run_utils import write_inputs
from spartans.mesh_utils import body_volume_integral, surface_area_integral
from spartans import transport

class Structure:
    """"
    A Structure is an object encompassing the mesh + material objects.
    """

    def __init__(self, communicator, config_dict, database_dict, structure_identifier,write_debug_outputs=False, debug_frequency=1):
        '''
        Each Structure has a unique identifier of the form
        structure(mesh({mesh_identifier})--material({material_identifier}))
        '''

        # Set up communicator
        self.comm = communicator

        # set up config/database dictionaries
        self.config_dict = config_dict
        self.database_dict = database_dict
        self.structure_identifier = structure_identifier
        self.write_debug_outputs = write_debug_outputs
        self.debug_frequency = debug_frequency
        
        # parse() is the opposite of format()
        mesh_identifier, material_identifier = parse('structure({}--{})',self.structure_identifier).fixed
        simple_mesh_id, simple_mat_id = parse('structure(mesh({})--material({}))',self.structure_identifier).fixed
        # strip parentheses for io
        self.structure_identifier_iosafe = simple_mesh_id + '--' + simple_mat_id

        # set up structure log file
        log_directory = self.config_dict['directories']['log_directory']
        setup_logger('{}_log'.format(self.structure_identifier_iosafe), '{}/{}_log.txt'.format(log_directory,self.structure_identifier_iosafe))
        self.structure_log = logging.getLogger('{}_log'.format(self.structure_identifier_iosafe))
        self.structure_log.info('Building structure object for {}'.format(self.structure_identifier))

        # Create Material and Mesh objects
        self.material = Material(self.database_dict, material_identifier)
        self.mesh = Mesh(self.database_dict, mesh_identifier)
        self.usr_surface_indexing = self.mesh.surface_triangle_indices
        self.inside_outside_mask = construct_inside_outside_mask(self.mesh.surface_triangles_normals,self.material.wavevectors)
        
        # Injections are optional, although we require at-least one injection in one of the structures.
        # Since we'll need an initialization run w/ only body_injection to get cpp_surface_indexing,
        # we default to an array of zeros for body_injection.
        self.body_injection=self.database_dict[self.structure_identifier].get('body_injection',
                np.zeros((self.material.num_carriers,self.mesh.vertices.shape[0])))

        # By contrast, we default to None for surface_injection to help with the logic below.
        self.surface_injection=self.database_dict[self.structure_identifier].get('surface_injection',None)

        # add structure object to database dictionary
        self.database_dict[self.structure_identifier]['structure_object'] = self
        
    def __str__(self):
        '''
        Prints structure identifier
        '''
        return self.structure_identifier


    def run_initialization(self):
        '''
        Submits a single initilization run
        and reads surface indexing and total number of carriers to normalize injection with.
        '''

        # read directories and parameters needed for run
        output_directory = self.config_dict['directories']['output_directory']
        self.structure_log.info('Preparing initialization run for {}.'.format(self.structure_identifier))

        # write cpp file with zero Si
        surface_in_cpp = np.zeros((self.material.num_carriers,*self.usr_surface_indexing.shape))
        write_inputs(output_directory, self.structure_identifier_iosafe, 
                     self.mesh.vertices, self.mesh.tetrahedra_indices, self.body_injection, 
                     surface_in_cpp, self.material.wavevectors, self.material.num_carriers)

        # perform transport
        self.structure_log.info('Running initialization run for {}.'.format(self.structure_identifier))
        inFile ="{outDir}/run_{structureId}.h5".format(outDir=output_directory,structureId=self.structure_identifier_iosafe)
        outFile="{outDir}/out_{structureId}.h5".format(outDir=output_directory,structureId=self.structure_identifier_iosafe)
        
        args   = {"inFile":inFile,"outFile":outFile}
        self.comm.Barrier()
        self.comm.bcast(args,root=0)
        
        transport(self.comm,inFile, outFile)

        # read initialization outputs
        with h5py.File('{outDir}/out_{structureId}.h5'.format(outDir=output_directory,structureId=self.structure_identifier_iosafe),'r') as f:
            self.cpp_surface_indexing = f['SoIndex'][:]
            self.carrier_count = np.sum(np.abs(f['intBSi'][:]))/2
        
        # run second initialization run if Si is non-zero
        if self.surface_injection is None:
            self.surface_injection =  np.zeros((self.material.num_carriers,*self.usr_surface_indexing.shape))
        else:
            self.run_initialization_with_surface_in()

        # normalize injections
        np.divide(self.body_injection,self.carrier_count,self.body_injection)
        np.divide(self.surface_injection,self.carrier_count,self.surface_injection)

        # initiallize empty interfaces dictionary
        self.interfaces = {}
        self.run_count = 0
        self.surface_out = None
        
        # initialize accumulated arrays
        self.body_in_accumulated = np.zeros_like(self.body_injection)
        self.body_out_accumulated = np.zeros_like(self.body_injection)
        self.surface_in_accumulated = np.zeros_like(self.surface_injection)
        self.surface_out_accumulated = np.zeros_like(self.surface_injection)
        
        with h5py.File('{outDir}/accumulated_{structureId}.h5'.format(
                outDir=output_directory,structureId=self.structure_identifier_iosafe),'w') as h5f:
            
            h5f.create_dataset('Bo_acc',data=self.body_out_accumulated)
            h5f.create_dataset('So_acc',data=self.surface_out_accumulated)
            h5f.create_dataset('Bi_acc',data=self.body_in_accumulated)
            h5f.create_dataset('Si_acc',data=self.surface_in_accumulated)
            h5f.create_dataset('Bi_initial',data=self.body_injection)
            h5f.create_dataset('Si_initial',data=self.surface_injection)
            h5f.create_dataset('verts',data=self.mesh.vertices)
            h5f.create_dataset('tets',data=self.mesh.tetrahedra_indices)
            h5f.create_dataset('SIndex',data=self.usr_surface_indexing)
            h5f.attrs['num_carriers'] = 0.

    
    def run_initialization_with_surface_in(self):
        '''
        Runs a second initialization run now that we know cpp_surface_indexing
        '''

        # read directories and parameters needed for run
        output_directory = self.config_dict['directories']['output_directory']
        self.structure_log.info('Preparing second initialization run for {}.'.format(self.structure_identifier))

        # reverse and write cpp file
        surface_in_cpp = map_triangles_from_user_to_cpp(self.cpp_surface_indexing, self.usr_surface_indexing, self.surface_injection)
        write_inputs(output_directory, self.structure_identifier_iosafe, 
                     self.mesh.vertices, self.mesh.tetrahedra_indices, self.body_injection, 
                     surface_in_cpp, self.material.wavevectors, self.material.num_carriers)

        # perform transport
        self.structure_log.info('Running second initialization run for {}.'.format(self.structure_identifier))
        inFile ="{outDir}/run_{structureId}.h5".format(outDir=output_directory,structureId=self.structure_identifier_iosafe)
        outFile="{outDir}/out_{structureId}.h5".format(outDir=output_directory,structureId=self.structure_identifier_iosafe)
        
        args   = {"inFile":inFile,"outFile":outFile}
        self.comm.Barrier()
        self.comm.bcast(args,root=0)
        
        transport(self.comm,inFile, outFile)

        # update initialization output
        with h5py.File('{outDir}/out_{structureId}.h5'.format(outDir=output_directory,structureId=self.structure_identifier_iosafe),'r') as f:
            self.carrier_count = np.sum(np.abs(f['intBSi'][:]))/2

    def run(self,weight):
        '''
        Submits a run, reads surface_out and body_out,
        and compute the next body_injection.
        '''

        # read directories and parameters needed for run
        output_directory = self.config_dict['directories']['output_directory']
        
        self.structure_log.info('Preparing run {runCount} for {structureId}.'.format(runCount=self.run_count,structureId=self.structure_identifier))
        
        # accumulate
        self.surface_in_accumulated += self.surface_injection
        self.body_in_accumulated += self.body_injection

        with h5py.File('{outDir}/accumulated_{structureId}.h5'.format(
                outDir=output_directory,structureId=self.structure_identifier_iosafe),'a') as h5f:

            h5f['Bi_acc'][...]=self.body_in_accumulated
            h5f['Si_acc'][...]=self.surface_in_accumulated
        
        # reverse and write cpp file
        surface_in_cpp = map_triangles_from_user_to_cpp(self.cpp_surface_indexing, self.usr_surface_indexing, self.surface_injection)
        write_inputs(output_directory, self.structure_identifier_iosafe, 
                     self.mesh.vertices, self.mesh.tetrahedra_indices, self.body_injection, 
                     surface_in_cpp, self.material.wavevectors, self.material.num_carriers)

        # perform transport
        self.structure_log.info('Running run {runCount} for {structureId}.'.format(runCount=self.run_count,structureId=self.structure_identifier))
        inFile ="{outDir}/run_{structureId}.h5".format(outDir=output_directory,structureId=self.structure_identifier_iosafe)
        outFile="{outDir}/out_{structureId}.h5".format(outDir=output_directory,structureId=self.structure_identifier_iosafe)
        
        args   = {"inFile":inFile,"outFile":outFile}
        self.comm.Barrier()
        self.comm.bcast(args,root=0)
        
        transport(self.comm,inFile, outFile)

        # read outputs
        with h5py.File('{outDir}/out_{structureId}.h5'.format(outDir=output_directory,structureId=self.structure_identifier_iosafe),'r') as f:
            body_out = f['Bo'][:]
            surface_out_cpp = f['So'][:]
            self.surface_out = map_triangles_from_cpp_to_user(self.cpp_surface_indexing, self.usr_surface_indexing, surface_out_cpp)
            self.carrier_count = np.sum(np.abs(f['intBSi'][:]))/2
        
        intBo=body_volume_integral(body_out,self.mesh.tetrahedra_indices,self.mesh.tetrahedra_volumes)
        intSo=surface_area_integral(self.surface_out,self.mesh.triangle_surface_areas)
        self.structure_log.info("Carrier integrals for {structureId}:    IntBo: {intBo}, IntSo: {intSo}".format(structureId=self.structure_identifier,intBo=intBo,intSo=intSo))

        # debugging outputs
        if self.write_debug_outputs and not self.run_count % self.debug_frequency:
            with h5py.File('{outDir}/debug_{structureId}_{runNumber}.h5'.format(outDir=output_directory,structureId=self.structure_identifier_iosafe,runNumber=self.run_count),'a') as h5f:
                h5f.create_dataset('Bo_after', data=body_out)
                h5f.create_dataset('So_after', data=self.surface_out)
                h5f.create_dataset('Bi_before',data=self.body_injection)
                h5f.create_dataset('Si_before',data=self.surface_injection)
        
        # Accumulate
        self.surface_out_accumulated += self.surface_out
        self.body_out_accumulated += body_out

        with h5py.File('{outDir}/accumulated_{structureId}.h5'.format(
                outDir=output_directory,structureId=self.structure_identifier_iosafe),'a') as h5f:

            h5f['Bo_acc'][...] = self.body_out_accumulated
            h5f['So_acc'][...] = self.surface_out_accumulated
            h5f.attrs['num_carriers'] = self.carrier_count

        # compute next body_injection
        self.body_injection = np.dot(self.material.mixing_matrix,body_out)
        self.body_injection = np.einsum('i,il->il',-self.material.inverse_diagonal,self.body_injection)
            
        # Only compute correction if weight <1
        #if weight < 1:
        #    np.multiply(body_out,(1-weight),out=body_out)
        #    np.add(self.body_injection,body_out,out=self.body_injection)

        # Average periodic boundary vertices
        #body_injection_tmp = np.zeros_like(self.body_injection)
        #body_injection_tmp[:,self.mesh.vertices_identity_ordering] += self.body_injection[:,self.mesh.vertices_periodic_ordering]/2
        #body_injection_tmp[:,self.mesh.vertices_identity_ordering] += self.body_injection[:,self.mesh.vertices_identity_ordering]/2
        #self.body_injection[:,self.mesh.vertices_identity_ordering] = body_injection_tmp[:,self.mesh.vertices_identity_ordering]
        
        # reset surface_injection
        self.surface_injection = np.zeros_like(self.surface_injection)

        self.run_count+=1
