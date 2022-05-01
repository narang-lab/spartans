from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank==0:

    import logging

    from spartans.import_database import Database
    from spartans.structure import Structure
    from spartans.mesh_utils import set_interfaces
    from spartans.logger_builder import setup_logger
    from spartans.parser_utils import print_spartans_logo

    config_dict = {
            'directories': {
                'database_directory': './in-files',
                'output_directory': './out-files',
                'log_directory': './log-files'
                }
            }

    # Set up log files
    log_directory = config_dict['directories']['log_directory']
    setup_logger('main_log', '{}/main_log.txt'.format(log_directory))
    main_log = logging.getLogger('main_log')
    
    # Print Spartans logo
    print_spartans_logo(main_log)

    # Set up Database object
    main_log.info('Setting up Database dictionary.')
    db=Database(config_dict)

    db.dry_run()         # Print dataset shapes to import
    db.import_datasets() # Actually import database datasets
    database_dict = db.database_dict

    # Set up Structure objects
    main_log.info('Setting up Structure objects.')
    structure_keys = [struct_key for struct_key in database_dict.keys() if 'structure(' in struct_key]
    structures = [Structure(comm,config_dict,database_dict,struct_key, write_debug_outputs=False,debug_frequency=1) for struct_key in structure_keys]

    # Assign meshes to structures
    mesh_identifier_to_structure_identifier = {
                database_dict[struct_key]['structure_object'].mesh.mesh_identifier : struct_key
                for struct_key in structure_keys
                }
    database_dict['mesh_to_structure'] = mesh_identifier_to_structure_identifier
    
    # Run initialization runs for all structures
    main_log.info('Running initialization runs.')
    for struct in structures:
        struct.run_initialization()

    # Set up interface objects
    main_log.info('Setting up Interface objects.')
    set_interfaces(config_dict,database_dict)

    # Set up Scheduler
    # Note: change this to fit your needs, and/or define new schedulers in scheduler.py
    main_log.info('Setting up Scheduler.')

    from spartans.scheduler import Unweighted_simple as Scheduler
    scheduler = Scheduler(
            min_condition=50, 
            structures=structures
            )

    # Start Scheduler while loop
    main_log.info('Starting main Scheduler loop.')
    while (ret := scheduler.next_job(structures))['arguments'] != 'stop' :
        ret['function'](*ret['arguments'])

    # Finalize MPI job
    comm.Barrier()
    args={'inFile':'stop','outFile':None}
    comm.bcast(args,root=0)
    
    main_log.info('Finished calculation.')
    MPI.Finalize()

else:
    from spartans import transport

    args = {'inFile':None, 'outFile':None}

    while True:
        comm.Barrier()
        args = comm.bcast(args,root=0)
        
        if args['inFile'] == 'stop':
            break

        transport(comm,args['inFile'], args['outFile'])

    MPI.Finalize()

