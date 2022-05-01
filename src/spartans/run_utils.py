import h5py

def write_inputs(output_directory, structureId, vertices, tetrahedra_indices, body_injection, surface_in_cpp, wavevectors, num_carriers):
    with h5py.File('{outDir}/run_{structureId}.h5'.format(
        outDir=output_directory,structureId=structureId),'w') as h5f:

        h5f.create_dataset('verts',data=vertices)
        h5f.create_dataset('tets',data=tetrahedra_indices)
        h5f.create_dataset('Bi',data=body_injection)
        h5f.create_dataset('Si',data=surface_in_cpp)
        h5f.create_dataset('k',data=wavevectors)
        h5f.attrs['nEvents'] = num_carriers
