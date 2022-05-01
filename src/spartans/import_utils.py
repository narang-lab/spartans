import numpy as np
import h5py

def h5_to_dict(h5_file, path="/", dry_run=True, collect_attributes=True, log_file=None):
    '''
    Recursively read h5 object into a nested dictionary.
    Code adapted from github.com/silx-kit/silx/blob/master/silx/io/dictdump.py
    '''
    with h5py.File(h5_file,'r') as h5f:
        ddict = {}
        for key in h5f[path]:
            dict_key = key
            node = h5f[path+"/"+key]

            # attributes
            attrs = node.attrs
            if bool(attrs):
                # avoid structure object ambiguity
                if 'mesh' in attrs and 'material' in attrs:
                    dict_key = 'structure(mesh({})--material({}))'.format(attrs['mesh'],attrs['material'])
                
                # avoid connectivity object ambiguity
                elif 'outgoing_structure' in attrs and 'incoming_structure' in attrs:
                    outgoing_mesh, outgoing_material = attrs['outgoing_structure'].split('--')
                    incoming_mesh, incoming_material = attrs['incoming_structure'].split('--')
                    dict_key = 'connectivity(structure(mesh({})--material({}))--structure(mesh({})--material({})))'.format(outgoing_mesh,outgoing_material,incoming_mesh,incoming_material)
                
                # manually add attributes as key - value pairs
                else:
                    for attr_key, attr in attrs.items():
                        if isinstance(attr,np.ndarray):
                            # decode attributes stored as numpy fixed length ASCII strings
                            attr=np.array(attr,dtype='str').tolist()

                        elif isinstance(attr,np.float64):
                            # yaml friendly
                            attr=float(attr)

                        ddict['{}.{}'.format(dict_key,attr_key)]=attr

            # datasets
            if isinstance(node,h5py.Dataset):
                if dry_run:
                    # yaml-friendly shape
                    data = '('+', '.join(map(str,node.shape))+')'

                else:
                    log_file.info('Importing dataset {}'.format(path+"/"+key))
                    data = node[:]

                ddict[dict_key] = data

            # groups
            else:
                # append key to path and call again
                ddict[dict_key] = h5_to_dict(h5_file,
                                      path + "/" + key,
                                      dry_run,
                                      collect_attributes,
                                      log_file)
    return ddict

def merge_nested_dicts(a,b,path=None):
    '''
    Merges nested dictionary b into nested dictionary a, mutating a.
    Code from https://stackoverflow.com/questions/7204805/how-to-merge-dictionaries-of-dictionaries
    '''
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge_nested_dicts(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a
