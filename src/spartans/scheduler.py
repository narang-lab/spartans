from spartans.mesh_utils import body_volume_integral, surface_area_integral
from collections import defaultdict

class Unweighted_simple:
    '''
    Default Scheduler.
    Pre-computes an evaluation order for the scheduler to cycle through:
    For each stucture, 
        we first push carriers from interfaces
        then run structure
        and finally pull carriers to interfaces
    rinse and repeat min_condition times
    '''
    def __init__(self, min_condition, structures):
        # Set values
        self.scheduler_name   = 'unweighted_simple'
        self.min_condition    = min_condition

        # Precompute fixed order of evaluations
        self.evaluation_order = []
        self.counter = 0

        # For each structure
        for struct in structures:
            # First push from structures interfaces
            for interface_id, interface_object in struct.interfaces.items():
                evaluation_dict = {
                        'function' : interface_object.push_to,
                        'arguments': (struct,1.0)
                        }
                self.evaluation_order.append(evaluation_dict)

            # Then, run structure
            evaluation_dict = {
                        'function' : struct.run,
                        'arguments': (1.0,)
                        }
            self.evaluation_order.append(evaluation_dict)

            # Finally, pull to structures interfaces
            for interface_id, interface_object in struct.interfaces.items():
                evaluation_dict = {
                        'function' : interface_object.pull_from,
                        'arguments': (struct,)
                        }
                self.evaluation_order.append(evaluation_dict)

        # Repeat for min_condition times
        self.evaluation_order *= self.min_condition
        
        # Append stopping critetion
        evaluation_dict= {
                'function' : lambda *args: None,
                'arguments': 'stop'
                }
        self.evaluation_order.append(evaluation_dict)

    def __str__(self):
        return self.scheduler_name

    def next_job(self, structures):
        '''
        Simply returns pre-computed evaluation dictionary
        '''
        evaluation_dict = self.evaluation_order[self.counter]
        self.counter += 1
        return evaluation_dict


class Weighted_simple:
    '''
    Simple scheduler accepting an iteration weight.
    Pre-computes an evaluation order for the scheduler to cycle through:
        First, we push carriers from all interfaces to all structures
        then,  we run all structures
        and finally, we pull carriers from all structures to all interfaces
    rinse and repeat min_condition times.
    The reason for this is ambiguity in surface swapping with non-unity weights
    '''
    def __init__(self, iteration_weight, min_condition, structures):

        # Set values
        self.scheduler_name   = 'weighted_simple'
        self.iteration_weight = iteration_weight
        self.min_condition    = min_condition

        # Precompute fixed order of evaluations
        self.evaluation_order = []
        self.counter = 0

        # First, push to all interfaces
        for struct in structures:
            for interface_id, interface_object in struct.interfaces.items():
                evaluation_dict = {
                        'function' : interface_object.push_to,
                        'arguments': (struct,self.iteration_weight)
                        }
                self.evaluation_order.append(evaluation_dict)

        # Then, run each structure
        for struct in structures:
            evaluation_dict = {
                        'function' : struct.run,
                        'arguments': (self.iteration_weight,)
                        }
            self.evaluation_order.append(evaluation_dict)

        # Finally, pull from all interfaces
        for struct in structures:
            for interface_id, interface_object in struct.interfaces.items():
                evaluation_dict = {
                        'function' : interface_object.pull_from,
                        'arguments': (struct,)
                        }
                self.evaluation_order.append(evaluation_dict)

        # Repeat for min_condition times
        self.evaluation_order *= self.min_condition
        
        # Append stopping critetion
        evaluation_dict= {
                'function' : lambda *args: None,
                'arguments': 'stop'
                }
        self.evaluation_order.append(evaluation_dict)

    def __str__(self):
        return self.scheduler_name

    def next_job(self, structures):
        '''
        Simply returns pre-computed evaluation dictionary
        '''
        evaluation_dict = self.evaluation_order[self.counter]
        self.counter += 1
        return evaluation_dict

class Unweighted_largest_injection:
    '''
    Simple example of how one would extend Scheduler to fit their needs.
    We precompute evaluation orders for each structure,
    and then cycle through the order for the stucture with the largest injection at the next run

    '''
    def __init__(self, max_condition, structures):
        # Set values
        self.scheduler_name   = 'Unweighted_largest_injection'
        self.max_condition    = max_condition

        # Make scheduler order template dictionary
        self.evaluation_order = defaultdict(list)
        # For each structure
        for struct in structures:
            structure_identifier = struct.structure_identifier
            # First push from structures interfaces
            for interface_id, interface_object in struct.interfaces.items():
                evaluation_dict = {
                        'function' : interface_object.push_to,
                        'arguments': (struct,1.0)
                        }
                self.evaluation_order[structure_identifier].append(evaluation_dict)

            # Then, run structure
            evaluation_dict = {
                        'function' : struct.run,
                        'arguments': (1.0,)
                        }
            self.evaluation_order[structure_identifier].append(evaluation_dict)

            # Finally, pull to structures interfaces
            for interface_id, interface_object in struct.interfaces.items():
                evaluation_dict = {
                        'function' : interface_object.pull_from,
                        'arguments': (struct,)
                        }
                self.evaluation_order[structure_identifier].append(evaluation_dict)
        
        # Scheduler counters
        # Initialize dictionary with 0 counters for all structures
        self.evaluation_order_counter  = dict.fromkeys(self.evaluation_order.keys(),0)
        # Dictionary to hold number of evaluations for each structure
        self.evaluation_order_lengths  = {key: len(value) for key, value in self.evaluation_order.items()}

        # Initialize int_carrier_counts to [1,0,0,...]
        self.int_carrier_counts        = [0] * len(structures)
        self.int_carrier_counts[0]     = 1

        # Initialize max run across all structures to 0
        self.current_max_run           = 0
    
    def __str__(self):
        return self.scheduler_name

    def _stopping_criterion_reached(self, structures):
        '''
        Example of how one would evaluate a non-predefined stopping criterion.
        Here, we check that no structure has more than max_condition runs.
        '''

        # Evaluate how many times each structure has ran
        self.current_num_runs  = []
        for struct in structures:
            self.current_num_runs.append(struct.run_count)

        # Assign min/max run and return stopping criterion
        self.current_max_run = max(self.current_num_runs)

        return self.current_max_run >= self.max_condition

    def _recompute_max_injections(self, structures):
        '''
        Helper function to compute the integrals of body_in and surface_in
        and update int_carrier_counts
        '''
        for istruct,struct in enumerate(structures):
            structure_id   = struct.structure_identifier
            int_body_in    = body_volume_integral(
                    struct.body_injection,
                    struct.mesh.tetrahedra_indices,
                    struct.mesh.tetrahedra_volumes
                    )
            int_surface_in = surface_area_integral(
                    struct.surface_injection,
                    struct.mesh.triangle_surface_areas
                    )
            int_in         = int_body_in + int_surface_in
            self.int_carrier_counts[istruct] = int_in

    def _set_max_injection_id(self, structures):
        '''
        Helper function to set identifier of structure with largest injection
        '''
        max_injection             = max(self.int_carrier_counts)
        structure_max_index       = self.int_carrier_counts.index(max_injection)
        current_structure         = structures[structure_max_index]
        self.current_structure_id = current_structure.structure_identifier

    def next_job(self, structures):
        '''
        Evaluates stopping criterion at each run.
        If not satisfied, loops through
        '''
        # return a lambda function doing nothing, with 'arguments': 'stop'
        if self._stopping_criterion_reached(structures):
            return {
                    'function' : lambda *args: None,
                    'arguments': 'stop'
                    }
        # main loop
        else:
            
            # First, check which structure has the most injection
            self._set_max_injection_id(structures)
            current_counter      = self.evaluation_order_counter[self.current_structure_id]

            # Cycle through that structures evaluations
            if current_counter < self.evaluation_order_lengths[self.current_structure_id]:
                evaluation_dict = self.evaluation_order[self.current_structure_id][current_counter]
                self.evaluation_order_counter[self.current_structure_id] +=1
                return evaluation_dict

            # If evaluations have been exhausted,
            # recompute max_injections and return first evaluation
            else:
                self._recompute_max_injections(structures)
                self._set_max_injection_id(structures)
                self.evaluation_order_counter[self.current_structure_id] = 0
                current_counter      = self.evaluation_order_counter[self.current_structure_id]
                evaluation_dict = self.evaluation_order[self.current_structure_id][current_counter]
                return evaluation_dict
