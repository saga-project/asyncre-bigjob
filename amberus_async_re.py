import os
from multiprocessing import Pool, cpu_count

from configobj import ConfigObj
from numpy import zeros, asarray

import amberio.ambertools as at
from amber_async_re import pj_amber_job, amber_states_from_configobj, \
    extract_amber_coordinates, DISANG_NAME, DUMPAVE_EXT, _exit

def _parse_state_params(paramline, state_delimiter=':'):
    """
    Return a list of restraint parameters defining a set of states given a 
    delimited string of those parameters. The parameters are assumed to be
    comma delimited, but the state delimiter can be changed.

    Example:
    >>> rstr_list = '1.00,1.00:1.00,2.00' # two states in two dimensions
    >>> _parse_state_params(rstr_list)
    [[1.0,1.0],[1.0,2.0]]
    """
    # If each state is only one dimension, make it a single element list. 
    # This makes iterating all parameter lists the same.
    if [paramline] == paramline.split(state_delimiter):
        params = [[float(param)] for param in paramline.split(',')]
    else:
        params = [[float(param) for param in state.split(',')] 
                  for state in paramline.split(state_delimiter)]
    return params

def _parse_bias_file(bias_filename):
    """
    Return a list of bias positions and force constants from a "bias file" of 
    the format:

    bias_position1 force_constant1 bias_position2 force_constant2 ...

    where any number of bias specifications are permitted. 
    """
    bias_positions = []
    force_constants = []
    for i,line in enumerate(open(bias_filename,'r')):
        try:
            line = line[:line.index('#')] # ignore comments
        except ValueError:
            pass
        if len(line) > 0: # ignore empty lines
            tokens = line.strip().split()
            if len(bias_positions) == len(force_constants) == 0:
                dim = len(tokens)/2
            else:
                if len(tokens)/2 != dim:
                    _exit('Bad harmonic bias specification on line %d of %s'
                          %(i+1,bias_filename))
            bias_positions.append([float(r0) for r0 in tokens[0::2]])
            force_constants.append([float(k0) for k0 in tokens[1::2]])
    return bias_positions,force_constants

def setup_us_states_from_configobj(states, keywords, verbose=False):
    """
    Augment a set of AMBER states to include umbrella sampling state
    information (i.e. restraint information) from a ConfigObj.
    """
    # Umbrella Sampling Parameters
    #
    nreplicas = len(states)
    bias_filename = keywords.get('BIAS_FILE')
    if bias_filename is not None:
        try:
            bias_positions,force_constants = _parse_bias_file(bias_filename)
        except IOError:
            _exit('Problem reading BIAS_FILE = %s'%bias_filename)
    elif (keywords.get('FORCE_CONSTANTS') is not None
          and keywords.get('BIAS_POSITIONS') is not None):
        force_constants = _parse_state_params(
            keywords.get('FORCE_CONSTANTS'))
        bias_positions = _parse_state_params(
            keywords.get('BIAS_POSITIONS'))
    else:
        _exit('No bias specifications! Either BIAS_FILE or BIAS_POSITIONS and '
              'FORCE_CONSTANTS must be specified.')

    if len(bias_positions) != len(force_constants):
        _exit('Number of FORCE_CONSTANTS not equal to number of BIAS_POSITIONS')
    if (len(bias_positions) != nreplicas
        or len(force_constants) != nreplicas):
        _exit('Expected %d umbrella parameter sets, but instead found %d '
              'FORCE_CONSTANTS and %d BIAS_POSITIONS'
              %(nreplicas,len(force_constants),len(bias_positions)))    

    # Look for a restraint template (try the basename?)
    basename = keywords.get('ENGINE_INPUT_BASENAME')
    if keywords.get('AMBER_RESTRAINT_TEMPLATE') is not None:
        restraint_template = keywords.get('AMBER_RESTRAINT_TEMPLATE')
    else:
        restraint_template = '%s.RST'%basename
    if verbose:
        print 'Using restraint template file: %s'%restraint_template

    # Read the restraint template and then modify the restraint objects. 
    for n,state in enumerate(states):
        state.add_restraints(restraint_template)
        state.mdin.nmr_vars['DISANG'] = DISANG_NAME
        state.rstr.set_restraint_params(r0=bias_positions[n],
                                        k0=force_constants[n])

class amberus_async_re_job(pj_amber_job):

    def _checkInput(self):
        pj_amber_job._checkInput(self)
        #    Check that all umbrellas are the same temperature. The reduced 
        # energies calculated in this module do not account for replicas 
        # running at different temperatures.
        #
        temp0 = self.states.state_params_are_same('cntrl','temp0')
        if not temp0:
            _exit('All temperatures MUST be the same when using AMBER-US.')
        self.beta = 1./(at.KB*temp0)

        # Umbrella sampling state information.
        #
        setup_us_states_from_configobj(self.states,self.keywords,self.verbose)
 
    def _computeSwapMatrix(self, replicas, states):
        """
        Compute the swap matrix U = (u_ij), where u_ij = u_i(x_j)
        
        Here it is assumed that u_i(x) = beta[U_0(x) + U_i(x)], so that 
        differences of the matrix elements only involve the bias potentials U_i:
        
        u_ii + u_jj - u_ij - u_ji 
                     = beta[U_0(x_i) + U_i(x_i)] + beta[U_0(x_j) + U_j(x_j)]
                       - beta[U_0(x_j) + U_i(x_j)] - beta[U_0(x_i) + U_j(x_i)]
                     =  beta[U_i(x_i) + U_i(x_j) - U_i(x_j) - U_j(x_i)]
        """ 
        cycles = [self.status[repl]['cycle_current'] for repl in replicas]
              
        nprocs = cpu_count()
        if nprocs >= 2*len(replicas):
            nprocs = 1
        pool = Pool(processes=nprocs)
        # Divide replicas evenly amongst processes. Add extra replicas to the
        # first few processes as needed to reach len(replicas). 
        avg_replicas_per_proc = int(len(replicas)/nprocs)
        replicas_per_proc = [avg_replicas_per_proc for n in range(nprocs)]
        for n in range(len(replicas)%nprocs):
            replicas_per_proc[n] += 1

        repl_cyc_pairs = []
        for n in range(nprocs):
            first = sum(replicas_per_proc[0:n])
            last = first + replicas_per_proc[n]
            repl_cyc_pairs.append(zip(replicas[first:last],
                                      cycles[first:last]))

        print 'Computing swap matrix on %d processor(s)...'%nprocs
        results = [pool.apply_async(_compute_columns,
                                    args=(repl_cyc_pairs[n],states,
                                          self.command_file))
                   for n in range(nprocs)]
        U = zeros([self.nreplicas,self.nreplicas])
        for result in results:
            U += asarray(result.get())
        pool.close()
        pool.join()
        return U.tolist()

    def _hasCompleted(self, repl, cyc):
        """Returns True if an umbrella sampling replica has completed a cycle.
        
        Basically checks if the trace file exists.
        """
        trace = 'r%d/%s_%d.%s'%(repl,self.basename,cyc,DUMPAVE_EXT)
        if os.path.exists(trace):
            return pj_amber_job._hasCompleted(self,repl,cyc)
        else:
            return False

def _compute_columns(replicas_and_cycles, states, command_file):
    keywords = ConfigObj(command_file)
    state_objs = amber_states_from_configobj(keywords)  
    setup_us_states_from_configobj(state_objs,keywords)
    nreplicas = len(state_objs) 
    temp0 = state_objs.state_params_are_same('cntrl','temp0')
    beta = 1./(at.KB*temp0)
    basename = keywords.get('ENGINE_INPUT_BASENAME')

    U = zeros([nreplicas,nreplicas])
    for repl_i,cyc_n in replicas_and_cycles:
        crds_i = extract_amber_coordinates(repl_i,cyc_n,basename)
        for sid_j in states:
            # energy of replica i in state j
            u_ji = state_objs[sid_j].rstr.energy(crds_i)
            U[sid_j,repl_i] = beta*u_ji
    return U

if __name__ == '__main__':
    import sys
    import time

    BIGJOB_VERBOSE=100

    start_time = time.time()
    try:
        command_file = sys.argv[1]
    except IndexError:
        print 'usage: amberus_async_re.py command_file'
        sys.exit(1)
    if len(sys.argv) > 2:
        print 'Ignoring extra option(s):' + ' '.join(sys.argv[2:])

    print 
    print '======================================'
    print 'AMBER-US Asynchronous Replica Exchange'
    print '======================================'
    print 
    print 'Started at: ' + str(time.asctime())
    print 'Input file:', command_file
    print 
    sys.stdout.flush()
    rx = amberus_async_re_job(command_file, options=None)
    rx.setupJob()
    rx.scheduleJobs()
    total_run_time = time.time() - start_time
    print 'Total Run Time: %f'%float(total_run_time)
