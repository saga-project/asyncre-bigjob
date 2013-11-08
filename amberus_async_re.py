import os

import amberio.ambertools as at
from amber_async_re import pj_amber_job, DISANG_NAME, DUMPAVE_EXT, _exit

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


class amberus_async_re_job(pj_amber_job):

    def _checkInput(self):
        pj_amber_job._checkInput(self)
        #    Check that all umbrellas are the same temperature. The reduced 
        # energies calculated in this module do not account for replicas 
        # running at different temperatures.
        #
        temp0 = self._state_params_are_same('temp0','cntrl')
        if not temp0:
            _exit('All temperatures MUST be the same when using AMBER-US.')
        self.beta = 1./(at.KB*temp0)

        
        # Umbrella Sampling Parameters
        #
        bias_filename = self.keywords.get('BIAS_FILE')
        if bias_filename is not None:
            try:
                bias_positions,force_constants = \
                    _parse_bias_file(bias_filename)
            except IOError:
                _exit('Problem reading BIAS_FILE = %s'%bias_filename)
        elif (self.keywords.get('FORCE_CONSTANTS') is not None
              and self.keywords.get('BIAS_POSITIONS') is not None):
            force_constants = _parse_state_params(
                self.keywords.get('FORCE_CONSTANTS'))
            bias_positions = _parse_state_params(
                self.keywords.get('BIAS_POSITIONS'))
        else:
            _exit('No bias specifications! Either BIAS_FILE or BIAS_POSITIONS '
                  'and FORCE_CONSTANTS must be specified.')

        if len(bias_positions) != len(force_constants):
            _exit('Number of FORCE_CONSTANTS not equal to number of '
                  'BIAS_POSITIONS')
        if (len(bias_positions) != self.nreplicas 
            or len(force_constants) != self.nreplicas):
            _exit('Expected %d umbrella parameter sets, but instead found '
                  '%d FORCE_CONSTANTS and %d BIAS_POSITIONS'
                  %(self.nreplicas,len(force_constants),len(bias_positions)))
        
        # Look for a restraint template (try the basename?)
        if self.keywords.get('AMBER_RESTRAINT_TEMPLATE') is None:
            restraint_template = '%s.RST'%self.basename
        else:
            restraint_template = self.keywords.get('AMBER_RESTRAINT_TEMPLATE')
        if self.verbose:
            print 'Using restraint template file: %s'%restraint_template

        # Read the restraint template and then modify the restraint objects 
        # based on the input. 
        for n,state in enumerate(self.states):
            state.add_restraints(restraint_template)
            state.mdin.set_namelist_value('DISANG',DISANG_NAME,None)
            state.rstr.set_restraint_params(r0=bias_positions[n],
                                            k0=force_constants[n])
  
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
        # U will be sparse matrix, but is convenient bc the indices of the
        # rows and columns will always be the same.
        U = [[ 0. for j in range(self.nreplicas)] 
             for i in range(self.nreplicas)]
        for repl_i in replicas:
            crds_i = self._extractLastCoordinates(repl_i)
            for sid_j in states:
                # energy of replica i in state j
                u_ji = self.states[sid_j].rstr.energy(crds_i)
                U[sid_j][repl_i] = self.beta*u_ji
        return U

#    def _computeSwapMatrix(self, replicas, states):
#        U_old = self._computeSwapMatrix_OLD(replicas,states)
#        U_new = pj_amber_job._computeSwapMatrix(self,replicas,states)
#        return U_new,U_old

    # def _reduced_energy(self, state_i, repl_j):
    #     """
    #     Return the reduced energy of replica j in state i. 
    #     NB: This is NOT the same as the current state of replica i.
    #     """
    #     crds_j = self._extractLastCoordinates(repl_j)
    #     u_ij = self.states[state_i].rstr.energy(crds_j)
    #     return self.beta*u_ij

    def _hasCompleted(self, repl, cyc):
        """Returns True if an umbrella sampling replica has completed a cycle.
        
        Basically checks if the trace file exists.
        """
        trace = 'r%d/%s_%d.%s'%(repl,self.basename,cyc,DUMPAVE_EXT)
        if os.path.exists(trace):
            return pj_amber_job._hasCompleted(self,repl,cyc)
        else:
            return False


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
