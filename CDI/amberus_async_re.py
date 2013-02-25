import random, math # Only used in the now deprected _doExchange_pair()
from pj_async_re import async_re_job
from amber_async_re import pj_amber_job, BOLTZMANN_CONSTANT
class amberus_async_re_job(pj_amber_job,async_re_job):

    def _checkInput(self):
        pj_amber_job._checkInput(self)
        #make sure AMBER umbrella sampling is wanted
        if self.keywords.get('RE_TYPE') != 'AMBERUS':
            self._exit("RE_TYPE is not AMBERUS")

        # Quick function to convert delimited state parameters to a 2d list
        def ParseStateParams(paramline, state_delimiter=':'):
            if [paramline] == paramline.split(state_delimiter):
                params = [ [item] for item in paramline.split(',') ]
            else:
                params = [ item.split(',') 
                           for item in paramline.split(state_delimiter)]
            return params

        # list of force constants
        if self.keywords.get('FORCE_CONSTANTS') is None:
            self._exit("FORCE_CONSTANTS needs to be specified")
        kbias = ParseStateParams(self.keywords.get('FORCE_CONSTANTS'))
        # list of bias positions
        if self.keywords.get('BIAS_POSITIONS') is None:
            self._exit("BIAS_POSITIONS needs to be specified")
        posbias = ParseStateParams(self.keywords.get('BIAS_POSITIONS'))
        # check that parameter dimensions match
        nR0 = len(posbias)
        nK0 = len(kbias)
        n   = self.nreplicas
        if nR0 != nK0:
            msg = ('Number of FORCE_CONSTANTS not equal to number of'
                   ' BIAS_POSITIONS')
            self._exit(msg)
        if nR0 != n or nK0 != n:
            msg = ('Expected %d umbrella parameter sets, but instead found %d'
                   ' FORCE_CONSTANTS and %d BIAS_POSITIONS'%(n,nK0,nR0))
            self._exit(msg)
                   
        # Check that all umbrellas are one temperature.
        # The reduced energies calculated in this module do not account for
        # replicas running at different temperatures.
        temp0 = self._checkStateParamsAreSame('temp0','cntrl')
        if not temp0:
            self._exit('All temperatures MUST be the same when using AMBERUS.')
        self.beta = 1./(BOLTZMANN_CONSTANT*temp0)

        # Look for a restraint template (try the basename?)
        if self.keywords.get('AMBER_RESTRAINT_TEMPLATE') is None:
            restraint_template = '%s.RST'%self.basename
        else:
            restraint_template = self.keywords.get('AMBER_RESTRAINT_TEMPLATE')
        if self.verbose:
            print 'Using restraint template file: %s'%restraint_template

        # Read the restraint template and then modify the restraint objects 
        # based on the input,
        # NB: For simplicity, hardcode all restraint files to have the same 
        #     name. Each replica will simply overwrite this file at each cycle.
        for state in self.states:
            state.AddRestraints(restraint_template)
            state.mdin.SetVariableValue('DISANG','US.RST',None)
        bias_dimensions = len(kbias[:][0])
        for n in range(self.nreplicas):
            for m in range(bias_dimensions): 
                k  = float(kbias[n][m])
                r0 = float(posbias[n][m])
                self.states[n].rstr[m].SetRestraintParameters(r0=r0,k0=k)

    def _buildInpFile(self, repl):
        """
        Builds input files for an AMBER umbrella sampling replica based on the
        current state and cycle. For simplicity, all restraint files are named
        US.RST and are remade at every cycle.
        """
        sid = self.status[repl]['stateid_current']
        cyc = self.status[repl]['cycle_current']

        # 1) Write a new restraint file for the current state
        # 2) Modify the input to print to a new output (trace) file
        title =  (' umbrella sampling restraint for replica %d in state %d'
                  ' during cycle %d'%(repl,sid,cyc))
        rst_file = 'r%d/US.RST'%repl
        self.states[sid].rstr.WriteAmberRestraintFile(rst_file,title)
        trace_file = '%s_%d.TRACE'%(self.basename,cyc)
        self.states[sid].mdin.SetVariableValue('DUMPAVE',trace_file,None)
        # NB! This needs to be done last since the mdin file is written by
        # this routine and the mdin object was modified here.
        pj_amber_job._buildInpFile(self,repl)

    def _doExchange_pair(self, repl_a, repl_b):
        """Given two replicas, swap the state ids according to the Metropolis
        criteria from swapping coordinates (proportional to bias energies). 
        """
        sid_a = self.status[repl_a]['stateid_current']
        sid_b = self.status[repl_b]['stateid_current']
           
        # do the energy evaluations
        u_aa = self._reduced_energy(sid_a,repl_a)
        u_bb = self._reduced_energy(sid_b,repl_b)
        u_ab = self._reduced_energy(sid_a,repl_b)
        u_ba = self._reduced_energy(sid_b,repl_a)

        # test for and perform the exchange
        Exchange = pj_amber_job.ReplicaExchange(self,u_aa,u_bb,u_ab,u_ba)
        if self.verbose: 
            self._printExchangePairReport(repl_a,repl_b,u_aa,u_bb,u_ab,u_ba)
        if Exchange: 
            self._swapStates(repl_a,repl_b)
            
    def _printStateReport(self, repl):
        """Report on a replica in relation to its current state (i.e. umbrella).
        """
        # Print a generic AMBER report
        pj_amber_job._printStateReport(self,repl)
        # Print an umbrella potential report
        sid = self.status[repl]['stateid_current']
        umbrella = self.states[sid].rstr
        crds = self._extractLastCoordinates(repl)
        print '\nUmbrella potential report:'
        umbrella.PrintRestraintReport(crds)
        umbrella.PrintRestraintEnergyReport(crds)
  
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
        U = [ [ 0. for j in range(self.nreplicas) ] 
              for i in range(self.nreplicas) ]
        for repl_i in replicas:
            crds_i = self._extractLastCoordinates(repl_i)
            for sid_j in states:
                # energy of replica i in state j
                u_ji = self.states[sid_j].rstr.Energy(crds_i)
                U[sid_j][repl_i] = self.beta*u_ji
        return U

    def _reduced_energy(self, state_i, repl_j):
        """
        Return the reduced energy of replica j in state i. 
        NB: This is NOT the same as the current state of replica i.
        """
        crds_j = self._extractLastCoordinates(repl_j)
        u_ij = self.states[state_i].rstr.Energy(crds_j)
        return self.beta*u_ij

if __name__ == '__main__':
    import sys, time

    BIGJOB_VERBOSE=100

    # Parse arguments:
    usage = "%prog <ConfigFile>"
    
    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)
    
    commandFile = sys.argv[1]

    print ""
    print "===================================="
    print "AMBER Asynchronous Replica Exchange "
    print "===================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    rx = amberus_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
