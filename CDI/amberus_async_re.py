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
        temp0 = self.states[0].mdin.GetVariableValue('temp0','cntrl')
        for state in self.states[1:]:
            this_temp0 = state.mdin.GetVariableValue('temp0','cntrl')
            if this_temp0 != temp0:
                self._exit('All temperatures MUST be the same when using'
                           ' AMBERUS.\n%s contains temp0=%5.1f, was expecting'
                           ' %5.1f'%(state.filenames['mdin'],this_temp0,temp0))
        self.beta = 1./(BOLTZMANN_CONSTANT*temp0)

        # Look for a restraint template (try the basename?)
        if self.keywords.get('AMBER_RESTRAINT_TEMPLATE') is None:
            restraint_template = '%s.RST'%self.basename
        else:
            restraint_template = self.keywords.get('AMBER_RESTRAINT_TEMPLATE')
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

    def _buildInpFile(self, replica):
        """
        Builds input files for an AMBER umbrella sampling replica based on the
        current state and cycle. For simplicity, all restraint files are named
        US.RST and are remade at every cycle.
        """
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']

        # 1) Write a new restraint file for the current state
        # 2) Modify the input to print to a new output (trace) file
        title =  (' umbrella sampling restraint for replica %d in state %d'
                  ' during cycle %d'%(replica,stateid,cycle))
        rst_file = 'r%d/US.RST'%replica
        self.states[stateid].rstr.WriteAmberRestraintFile(rst_file,title)
        trace_file = '%s_%d.TRACE'%(self.basename,cycle)
        self.states[stateid].mdin.SetVariableValue('DUMPAVE',trace_file,None)
        # NB! This needs to be done last since the mdin file is written by
        # this routine and the mdin object was modified here.
        pj_amber_job._buildInpFile(self,replica)

    def _doExchange_pair(self,repl_a,repl_b):
        """Given two replicas, swap the state ids according to the Metropolis
        criteria from swapping coordinates (proportional to bias energies). 
        """
        sid_a = self.status[repl_a]['stateid_current']
        sid_b = self.status[repl_b]['stateid_current']
            
        # extract the latest configuration and state information 
        crds_a = self._extractLastCoordinates(repl_a)
        crds_b = self._extractLastCoordinates(repl_b)
        umbrella_a = self.states[sid_a].rstr
        umbrella_b = self.states[sid_b].rstr

        # do the energy evaluations
        u_aa = umbrella_a.Energy(crds_a)
        u_ab = umbrella_a.Energy(crds_b)
        u_ba = umbrella_b.Energy(crds_a)
        u_bb = umbrella_b.Energy(crds_b)
        delta = (u_ab + u_ba) - (u_aa + u_bb)
        u = self.beta*delta
       
        # test for and perform the exchange
        Exchange = True
        P_ab = 1.
        csi = 0.
        if u > 0.:
            csi = random.random()
            P_ab = math.exp(-u)
            if csi > P_ab: Exchange = False
        if Exchange:
            sid_a = self.status[repl_a]['stateid_current']
            sid_b = self.status[repl_b]['stateid_current']
            self.status[repl_a]['stateid_current'] = sid_b
            self.status[repl_b]['stateid_current'] = sid_a

        if self.keywords.get('VERBOSE') == 'yes':
            # extract the actual coordinates for reporting purposes
            print ('======================================================='
                   '=========================')
            print 'Exchange Attempt : Replicas (%d,%d)'%(repl_a,repl_b)
            print 'Replica %d Coordinate/Umbrella Info:'%repl_a
            umbrella_a.PrintRestraintReport(crds_a)
            umbrella_a.PrintRestraintEnergyReport(crds_a)
            print 'Replica %d Coordinate/Umbrella Info:'%repl_b
            umbrella_b.PrintRestraintReport(crds_b)
            umbrella_b.PrintRestraintEnergyReport(crds_b)
            print ('======================================================='
                   '=========================')
            print 'U_a(x_b) - U_a(x_a) + U_b(x_a) - U_b(x_b) = %f kcal/mol'%delta
            if Exchange:
                print 'Accepted! P(a<->b) = %f >= %f'%(P_ab,csi)
                print ('New states for replicas (%s,%s) are (%s,%s)'
                       %(repl_a,repl_b,
                         self.status[repl_a]['stateid_current'], 
                         self.status[repl_b]['stateid_current']))
            else:
                print 'Rejected! P(a<->b) = %f < %f'%(P_ab,csi)
            print ('======================================================='
                   '=========================')
         
    def _computeSwapMatrix(self,replicas):
        """
        Compute the swap matrix U = (u_ij):
        u_ij = u_i(x_i) (i == j) and u_i(x_j) + u_j(x_i) (i != j)
       
        Here it is assumed that u_i(x) = beta[U_0(x) + U_i(x)], so that 
        differences of the matrix elements only involve the bias potentials U_i:
        
        u_ii + u_jj - u_ij = beta[U_0(x_i) + U_i(x_i)] 
                             + beta[U_0(x_j) + U_j(x_j)]
                             - beta[U_0(x_j) + U_i(x_j)]
                             - beta[U_0(x_i) + U_j(x_i)]
                           =  beta[U_i(x_i) + U_i(x_j) - U_i(x_j) - U_j(x_i)]
        """
        # Compute all energies needed for permutation of replicas and states
        nreplicas = len(replicas)
        ee = [ [ 0. for j in range(nreplicas) ] for i in range(nreplicas) ]
        for i,repl_i in enumerate(replicas):
            crds_i = self._extractLastCoordinates(repl_i)
            # case repl_i = repl_j:
            sid_i = self.status[repl_i]['stateid_current'] 
            eii = self.beta*self.states[sid_i].rstr.Energy(crds_i)
            ee[i][i] = eii
            for j,repl_j in enumerate(replicas[0:i]): 
                # case repl_i != repl_j:
                sid_j = self.status[repl_j]['stateid_current'] 
                eij = self.beta*self.states[sid_j].rstr.Energy(crds_i)
                ee[i][j] += eij
                ee[j][i] += eij
        return ee

    def _reduced_energy(self,state_i,replica_j):
        # Return the reduced energy in state_i of crds from replica_j
        crds_j = self._extractLastCoordinates(replica_j)
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
