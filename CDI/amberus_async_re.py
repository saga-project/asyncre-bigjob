import random, math
from pj_async_re import async_re_job
from amber_async_re import pj_amber_job
# AMBER plugins
from amberio.AmberRestraint import ReadAmberRestraintFile
# uses values from sander/src/constants.F90
BOLTZMANN_CONSTANT = 1.380658*6.0221367/4184 # in kcal/mol-K

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
        self.nreplicas = len(kbias)
        bias_dimensions = len(kbias[:][0])
        if len(posbias) != self.nreplicas:
            msg = ('Number of FORCE_CONSTANTS not equal to number of'
                   ' BIAS_POSITIONS')
            self._exit(msg)

        #simulation temperature
        if self.keywords.get('TEMPERATURE') is None:
            self._exit("TEMPERATURE is a required parameter")
        temperature = float(self.keywords.get('TEMPERATURE'))
        self.beta = 1./(BOLTZMANN_CONSTANT*temperature)

        # 1) Read the AmberRestraint template
        # 2) Modify a new object for each state/replica
        # 3) Store a list of the umbrella objects
        restraint_template = '%s.RST'%self.basename
        self.umbrellas = [ ReadAmberRestraintFile(restraint_template)
                           for n in range(self.nreplicas) ]
        for n in range(self.nreplicas):
            for m in range(bias_dimensions): 
                k  = float(kbias[n][m])
                r0 = float(posbias[n][m])
                self.umbrellas[n][m].SetRestraintParameters(r0=r0,k0=k)

    def _buildInpFile(self, replica):
        """
        Builds input file for a AMBER umbrella sampling replica based on a
        template input file, BASENAME.inp, for the specified replica and cycle.
        """
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']

        # Write a new restraint file for the current state
        restraint_file = 'r%d/%s_%d.RST'%(replica,self.basename,cycle)
        self.umbrellas[stateid].WriteAmberRestraintFile(restraint_file)

        input_template = '%s.inp'%self.basename
        input_file = 'r%d/%s_%d.inp'%(replica,self.basename,cycle)
        # read template buffer
        tfile = self._openfile(input_template,'r')
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace('@n@',str(cycle))
        # write out
        ofile = self._openfile(input_file,'w')
        ofile.write(tbuffer)
        ofile.close()
      
    def _doExchange_pair(self,repl_a,repl_b):
        """Perform exchange of bias parameters.        
        """
        if self.do_exchanges:
            sid_a = self.status[repl_a]['stateid_current']
            sid_b = self.status[repl_b]['stateid_current']
            
            # extract the latest configuration and state information 
            crds_a = self._extractLastCoordinates(repl_a)
            crds_b = self._extractLastCoordinates(repl_b)
            umbrella_a = self.umbrellas[sid_a]
            umbrella_b = self.umbrellas[sid_b]

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
        # Compute all energies needed for permutation of replicas and states
        nreplicas = len(replicas)
        ee = [ [ 0. for j in range(nreplicas) ] for i in range(nreplicas) ]
        for i,repl_i in enumerate(replicas):
            crds_i = self._extractLastCoordinates(repl_i)
            # case repl_i = repl_j:
            sid_i = self.status[repl_i]['stateid_current']
            eii =  self.beta*self.umbrellas[sid_i].Energy(crds_i) 
            ee[i][i] = eii
            for j,repl_j in enumerate(replicas[0:i]): 
                # case repl_i != repl_j:
                sid_j = self.status[repl_j]['stateid_current']
                eij = self.beta*self.umbrellas[sid_j].Energy(crds_i) 
                ee[i][j] += eij
                ee[j][i] += eij
        return ee

    def _reduced_energy(self,state_i,replica_j):
        # Return the reduced energy in state_i of crds from replica_j
        crds_j = self._extractLastCoordinates(replica_j)
        umbrella_i = self.umbrellas[state_i]
        return self.beta*umbrella_i.Energy(crds_j)

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
