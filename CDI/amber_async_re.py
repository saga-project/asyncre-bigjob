import os, sys, random
from pj_async_re import async_re_job
import numpy as np

import copy, math

# Check that AMBER is installed and is an appropriate version
AMBERHOME = os.getenv('AMBERHOME')
if AMBERHOME is None:
    raise Exception('AMBERHOME is not set.')
sys.path.append(os.path.join(AMBERHOME,'bin'))
try:
    from chemistry.amber.readparm import rst7
except:
    raise Exception('Could not load AMBER python libraries. These are only'
                    ' available in AmberTools12 and later.')

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
            
        # Check which AMBER MD engine to use, default to sander
        engine = self.keywords.get('ENGINE').upper()
        sander_flags = [ 'AMBER', 'SANDER', 'AMBER-SANDER' ]
        pmemd_flags = [ 'PMEMD', 'AMBER-PMEMD' ]
        engine_name = ''
        if engine in sander_flags:  engine_name = 'sander'
        elif engine in pmemd_flags: engine_name = 'pmemd'
        else:                       self._exit('ENGINE is not AMBER')
        if self.spmd == 'mpi' or int(self.keywords.get('SUBJOB_CORES')) > 1:
            self.spmd = 'mpi'
            engine_name += '.MPI'
        # else just assume that a serial executable is desired
        # TODO?: Cuda

        # Check that this executable exists, etc.
        self.exe = os.path.join(AMBERHOME,'bin',engine_name)
        if not os.path.exists(self.exe) or not os.access(self.exe,os.X_OK):
            self._exit('Could not find an executable called %s, expected it to'
                       ' be at %s'%(engine_name,self.exe))

        #flag for turning off exchange
        do_exchanges = self.keywords.get('DO_EXCHANGES')
        if do_exchanges is None:
            self.do_exchanges = True
        else:
            if do_exchanges.capitalize() in ['False','No']:
                self.do_exchanges = False
            else:
                self.do_exchanges = True

        #input files
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if not (self.extfiles is None):
            if self.extfiles != '':
                self.extfiles = self.extfiles.split(',')

    def _launchReplica(self,replica,cycle):
        """Launches Amber sub-job using pilot-job
        """
        input_file = "%s_%d.inp" % (self.basename, cycle)
        out_file = "%s_%d.out" % (self.basename, cycle)
        prm_file = "%s.parm7" % self.basename
        crd_file = "%s_%d.rst7" % (self.basename, cycle-1)
        rst_file = "%s_%d.rst7" % (self.basename, cycle)
        xyz_file = "%s_%d.nc" % (self.basename, cycle)
        info_file = "%s_%d.info" % (self.basename, cycle)
        log_file = "%s_%d.log" % (self.basename, cycle)
        err_file = "%s_%d.err" % (self.basename, cycle)

        arguments = ["-O",
                     "-i", input_file,
                     "-o", out_file, 
                     "-p", prm_file, 
                     "-c", crd_file, 
                     "-r", rst_file, 
                     "-x", xyz_file, 
                     "-inf", info_file]

        #pilotjob: Compute Unit (i.e. Job) description
        compute_unit_description = {
            "executable": self.exe,
            "environment": [],
            "arguments": arguments,
            "output": log_file,
            "error": err_file,   
            "working_directory": os.getcwd()+"/r"+str(replica),
            "number_of_processes": int(self.keywords.get('SUBJOB_CORES')),
            "spmd_variation": self.spmd,
            }

        if self.keywords.get('VERBOSE') == "yes":
            print ( "Launching %s in directory %s (cycle %d)" % 
                    (self.exe.split('/')[-1], os.getcwd()+"/r"+str(replica), 
                     cycle) )

#        compute_unit=self.cds.submit_compute_unit(compute_unit_description)
        compute_unit=self.pilotcompute.submit_compute_unit(compute_unit_description)
        return compute_unit
        
    def _hasCompleted(self,replica,cycle):
        """
        Returns true if an Amber replica has completed a cycle. Basically 
        checks if the restart file exists.
        """
        rstfile = "r%d/%s_%d.rst7" % (replica, self.basename,cycle)
        if os.path.exists(rstfile):
            return True
        else:
            return False

    def _extractLastCoordinates(self,repl):
        """Returns a 3N list of coordinates from the last restart (rst7) file
        """
        cycle = self.status[repl]['cycle_current']
        rst_file = 'r%d/%s_%d.rst7'%(repl,self.basename,cycle)
        return rst7(rst_file).coords

###########################################################################
#
# Work in Progress: Gibbs sampling style exchanges (see impact_async_re.py)
#
###########################################################################
# gives random choice from a set with weight probabilities
    def _weighted_choice_sub(self,weights):
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
   
    def _gibbs_re_j(self,i,replicas):
        ee = self._computeSwapMatrix(replicas)
        ps = np.zeros(len(replicas))
        for j in range(len(replicas)): ps[j] = -ee[i][j] + ee[i][i] + ee[j][j]
        ps = np.exp(ps)
        return replicas[self._weighted_choice_sub(ps)]


    def doExchanges(self):
        """
        Perform nreps rounds of exchanges among waiting replicas using 
        some form of Gibbs sampling.
        """
        # find out which replicas are waiting
        self._update_running_no()
        if self.waiting > 1:
            replicas_waiting = []
            for k in range(self.nreplicas):
                if ( self.status[k]['running_status'] == "W" 
                     and self.status[k]['cycle_current'] > 1 ):
                    replicas_waiting.append(k)

        # backtrack cycle
        for k in replicas_waiting:
            self.status[k]['cycle_current'] -= 1
            self.status[k]['running_status'] = "E"

        nreps = 1000
        npermt = {}
        permt = {}
        for reps in range(nreps):
            for i,repl_i in enumerate(replicas_waiting):
                repl_j = self._gibbs_re_j(i,replicas_waiting)
                if repl_i != repl_j:
                    #swap state id's
                    sid_i = self.status[repl_i]['stateid_current'] 
                    sid_j = self.status[repl_j]['stateid_current']
                    self.status[repl_i]['stateid_current'] = sid_j
                    self.status[repl_j]['stateid_current'] = sid_i
                    
        ###### DEBUG
        #     curr_states = [ self.status[i]['stateid_current'] 
        #                     for i in range(self.nreplicas) ]
        #     curr_perm = str(zip(range(self.nreplicas),curr_states))
        #     # e.g. '[(0,1), (1,2), (2,0)]'

        #     if npermt.has_key(curr_perm):
        #         npermt[curr_perm] += 1
        #     else:
        #         npermt[curr_perm] = 1
        #         permt[curr_perm] = copy.copy(curr_states)

        # print ('Report for %d rounds of independence sampling of %d'
        #        ' replicas'%(nreps,len(replicas_waiting)))
        # print 'Swaps among replicas',replicas_waiting
        # #for k in replicas_waiting:
        # #    state = self.status[k]['stateid_current']
        # #    print "Replica %d is now in state %d"%(k,state)
        # #    self.umbrellas[state].PrintRestraintReport()

        # ss = 0 # total number of non-unique permutations
        # for k in npermt.keys(): ss += npermt[k]
        # ps = []
        # sumps = 0
        # for k in npermt.keys():
        #     b = permt[k] # the state list of all replicas
        #     e = 0
        #     for i in replicas_waiting: e += self._reduced_energy(b[i],i)
        #     p = math.exp(-e)
        #     ps.append(p)
        #     sumps += p
        # i = 0
        # sum1 = 0.
        # sum2 = 0.
        # DKL = 0.
        # print 'empirical exact   %diff state permutation'
        # print '-----------------------------------'
        # for k in npermt.keys():
        #     emp = npermt[k] / float(ss)
        #     exa = ps[i]/sumps
        #     pdiff = (emp-exa)*100/exa
        #     print '%8.6f %8.6f %5.1f %s'%(emp,exa,pdiff,str(k))
        #     i += 1
        #     sum1 += emp
        #     sum2 += exa
        #     DKL += emp*math.log(emp/exa)
        # print '-----------------------------------'
        # print '%8.6f %8.6f (sum)'%(sum1,sum2)
        # print 'Kullback-Liebler Divergence =',DKL
        #######


        # write input files
        for k in replicas_waiting:
            # Creates new input file for the next cycle
            # Places replica back into "W" (wait) state 
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = "W"
