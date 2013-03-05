import os, sys, random
from pj_async_re import async_re_job
import numpy as np
import math
#import copy  # only needed for debug
from amberio.ambertools import AMBERHOME,KB,rst7
from amberio.amberrun import ReadAmberGroupfile, AmberRun

__all__ = 'pj_amber_job'

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        
        # TODO: Move this to async_re_job?
        self.verbose = False
        if self.keywords.get('VERBOSE') == "yes": self.verbose = True
        
        # ===========================
        # Set up the AMBER executable
        # ===========================
        engine = self.keywords.get('ENGINE').upper()
        # TODO?: Cuda
        supported_amber_engines = {'AMBER'        : 'sander', 
                                   'SANDER'       : 'sander', 
                                   'AMBER-SANDER' : 'sander',
                                   'PMEMD'        : 'pmemd', 
                                   'AMBER-PMEMD'  : 'pmemd'}
        if supported_amber_engines.has_key(engine):
            engine = supported_amber_engines[engine]
        else:
            self._exit('Requested ENGINE is not from AMBER (sander or pmemd)')

        if self.spmd == 'mpi' or int(self.keywords.get('SUBJOB_CORES')) > 1:
            self.spmd = 'mpi'
            engine += '.MPI'
        # else just assume that a serial executable is desired

        # Check that this executable exists, etc.
        self.exe = os.path.join(AMBERHOME,'bin',engine)
        if not os.path.exists(self.exe) or not os.access(self.exe,os.X_OK):
            self._exit('Could not find an executable: %s\nExpected it to'
                       ' be at %s'%(engine,self.exe))

        # ========================================================
        # Set up the general state/replica information - 2 methods
        # ========================================================
        # (1) If present, read the AMBER groupfile and define the states,
        if self.keywords.get('AMBER_GROUPFILE') is not None:
            groupfile = self.keywords.get('AMBER_GROUPFILE')
            self.states = ReadAmberGroupfile(groupfile,engine)
            self.nreplicas = len(self.states)
            if self.verbose:
                print ('Creating %d replicas from AMBER groupfile: %s'
                       %(self.nreplicas,groupfile))
        # (2) otherwise assume that the states can be inferred from the
        # extfiles and input from a specific application (e.g. umbrella 
        # sampling).
        else:
            if self.nreplicas is None:
                self._exit('Could not determine the replica count from the'
                           ' input provided (set NREPLICAS directly or provide'
                           ' an AMBER groupfile)')
                
            # These are the bare minimum files that can define a state 
            files = {'mdin' : None, 'prmtop' : None, 'inpcrd' : None}
            # First, try to match against the known extfiles based off of 
            # their file extensions.
            exts = { 'mdin'   : ['.in','.inp','.mdin'],
                     'prmtop' : ['.prmtop','.parm7','.parm'],
                     'inpcrd' : ['.inpcrd','.rst7','.mdcrd','.crd','rst'] }
            if self.extfiles is not None:
                for file in self.extfiles:
                    basename,ext = os.path.splitext(file)
                    if ext in exts['mdin']:
                        files['mdin'] = file
                    elif ext in exts['prmtop']:
                        files['prmtop'] = file
                    elif ext in exts['inpcrd']:
                        files['inpcrd'] = file
            # Next, check in the current directory using the basename.
            for input in exts.keys():
                for ext in exts[input]:
                    file = '%s%s'%(self.basename,ext)
                    if os.path.exists(file) and files[input] is None:
                        files[input] = file
            # There's no way to delineate reference and input coordinates w/o
            # using a groupfile. Best guess is to assume they are the same.
            files['ref'] = files['inpcrd']
            # Give an error if no match of required files could be made.
            for input in exts.keys():
                if files[input] is None:
                    self._exit('Could not identify a(n) %s file based on the'
                               ' given input!'%input)
            # Define states with these file names
            if self.verbose:
                print ('Creating %d replicas using the provided'
                       ' ENGINE_INPUT_EXTFILES and ENGINE_INPUT_BASENAME'
                       %self.nreplicas)
            self.states = [ AmberRun(mode='-O',engine=engine,**files)
                            for n in range(self.nreplicas) ]

    def _buildInpFile(self, repl):
        """
        For a given replica:
        1) determine the current state 
        2) write a new mdin file (change to a restart input if cycle > 1)
        3) link to a new prmtop 
        4) link to a new ref file (as needed)
        5) link to the inpcrd from cycle = 0 if cycle = 1
        """
        # convenient shorthand for current information
        sid = self.status[repl]['stateid_current']
        cyc = self.status[repl]['cycle_current']
        new_state = self.states[sid] # AmberRun object defining the state
        new_mdin = 'r%d/mdin'%repl
        new_prmtop = new_state.filenames['prmtop']
        new_refcrd = new_state.filenames['ref']
        new_inpcrd = new_state.filenames['inpcrd']
        
        # Write or link to new state inputs, this includes:
        # 1) a new mdin file (a restart if cycle > 1)
        # 2) a new prmtop file 
        # 4) new ref coordinates (depends on mdin contents) 
        # 5) input coordinates (only if the cycle = 1)
        if cyc > 1: new_state.SetRestart()
        new_state.mdin.WriteAmberMdinFile(new_mdin)

        self._linkReplicaFile('prmtop',new_prmtop,repl)

        if new_state.HasReferenceCoordinates():
            self._linkReplicaFile('refc',new_refcrd,repl)

        if cyc == 1:
            inpcrd = '%s_0.rst7'%self.basename
            self._linkReplicaFile(inpcrd,new_inpcrd,repl)

    def _launchReplica(self,repl,cyc):
        """
        Launch an AMBER sub-job using pilot-job. 

        The input files for AMBER that define a state are assumed to be the 
        default names mdin, prmtop, and refc. These files are always re-written 
        or symlinked to in _buildInpFile().
        """
        # Working directory for this replica
        wdir = '%s/r%d'%(os.getcwd(),repl)

        # Cycle dependent input and output file names
        inpcrd = '%s_%d.rst7'%(self.basename,cyc-1)
        mdout  = '%s_%d.out'%(self.basename,cyc)
        mdcrd  = '%s_%d.nc'%(self.basename,cyc)
        rstrt  = '%s_%d.rst7'%(self.basename,cyc)
        stdout = '%s_%d.log'%(self.basename,cyc)
        stderr = '%s_%d.err'%(self.basename,cyc)

        args = ['-O','-c',inpcrd,'-o',mdout,'-x',mdcrd,'-r',rstrt]

        #pilotjob: Compute Unit (i.e. Job) description
        cpt_unit_desc = {
            "executable": self.exe,
            "environment": ['AMBERHOME=%s'%AMBERHOME],
            "arguments": args,
            "output": stdout,
            "error": stderr,   
            "working_directory": wdir,
            "number_of_processes": int(self.keywords.get('SUBJOB_CORES')),
            "spmd_variation": self.spmd,
            }

        compute_unit = self.pilotcompute.submit_compute_unit(cpt_unit_desc)
        return compute_unit
        
    def _hasCompleted(self,repl,cyc):
        """
        Returns true if an Amber replica has completed a cycle. Basically 
        checks if the restart file exists.
        """
        # TODO: Parse the output file and look for more sure signs of 
        #       completion?
        rst = 'r%d/%s_%d.rst7'%(repl,self.basename,cyc)
        if os.path.exists(rst):
            return True
        else:
            return False

    def _linkReplicaFile(self, link_filename, real_filename, repl):
        """
        Link the file at real_filename to the name at link_filename in the
        directory belonging to the given replica. If a file is already linked
        to this name (e.g. from a previous cycle), remove it first.
        """
        # TODO: move this to pj_async_re?
        real_filename = '../%s'%real_filename
        os.chdir('r%d'%repl)
        if os.path.exists(link_filename): os.remove(link_filename)
        os.symlink(real_filename,link_filename)
        os.chdir('..')

    def _extractLastCoordinates(self,repl):
        """
        Returns a 3N list of coordinates from the last restart (rst7) file of a
        given replica.
        """
        cyc = self.status[repl]['cycle_current']
        rst = 'r%d/%s_%d.rst7'%(repl,self.basename,cyc)
        return rst7(rst).coords

    def _checkStateParamsAreSame(self, variable, namelist):
        """
        Returns False if any two states have different values of a variable in
        the specified namelist. If all states have the same value, then that 
        value is returned.

        This routine can be useful if a particular exchange protocol assumes
        that certain state parameters (e.g. temperature) are the same in all 
        states.
        """
        value = self.states[0].mdin.GetVariableValue(variable,namelist)
        for state in self.states[1:]:
            this_value = state.mdin.GetVariableValue(variable,namelist)
            if this_value != value: return False
        return value

###########################################################################
#
# Work in Progress: Gibbs sampling style exchanges (see impact_async_re.py)
#
###########################################################################
    def _replicasWaiting(self):
        # find out which replicas are waiting
        self._update_running_no()
        # Return all replicas that are waiting AND have run at least once.
        replicas_waiting = [ k for k in range(self.nreplicas)
                             if ( self.status[k]['running_status'] == "W" 
                                  and self.status[k]['cycle_current'] > 1 ) ]
        return replicas_waiting

    def _statesWaiting(self):
        # find out which replicas are waiting
        self._update_running_no()
        # Return all states occupied by replicas that are waiting AND have run
        # at least once.
        states_waiting = [ self.status[k]['stateid_current']
                           for k in range(self.nreplicas)
                           if ( self.status[k]['running_status'] == "W" 
                                and self.status[k]['cycle_current'] > 1 ) ]
        return states_waiting

    def _weighted_choice_sub(self,weights):
        # gives random choice from a set with weighted probabilities
        rnd = random.random() #* sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
   
    # def _gibbs_re_j(self,repl_i,replicas,U):
    #     # Select a replica swap partner based off independence sampling
    #     #
    #     # Pswap = exp(-duij) / sum(exp(-duij))
    #     ee = U
    #     ps = np.zeros(len(replicas))
    #     for repl_j in replicas: 
    #         ps[j] = ee[i][i] + ee[j][j] - ee[i][j] - e[j][i]
    #     ps = np.exp(-ps)
    #     ps /= ps.sum()
    #     return replicas[self._weighted_choice_sub(ps)]
      
    @staticmethod
    def MetropolisAcceptance(du):
        """
        Return whether a Monte Carlo move is accepted or not based off of the
        Metropolis criteria:

        P_accept = min{1,exp(-du)}

        where du = u_new - u_old is the difference in the new and old
        REDUCED (i.e. in kT units) energy difference.
        """
        acceptMove = True
        if du > 0.:
            Pacc = math.exp(-du)
            rand = random.random()
            if rand > Pacc: 
                acceptMove = False
        return acceptMove

    @staticmethod
    def ReplicaExchange(u11,u22,u12,u21):
        du = u12 + u21 - u11 - u22
        return pj_amber_job.MetropolisAcceptance(du)

    def _computeSwapMatrix(self,replicas,states):
        # TODO: Implement this with single point energy calls to AMBER
        pass

    def _swapStates(self, repl_a, repl_b):
        sid_a = self.status[repl_a]['stateid_current']
        sid_b = self.status[repl_b]['stateid_current']
        self.status[repl_a]['stateid_current'] = sid_b
        self.status[repl_b]['stateid_current'] = sid_a
        if self.verbose:
            print ('Exchange accepted!\nNew states for replicas (%s,%s) are'
                   ' (%s,%s)'%(repl_a,repl_b,sid_a,sid_b))

    def _printStateReport(self, repl):
        """Report on a replica in relation to its current state.
        """
        sid = self.status[repl]['stateid_current']
        cyc   = self.status[repl]['cycle_current']
        print 'Replica : %5d State : %5d Last Cycle : %d'%(repl,sid,cyc)
        print 'Input files:'
        print 'mdin   : %s'%self.states[sid].filenames['mdin']
        print 'prmtop : %s'%self.states[sid].filenames['prmtop']
        print 'refc   : %s'%self.states[sid].filenames['ref']

    def _printExchangePairReport(self, repl_a, repl_b, u_aa, u_bb, u_ab, u_ba):
        delta = (u_ab + u_ba) - (u_aa + u_bb)
        # print '================================================================'
        print 'Exchange Attempt : Replicas (%d,%d)'%(repl_a,repl_b)
        # print '================================================================'
        # self._printStateReport(repl_a)
        # print
        # self._printStateReport(repl_b)
        # print '================================================================'
        print 'u_a(x_b) - u_a(x_a) + u_b(x_a) - u_b(x_b) = %f kT'%delta
    
    def doExchanges(self):
        """
        Perform nreps rounds of exchanges among waiting replicas using 
        some form of Gibbs sampling.
        """
        replicas_waiting = self._replicasWaiting()
        states_waiting = self._statesWaiting()
        nwaiting = len(replicas_waiting)
        if nwaiting != len(states_waiting):
            print ('Something weird happened, number of waiting replicas and'
                   ' states are different!')
        # If not enough replicas are available, exit (this shouldn't happen
        # a lot in practice if lots of replicas are running).
        # if self.waiting < 2: return 0
        if len(replicas_waiting) < 2: return 0

        # backtrack cycle
        for k in replicas_waiting:
            self.status[k]['cycle_current'] -= 1
            self.status[k]['running_status'] = "E"

        # Chodera and Shirts suggested that K^3-K^5 iterations of random pairs
        # should approach independence sampling. For tests of up to 4 replicas,
        # the similarity of the empirical and exact distribution are
        # similar with a Kullback-Liebler divergence on the order of 0.1 or 
        # less.
        nreps = nwaiting**3
        # npermt = {}
        # permt = {}
        U = self._computeSwapMatrix(replicas_waiting,states_waiting)
        # Remember that U is nreplicas x nreplicas with rows corresponding
        # to the replica ids and columns corresponding to the STATIC state
        # ids (state id and replica id can and will be different)
        for reps in range(nreps):
            #
            # Independence sampling-type swaps
            #
            # for repl_i in replicas_waiting:
            #     sid_i = self.status[repl_i]['stateid_current']
            #     # Compute the swap probablility (ps) for all states
            #     ps = np.zeros(nwaiting)
            #     for j,repl_j in zip(range(nwaiting),replicas_waiting):
            #         sid_j = self.status[repl_j]['stateid_current']
            #         ps[j] = ( U[sid_i][repl_j] + U[sid_j][repl_i] 
            #                   - U[sid_i][repl_i] - U[sid_j][repl_j] )
            #     ps = np.exp(-ps)
            #     ps /= ps.sum()
            #     repl_j = replicas_waiting[self._weighted_choice_sub(ps)]
            #     if repl_i != repl_j:
            #         if self.verbose: 
            #             self._printExchange_PairReport(repl_i,repl_j,
            #                                           U[sid_i][repl_i],
            #                                           U[sid_j][repl_j],
            #                                           U[sid_i][repl_j],
            #                                           U[sid_j][repl_i])
            #         self._swapStates(repl_i,repl_j)
            #
            # Traditional Metropolis-type swaps
            #
            repl_i,repl_j = random.sample(replicas_waiting,2)
            sid_i = self.status[repl_i]['stateid_current']
            sid_j = self.status[repl_j]['stateid_current']
            exchange = pj_amber_job.ReplicaExchange(U[sid_i][repl_i],
                                                    U[sid_j][repl_j],
                                                    U[sid_i][repl_j],
                                                    U[sid_j][repl_i])
            # if self.verbose: self._printExchangePairReport(repl_i,repl_j,
            #                                                U[sid_i][repl_i],
            #                                                U[sid_j][repl_j],
            #                                                U[sid_i][repl_j],
            #                                                U[sid_j][repl_i])
            if exchange: self._swapStates(repl_i,repl_j)
        ###### DEBUG
        #     # list of current states of ALL replicas
        #     curr_states = [ self.status[i]['stateid_current'] 
        #                     for i in replicas_waiting ]
        #     # list of tuples of replicas and sids
        #     # e.g. '[(0,1), (1,2), (2,0)]'
        #     curr_perm = str(zip(replicas_waiting,curr_states))
           
        #     # If the permutation has been observed, add a count 
        #     if npermt.has_key(curr_perm):
        #         npermt[curr_perm] += 1
        #     # otherwise add this permutation and store it as a possibility
        #     else:
        #         npermt[curr_perm] = 1
        #         permt[curr_perm] = copy.copy(curr_states)

        # print '================================================================'
        # print ('Report for %d rounds of independence sampling of %d'
        #        ' replicas'%(nreps,len(replicas_waiting)))
        # print ('Swaps among replicas %s in states %s N! = %d permutations'
        #        %(str(replicas_waiting),str(curr_states),
        #          math.factorial(len(replicas_waiting))))

        # from itertools import permutations
        # emp  = []
        # exact = []
        # sumps = 0
        # ss = 0
        # for state_perm in permutations(curr_states):
        #     perm = str(zip(replicas_waiting,state_perm))
        #     # emperical distribution observed here
        #     if npermt.has_key(perm):
        #         emp.append(npermt[perm])
        #         ss += npermt[perm]
        #     else:
        #         emp.append(0.)
        #     # exact boltzmann weight of all permutations
        #     e = 0
        #     for i,j in zip(replicas_waiting,state_perm):
        #         e += self._reduced_energy(j,i)
        #     p = math.exp(-e)
        #     exact.append(p)
        #     sumps += p
        # exact = [ p / sumps for p in exact ]
        # emp   = [ float(p) / ss for p in emp ]
        # sum1 = 0.
        # sum2 = 0.
        # DKL = 0.
        # print '   empirical exact   state permutation'
        # print '--------------------------------------'
        # dP = 1.e-4 # this will show up as 0 but contribute a lot to DKL
        # for k,state_perm in enumerate(permutations(curr_states)):
        #     perm = str(zip(replicas_waiting,state_perm))
        #     print '%3d %8.3f %8.3f %s'%(k+1,emp[k],exact[k],perm)
        #     sum1 += emp[k]
        #     sum2 += exact[k]
        #     if emp[k] > 0.: empk = emp[k]
        #     else:           empk = dP
        #     if exact[k] > 0.: exactk = exact[k]
        #     else:             exactk = dP
        #     DKL += empk*math.log(empk/exactk)
        # print '--------------------------------------'
        # print ('    %8.3f %8.3f (sum) Kullback-Liebler Divergence = %f'
        #        %(sum1,sum2,DKL))
        # print '================================================================'
        #######

        # write input files
        for k in replicas_waiting:
            # Creates new input file for the next cycle
            # Places replica back into "W" (wait) state 
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = "W"
            
        return 1
