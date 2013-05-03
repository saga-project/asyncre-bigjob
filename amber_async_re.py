import os
import sys
import random
import math
import copy

import numpy as np

from pj_async_re import async_re_job
from mcmc import *
from amberio.ambertools import AMBERHOME,KB
from amberio.amberrun import *
from amberio.inpcrd import *

__all__ = ['pj_amber_job','SUPPORTED_AMBER_ENGINES','DISANG_NAME','DUMPAVE_EXT']

# TODO?: Cuda
SUPPORTED_AMBER_ENGINES = {'AMBER': 'sander', 'SANDER': 'sander', 
                            'AMBER-SANDER': 'sander', 'PMEMD': 'pmemd', 
                            'AMBER-PMEMD': 'pmemd'}
DISANG_NAME = 'restraint.RST' # hardcoded AMBER restraint file name
DUMPAVE_EXT = 'TRACE' # hardcoded file extension for restraint coordinates

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        
        # ===========================
        # Set up the AMBER executable
        # ===========================
        engine = self.keywords.get('ENGINE').upper()
        if SUPPORTED_AMBER_ENGINES.has_key(engine):
            engine = SUPPORTED_AMBER_ENGINES[engine]
        else:
            self._exit('Requested ENGINE (%s) is either invalid or not '
                       'currently supported.'%engine)
        if self.spmd == 'mpi' or int(self.keywords.get('SUBJOB_CORES')) > 1:
            self.spmd = 'mpi' # Always use MPI if SUBJOB_CORES > 1.
            engine += '.MPI'
        # Check that this file exists and is exectuable.
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
            self.states = read_amber_groupfile(groupfile,engine)
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
            print 'basename',self.basename
            print 'extfiles',self.extfiles
            print 'nreplicas',self.nreplicas
            self.states = amberrun_from_files(self.basename,self.extfiles,
                                              self.nreplicas,'-O',engine)
            if self.verbose:
                print ('Creating %d replicas using the provided'
                       ' ENGINE_INPUT_EXTFILES and ENGINE_INPUT_BASENAME'
                       %self.nreplicas)

    def _buildInpFile(self, repl):
        """
        For a given replica:
        1) determine the current state 
        2) write a new mdin file (change to a restart input if cycle > 1)
        3) link to a new prmtop 
        4) link to a new ref file (as needed)
        5) link to the inpcrd from cycle = 0 if cycle = 1
        """
        sid = self.status[repl]['stateid_current']
        cyc = self.status[repl]['cycle_current']
        # Make a copy of one of the existing AmberRun state templates.
        title = ' replica %d : state %d : cycle %d'%(repl,sid,cyc)
        self.states[sid].mdin.title = title
        # Modify the template as appropriate.
        if cyc > 1: 
            self.states[sid].restart()
        if self.states[sid].has_restraints():
            rstr_title = title
            rstr_file = 'r%d/%s'%(repl,DISANG_NAME)
            self.states[sid].rstr.write_amber_restraint_file(rstr_file,title)
            trace_file = '%s_%d.%s'%(self.basename,cyc,DUMPAVE_EXT)
            self.states[sid].mdin.set_namelist_value('DUMPAVE',trace_file,None)
        self.states[sid].mdin.write_amber_mdin('r%d/mdin'%repl)
        # Links
        prmtop = self.states[sid].filenames['prmtop']
        self._linkReplicaFile('prmtop',prmtop,repl)
        if self.states[sid].has_refc():
            refc = self.states[sid].filenames['ref']
            self._linkReplicaFile('refc',refc,repl)
        if cyc == 1:
            inpcrd = self.states[sid].filenames['inpcrd']
            self._linkReplicaFile('%s_0.rst7'%self.basename,inpcrd,repl) 

    def _launchReplica(self, repl, cyc):
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
        restrt = '%s_%d.rst7'%(self.basename,cyc)
        stdout = '%s_%d.log'%(self.basename,cyc)
        stderr = '%s_%d.err'%(self.basename,cyc)

        args = ['-O','-c',inpcrd,'-o',mdout,'-x',mdcrd,'-r',restrt]

        # Compute Unit (i.e. Job) description
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
        return os.path.exists(rst)

    def _extractLastCoordinates(self, repl):
        """
        Return a 3N list of coordinates from the last restart (rst7) file of a
        given replica.
        """
        cyc = self.status[repl]['cycle_current']
        rst = 'r%d/%s_%d.rst7'%(repl,self.basename,cyc)
        return rst7(rst,self.states[repl].restrt_is_binary()).coords

    def _checkStateParamsAreSame(self, variable, namelist):
        """
        Returns False if any two states have different values of a variable in
        the specified namelist. If all states have the same value, then that 
        value is returned.

        This routine can be useful if a particular exchange protocol assumes
        that certain state parameters (e.g. temperature) are the same in all 
        states.
        """
        value = self.states[0].mdin.namelist_value(variable,namelist)
        for state in self.states[1:]:
            this_value = state.mdin.namelist_value(variable,namelist)
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

    # def _weighted_choice_sub(self,weights):
    #     # gives random choice from a set with weighted probabilities
    #     rnd = random.random() #* sum(weights)
    #     for i, w in enumerate(weights):
    #         rnd -= w
    #         if rnd < 0:
    #             return i
   
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

    def _computeSwapMatrix(self, replicas, states):
        # U will be sparse matrix, but is convenient bc the indices of the
        # rows and columns will always be the same.
        U = [[ 0. for j in range(self.nreplicas)] 
             for i in range(self.nreplicas)]
        for repl_i in replicas:
            cyc = self.status[repl_i]['cycle_current']  
            inpcrd_i = '%s/r%d/%s_%d.rst7'%(os.getcwd(),repl_i,self.basename,
                                            cyc)
            for sid_j in states:
                state_dir = '%s/r%d'%(os.getcwd(),sid_j)
                # energy of replica i in state j
                os.chdir('r%d'%sid_j)
                u_ji = self.states[sid_j].reduced_energy(inpcrd_i)
                U[sid_j][repl_i] = u_ji
                os.chdir('..')
        return U

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
        # U = self._computeSwapMatrix(replicas_waiting,states_waiting)

        Unew,Uold = self._computeSwapMatrix(replicas_waiting,states_waiting)
        for i in range(self.nreplicas):
            for j in range(self.nreplicas):
                exp_OLD = Uold[i][i] + Uold[j][j] - Uold[i][j] - Uold[j][i]
                exp_NEW = Unew[i][i] + Unew[j][j] - Unew[i][j] - Unew[j][i]
                diff = exp_NEW - exp_OLD
                print 'old % 20.8f new % 20.8f diff % 20.8f'%(exp_OLD,exp_NEW,diff) 
        U = Uold

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
            exchange = replica_exchange_acceptance(U[sid_i][repl_i],
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
