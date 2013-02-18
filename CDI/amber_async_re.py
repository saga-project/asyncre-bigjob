import os, sys, random
from pj_async_re import async_re_job
import numpy as np
from amberio.amberrun import ReadAmberGroupfile, AmberRun
import copy, math # only needed for debug

__all__ = 'pj_amber_job, BOLTZMANN_CONSTANT'

# =============================================
# AMBER environment check before doing anything
# =============================================
# Is AMBER installed and an appropriate version?
AMBERHOME = os.getenv('AMBERHOME')
if AMBERHOME is None:
    raise Exception('AMBERHOME is not set.')
sys.path.append(os.path.join(AMBERHOME,'bin'))
try:
    from chemistry.amber.readparm import rst7
except:
    raise Exception('Could not load AMBER python libraries. These are only'
                    ' available in AmberTools12 and later.')

# This definition uses values from sander/src/constants.F90
BOLTZMANN_CONSTANT = 1.380658*6.0221367/4184 # in kcal/mol-K

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        
        self.verbose = False
        if self.keywords.get('VERBOSE') == "yes": self.verbose = True
        
        # ===========================
        # Set up the AMBER executable
        # ===========================
        # Check which AMBER MD engine to use, default to sander
        engine = self.keywords.get('ENGINE').upper()
        sander_flags = [ 'AMBER', 'SANDER', 'AMBER-SANDER' ]
        pmemd_flags = [ 'PMEMD', 'AMBER-PMEMD' ]
        engine_name = ''
        if engine in sander_flags:  
            engine_name = 'sander'
        elif engine in pmemd_flags: 
            engine_name = 'pmemd'
        else:
            self._exit('ENGINE is not AMBER')
        if self.spmd == 'mpi' or int(self.keywords.get('SUBJOB_CORES')) > 1:
            self.spmd = 'mpi'
            engine_name += '.MPI'
        # else just assume that a serial executable is desired
        # TODO?: Cuda

        # Check that this executable exists, etc.
        self.exe = os.path.join(AMBERHOME,'bin',engine_name)
        if not os.path.exists(self.exe) or not os.access(self.exe,os.X_OK):
            self._exit('Could not find an executable: %s\nExpected it to'
                       ' be at %s'%(engine_name,self.exe))

        # ========================================================
        # Set up the general state/replica information - 2 methods
        # ========================================================
        # (1) If present, read the AMBER groupfile and define the states,
        if self.keywords.get('AMBER_GROUPFILE') is not None:
            groupfile = self.keywords.get('AMBER_GROUPFILE')
            self.states = ReadAmberGroupfile(groupfile)
            self.states.SetEngine(engine_name)
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
                     'inpcrd' : ['.inpcrd','.rst7','.mdcrd','.crd'] }
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
            for file in exts.keys():
                if files[file] is None:
                    self._exit('Could not identify a(n) %s file based on the'
                               ' given input!'%file)
            # Define states with these file names
            if self.verbose:
                print ('Creating %d replicas using the provided'
                       ' ENGINE_INPUT_EXTFILES'%self.nreplicas)
            self.states = [ AmberRun(mode='-O',basename=self.basename,
                                     engine=engine_name,**files) 
                            for n in range(self.nreplicas) ]

    def _buildInpFile(self, replica):
        """
        For a given replica, determine the state and Write a new AMBER mdin 
        file and link to new prmtop and ref files as needed (this uniquely
        defines a state in AMBER). If this is the first cycle for this replica,
        also make sure that input coordinates (cycle=0) are available.
        """
        # convenient shorthand for current information
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']
        # Write a new input file (make sure it is a restart for cycle > 1)
        input_file = 'r%d/mdin'%replica
        if cycle > 1:  self.states[stateid].mdin.SetRestart()
        self.states[stateid].mdin.WriteAmberMdinFile(input_file)
        # Link to a new prmtop
        self._linkReplicaFile('prmtop',
                              self.states[stateid].filenames['prmtop'],
                              replica)
        # If restraints with reference coordinates were specificed, then update
        # that coordinate file as well.
        if self.states[stateid].mdin.GetVariableValue('ntr','cntrl') == 1:
            self._linkReplicaFile('refc',
                                  self.states[stateid].filenames['ref'],
                                  replica)
        # If this is the first cycle then also link to input coordinates
        if cycle == 1:
            self._linkReplicaFile('%s_0.rst7'%self.basename,
                                  self.states[stateid].filenames['inpcrd'],
                                  replica)

    def _launchReplica(self,replica,cycle):
        """
        Launch an AMBER sub-job using pilot-job. 

        The input files for AMBER that define a state are assumed to be the 
        default names mdin, prmtop, and refc. These files are always re-written 
        or symlinked to in _buildInpFile().
        """
        # Working directory for this replica
        wdir = '%s/r%d'%(os.getcwd(),replica)

        # Cycle dependent input and output file names
        inpcrd = '%s_%d.rst7'%(self.basename,cycle-1)
        mdout  = '%s_%d.out'%(self.basename,cycle)
        mdcrd  = '%s_%d.nc'%(self.basename,cycle)
        rstrt  = '%s_%d.rst7'%(self.basename,cycle)
        stdout = '%s_%d.log'%(self.basename,cycle)
        stderr = '%s_%d.err'%(self.basename,cycle)

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

        if self.verbose:
            engine_name = self.exe.split('/')[-1]
            print 'Launching %s in %s (cycle %d)'%(engine_name,wdir,cycle)

        compute_unit = self.pilotcompute.submit_compute_unit(cpt_unit_desc)
        return compute_unit
        
    def _hasCompleted(self,replica,cycle):
        """
        Returns true if an Amber replica has completed a cycle. Basically 
        checks if the restart file exists.
        """
        # TODO: Parse the output file and look for more sure signs of 
        #       completion?
        rst = 'r%d/%s_%d.rst7'%(replica,self.basename,cycle)
        if os.path.exists(rst):
            return True
        else:
            return False

    def _linkReplicaFile(self, link_filename, real_filename, replica):
        """
        Link the file at real_filename to the name at link_filename in the
        directory belonging to the given replica. If a file is already linked
        to this name (e.g. from a previous cycle), remove it first.
        """
        # TODO: move this to pj_async_re?
        real_filename = '../%s'%real_filename
        os.chdir('r%d'%replica)
        if os.path.exists(link_filename): os.remove(link_filename)
        os.symlink(real_filename,link_filename)
        os.chdir('..')

    def _extractLastCoordinates(self,repl):
        """
        Returns a 3N list of coordinates from the last restart (rst7) file of a
        given replica.
        """
        cycle = self.status[repl]['cycle_current']
        rst = 'r%d/%s_%d.rst7'%(repl,self.basename,cycle)
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
   
    def _gibbs_re_j(self,repl_i,replicas,U):
        # Select a replica swap partner based off independence sampling
        #
        # Pswap = exp(-duij) / sum(exp(-duij))
        ee = U
        ps = np.zeros(len(replicas))
        for repl_j in replicas: 
            ps[j] = ee[i][i] + ee[j][j] - ee[i][j] - e[j][i]
        ps = np.exp(-ps)
        ps /= ps.sum()
        return replicas[self._weighted_choice_sub(ps)]
      
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

    def _swapStates(self,repl_a,repl_b):
        sid_a = self.status[repl_a]['stateid_current']
        sid_b = self.status[repl_b]['stateid_current']
        self.status[repl_a]['stateid_current'] = sid_b
        self.status[repl_b]['stateid_current'] = sid_a

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
        nreps = nwaiting**4
        npermt = {}
        permt = {}
        U = self._computeSwapMatrix(replicas_waiting,states_waiting)
        # Remember that U is nreplicas x nreplicas with rows corresponding
        # to the replica ids and columns corresponding to the STATIC state
        # ids (state id and replica id can and will be different)
        for reps in range(nreps):
            #
            # Independence sampling type swaps
            #
            for repl_i in replicas_waiting:
                sid_i = self.status[repl_i]['stateid_current']
                ps = np.zeros(nwaiting)
                for j,repl_j in zip(range(nwaiting),replicas_waiting):
                    sid_j = self.status[repl_j]['stateid_current']
                    ps[j] = ( U[sid_i][repl_i] + U[sid_j][repl_j] 
                              - U[sid_i][repl_j] - U[sid_j][repl_i] )
                    ps = np.exp(-ps)
                    ps /= ps.sum()
                    repl_j = replicas_waiting[self._weighted_choice_sub(ps)]
                    if repl_i != repl_j:
                        self._swapStates(repl_i,repl_j)
            #
            # Traditional Metropolis type swaps
            #
            # repl_i,repl_j = random.sample(replicas_waiting,2)
            # sid_i = self.status[repl_i]['stateid_current']
            # sid_j = self.status[repl_j]['stateid_current']
            # exchange = pj_amber_job.ReplicaExchange(U[sid_i][repl_i],
            #                                         U[sid_j][repl_j],
            #                                         U[sid_i][repl_j],
            #                                         U[sid_j][repl_i])
            # if exchange: self._swapStates(repl_i,repl_j)
        ###### DEBUG
            # list of current states of ALL replicas
            curr_states = [ self.status[i]['stateid_current'] 
                            for i in replicas_waiting ]
            # list of tuples of replicas and stateids
            # e.g. '[(0,1), (1,2), (2,0)]'
            curr_perm = str(zip(replicas_waiting,curr_states))
           
            # If the permutation has been observed, add a count 
            if npermt.has_key(curr_perm):
                npermt[curr_perm] += 1
            # otherwise add this permutation and store it as a possibility
            else:
                npermt[curr_perm] = 1
                permt[curr_perm] = copy.copy(curr_states)

        print ('Report for %d rounds of independence sampling of %d'
               ' replicas'%(nreps,len(replicas_waiting)))
        print ('Swaps among replicas %s in states %s N! = %d permutations'
               %(str(replicas_waiting),str(curr_states),
                 math.factorial(len(replicas_waiting))))

        from itertools import permutations
        emp  = []
        exact = []
        sumps = 0
        ss = 0
        for state_perm in permutations(curr_states):
            perm = str(zip(replicas_waiting,state_perm))
            # emperical distribution observed here
            if npermt.has_key(perm):
                emp.append(npermt[perm])
                ss += npermt[perm]
            else:
                emp.append(0.)
            # exact boltzmann weight of all permutations
            e = 0
            for i,j in zip(replicas_waiting,state_perm):
                e += self._reduced_energy(j,i)
            p = math.exp(-e)
            exact.append(p)
            sumps += p
        exact = [ p / sumps for p in exact ]
        emp   = [ float(p) / ss for p in emp ]
        sum1 = 0.
        sum2 = 0.
        DKL = 0.
        print '   empirical exact   %diff state permutation'
        print '--------------------------------------'
        for k,state_perm in enumerate(permutations(curr_states)):
            perm = str(zip(replicas_waiting,state_perm))
            pdiff = (emp[k]-exact[k])*100/exact[k]
            print '%2d %8.6f %8.6f %6.1f %s'%(k+1,emp[k],exact[k],pdiff,perm)
            sum1 += emp[k]
            sum2 += exact[k]
            if emp[k] > 0.: DKL += emp[k]*math.log(emp[k]/exact[k])
        print '--------------------------------------'
        print ('   %8.6f %8.6f (sum) Kullback-Liebler Divergence = %f'
               %(sum1,sum2,DKL))
        #######

        # write input files
        for k in replicas_waiting:
            # Creates new input file for the next cycle
            # Places replica back into "W" (wait) state 
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = "W"
            
        return 1
