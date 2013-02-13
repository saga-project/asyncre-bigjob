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

        # ============================================
        # Set up the general state/replica information
        # ============================================
        # If present, read the AMBER groupfile and define the states,
        if self.keywords.get('AMBER_GROUPFILE') is not None:
            groupfile = self.keywords.get('AMBER_GROUPFILE')
            self.states = ReadAmberGroupfile(groupfile)
            self.nreplicas = len(self.states)
            if self.keywords.get('VERBOSE') == 'yes':
                print ('Creating %d replicas from AMBER groupfile: %s'
                       %(self.nreplicas,groupfile))
        # otherwise assume that the states can be inferred from the extfiles
        # and input from a specific application (e.g. umbrella sampling).
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
            if self.keywords.get('VERBOSE') == 'yes':
                print ('Creating %d replicas using the provided'
                       ' ENGINE_INPUT_EXTFILES'%self.nreplicas)
            self.states = [ AmberRun(**files) for n in range(self.nreplicas) ]
        # Set the inputfile defaults based on the MD engine
        for state in self.states: state.mdin.SetDefaults(engine_name)

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
        """Launch an AMBER sub-job using pilot-job. The input files for AMBER
        that define a state are assumed to be the default names mdin, prmtop, 
        and refc. Those files are simply re-written or symlinked to.

        (see _buildInpFile)
        """
        stateid = self.status[replica]['stateid_current']
        self.states[stateid].SetBasename('%s_%d'%(self.basename,cycle))
        self.states[stateid].filenames['inpcrd'] = '%s_%d.rst7'%(self.basename,
                                                                 cycle-1)    
        args = self.states[stateid].GetArgs()
        # delete specific names for mdin, prmtop, and ref files, these are 
        # always re-written or else linked to by the default names.
        for flag in ['-i','-p','-ref']:
            try:
                i = args.index(flag)
                args.pop(i)
                args.pop(i)
            except ValueError:
                pass

        log_file = '%s_%d.log'%(self.basename,cycle)
        err_file = '%s_%d.err'%(self.basename,cycle)

        #pilotjob: Compute Unit (i.e. Job) description
        compute_unit_description = {
            "executable": self.exe,
            "environment": [],
            "arguments": args,
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

        compute_unit=self.pilotcompute.submit_compute_unit(
            compute_unit_description)
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

###########################################################################
#
# Work in Progress: Gibbs sampling style exchanges (see impact_async_re.py)
#
###########################################################################
    def _weighted_choice_sub(self,weights):
        # gives random choice from a set with weighted probabilities
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
   
    def _gibbs_re_j(self,i,replicas):
        # Select a replica swap partner based off independence sampling
        ee = np.asarray(self._computeSwapMatrix(replicas),'float64')
        ps = np.zeros(len(replicas))
        for j in range(len(replicas)): ps[j] = -ee[i][j] + ee[i][i] + ee[j][j]
        min = ps.min()
        ps -= min
        ps = np.exp(ps)*np.exp(min)
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

        nreps = 10
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
            curr_states = [ self.status[i]['stateid_current'] 
                            for i in range(self.nreplicas) ]
            curr_perm = str(zip(range(self.nreplicas),curr_states))
            # e.g. '[(0,1), (1,2), (2,0)]'

            if npermt.has_key(curr_perm):
                npermt[curr_perm] += 1
            else:
                npermt[curr_perm] = 1
                permt[curr_perm] = copy.copy(curr_states)

        print ('Report for %d rounds of independence sampling of %d'
               ' replicas'%(nreps,len(replicas_waiting)))
        print 'Swaps among replicas',replicas_waiting
        #for k in replicas_waiting:
        #    state = self.status[k]['stateid_current']
        #    print "Replica %d is now in state %d"%(k,state)
        #    self.umbrellas[state].PrintRestraintReport()

        ss = 0 # total number of non-unique permutations
        for k in npermt.keys(): ss += npermt[k]
        ps = []
        sumps = 0
        for k in npermt.keys():
            b = permt[k] # the state list of all replicas
            e = 0
            for i in replicas_waiting: e += self._reduced_energy(b[i],i)
            p = math.exp(-e)
            ps.append(p)
            sumps += p
        i = 0
        sum1 = 0.
        sum2 = 0.
        DKL = 0.
        print 'empirical exact   %diff state permutation'
        print '-----------------------------------'
        for k in npermt.keys():
            emp = npermt[k] / float(ss)
            exa = ps[i]/sumps
            pdiff = (emp-exa)*100/exa
            print '%8.6f %8.6f %5.1f %s'%(emp,exa,pdiff,str(k))
            i += 1
            sum1 += emp
            sum2 += exa
            DKL += emp*math.log(emp/exa)
        print '-----------------------------------'
        print '%8.6f %8.6f (sum)'%(sum1,sum2)
        print 'Kullback-Liebler Divergence =',DKL
        #######

        # write input files
        for k in replicas_waiting:
            # Creates new input file for the next cycle
            # Places replica back into "W" (wait) state 
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = "W"
