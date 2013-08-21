# File Based Replica Exchange class
"""A module to prepare and run file-based asynchronous RE jobs
See documentation in doc/ directory.

Contributors: 

Emilio Gallicchio <emilio@biomaps.rutgers.edu>
Brian Radak <radakb@biomaps.rutgers.edu>
Melissa Romanus <melissa.romanus@rutgers.edu>
"""
import os
import sys
import time
import pickle
import random
import copy

from numpy import *
from configobj import ConfigObj

from pilot import PilotComputeService, ComputeDataService, State

__version__ = '0.2.1'

class async_re_job(object):
    """
    Class to set up and run asynchronous file-based RE calculations
    """
    def __init__(self, command_file, options):
        self.command_file = command_file
        self.cus = {}
        self.jobname = os.path.splitext(os.path.basename(command_file))[0]
        self._parseInputFile()
        self._checkInput()
        self._printStatus()

    def __getattribute__(self, name):
        if name == 'replicas_waiting':
            # Return a list of replica indices of replicas in a wait state.
            self.updateStatus()
            return [k for k in range(self.nreplicas) 
                    if self.status[k]['running_status'] == 'W']
        elif name == 'states_waiting':
            # Return a list of state ids of replicas in a wait state.
            return [self.status[k]['stateid_current'] 
                    for k in self.replicas_waiting]
        elif name == 'replicas_waiting_to_exchange':
            # Return a list of replica indices of replicas in a wait state that
            # have ALSO completed at least one cycle.
            self.updateStatus()
            return [k for k in range(self.nreplicas) 
                    if (self.status[k]['running_status'] == 'W' and
                        self.status[k]['cycle_current'] > 1)]
        elif name == 'states_waiting_to_exchange':
            # Return a list of state ids of replicas in a wait state that have 
            # ALSO completed at least one cycle.
            return [self.status[k]['stateid_current'] 
                    for k in self.replicas_waiting_to_exchange]
        else:
            return object.__getattribute__(self,name)

    def _error(self, text):
        """Print and flush an error message to stdout."""
        print text
        sys.stdout.flush()
    
    def _exit(self, text):
        """Print an error message and exit."""
        self._error(text)
        print 'exiting...'
        sys.exit(1)

    def _parseInputFile(self):
        """Read keywords from a configure file."""
        self.keywords = ConfigObj(self.command_file)
        
    def _printStatus(self):
        """Print a report of the input parameters."""
        print 'command_file =',self.command_file
        print 'jobname =',self.jobname
        for k,v in self.keywords.iteritems():
            print k,v

    def _openfile(self, name, mode, max_attempts = 100):
        attempts = 0
        f = None
        while not f and attempts <= max_attempts:
            try:
                f = open(name,mode)
            except IOError:
                print ('Warning: unable to access file %s, re-trying in 1 '
                       'second ...')%name
                f = None
                attempts += 1
                time.sleep(1)
        if attempts > max_attempts:
            self._exit('Too many failures accessing file %s: quitting.'%name)
        return f

    def _checkInput(self):
        """ 
        Check that required parameters are specified. Parse these and other 
        optional settings. 
        """ 
        # Required Options
        #
        # basename for the job
        self.basename = self.keywords.get('ENGINE_INPUT_BASENAME')
        if self.basename is None:
            self._exit('ENGINE_INPUT_BASENAME needs to be specified')
        # execution time in minutes
        self.walltime = float(self.keywords.get('WALL_TIME'))
        if self.walltime is None:
            self._exit('WALL_TIME needs to be specified')
        # variables required for PilotJob
        if self.keywords.get('COORDINATION_URL') is None:
            self._exit('COORDINATION_URL needs to be specified')
        if self.keywords.get('RESOURCE_URL') is None:
            self._exit('RESOURCE_URL needs to be specified')
        if self.keywords.get('QUEUE') is None:
            self._exit('QUEUE needs to be specified')
        if self.keywords.get('BJ_WORKING_DIR') is None:
            self._exit('BJ_WORKING_DIR needs to be specified')
        if self.keywords.get('TOTAL_CORES') is None:
            self._exit('TOTAL_CORES needs to be specified')
        if self.keywords.get('SUBJOB_CORES') is None:
            self._exit('SUBJOB_CORES needs to be specified')

        # Optional variables
        #
        # processors per node on this machine (can be auto-detected)
        self.ppn = 1
        if self.keywords.get('PPN') is not None: 
            self.ppn = int(self.keywords.get('PPN'))
        # spmd_variation for PilotJob (may override this later)
        self.spmd = 'single'
        if self.keywords.get('SPMD') is not None: 
            self.spmd = self.keywords.get('SPMD')
        # number of replicas (may be determined by other means)
        self.nreplicas = None
        
        self.nexchg_rounds = 1
        if self.keywords.get('NEXCHG_ROUNDS') is not None:
            self.nexchg_rounds = int(self.keywords.get('NEXCHG_ROUNDS'))

        #examine RESOURCE_URL to see if it's remote (file staging)
#        self.remote = self._check_remote_resource(self.keywords.get('RESOURCE_URL'))
#        if self.remote:
#            print "Use remote execution and file staging"
#            if self.keywords.get('REMOTE_WORKING_DIR') is None:
#                self._exit("REMOTE_WORKING_DIR needs to be specified")
#            if self.keywords.get('REMOTE_DATA_SERVICE') is None: #something like ssh://<user>@<machine>/<datadir>
#                self._exit("REMOTE_DATA_SERVICE needs to be specified")
#            if self.keywords.get('REMOTE_DATA_SIZE') is None:
#                self.remote_data_size = 2048 # 2GB by default
#            else:
#                self.remote_data_size = self.keywords.get('REMOTE_DATA_SIZE')

        if self.keywords.get('NREPLICAS') is not None:
            self.nreplicas = int(self.keywords.get('NREPLICAS'))
        # extfiles variable for 'setupJob'
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if self.extfiles is not None and self.extfiles != '':
            self.extfiles = self.extfiles.split(',')
        else:
            self.extfiles = None
        # verbose printing
        self.verbose = False
        if self.keywords.get('VERBOSE').lower() == 'yes': 
            self.verbose = True

    def _linkReplicaFile(self, link_filename, real_filename, repl):
        """
        Link the file at real_filename to the name at link_filename in the
        directory belonging to the given replica. If a file is already linked
        to this name (e.g. from a previous cycle), remove it first.
        """
        os.chdir('r%d'%repl)
        # Check that the file to be linked actually exists.
        # TODO: This is not robust to absolute path specifications.
        real_filename = '../%s'%real_filename
        if not os.path.exists(real_filename):
            self._exit('No such file: %s'%real_filename)
        # Make/re-make the symlink.
        if os.path.exists(link_filename): 
            os.remove(link_filename)
        os.symlink(real_filename,link_filename)
        os.chdir('..')

    def setupJob(self):
        """
        If RE_SETUP='yes' creates and populates subdirectories, one for each 
        replica called r0, r1, ..., rN in the working directory. Otherwise 
        reads saved state from the ENGINE_BASENAME.stat file.
        
        To populate each directory calls _buildInpFile(k) to prepare the MD 
        engine input file for replica k. Also creates soft links to the working 
        directory for the accessory files specified in ENGINE_INPUT_EXTFILES.
        """
	#pilotjob: Initialize PilotJob at given COORDINATION_URL (CU)
        self.pj = PilotComputeService(self.keywords.get('COORDINATION_URL'))
	#pilotjob: Initialize PilotJob Data service (DU)
        self.cds=ComputeDataService()
	#pilotjob: Launch the PilotJob at the given COORDINATION_URL
        self.launch_pilotjob()

        if (self.keywords.get('RE_SETUP') is not None and 
            self.keywords.get('RE_SETUP').lower() == 'yes'):
            # create replicas directories r1, r2, etc.
            for k in range(self.nreplicas):
                repl_dir = 'r%d'%k
                if os.path.exists(repl_dir):
                    self._exit('Replica directories already exist. Either '
                               'turn off RE_SETUP or remove the directories.')
                else:
                    os.mkdir('r%d'%k)
            # create links for external files
            if self.extfiles is not None:
                for file in self.extfiles:
                    for k in range(self.nreplicas):
                        self._linkReplicaFile(file,file,k)
            # create status table
            self.status = [{'stateid_current': k, 'running_status': 'W', 
                            'cycle_current': 1} for k in range(self.nreplicas)]
            # save status tables
            self._write_status()
            # create input files no. 1
            for k in range(self.nreplicas):
                self._buildInpFile(k)
            self.updateStatus()
        else:
            self._read_status()
            self.updateStatus(restart=True)

#        if self.remote:
#            self._setup_remote_workdir()

        self.print_status()
        #at this point all replicas should be in wait state
        for k in range(self.nreplicas):
            if self.status[k]['running_status'] != 'W':
                self._exit('Internal error after restart. Not all jobs are in '
                           'wait state.')

    def scheduleJobs(self):
        # wait until bigjob enters executing
        while self.pilotcompute.get_state() != 'Running':
            time.sleep(10)

        # Gets the wall clock time for a replica to complete a cycle
        # If unspecified it is estimated as 10% of job wall clock time
        # Note  
        replica_run_time = self.keywords.get('REPLICA_RUN_TIME')
        if self.keywords.get('REPLICA_RUN_TIME') is None:
            replica_run_time = int(round(self.walltime/10.))
        else:
            replica_run_time = int(self.keywords.get('REPLICA_RUN_TIME'))
        # double it to give time for current running processes 
        # and newly submitted processes to complete
        replica_run_time *= 2

        # Time in between cycles in seconds
        # If unspecified it is set as 30 secs
        if self.keywords.get('CYCLE_TIME') is None:
            cycle_time = 30.0
        else:
            cycle_time = float(self.keywords.get('CYCLE_TIME'))

        start_time = time.time()
        end_time = (start_time + 60*(self.walltime - replica_run_time) - 
                    cycle_time - 10)
        while time.time() < end_time:
            time.sleep(1)

            self.updateStatus()
            self.print_status()
            self.launchJobs()
            self.updateStatus()
            self.print_status()        

            time.sleep(cycle_time)

            self.updateStatus()
            self.print_status()        
            self.doExchanges()
        
        self.updateStatus()
        self.print_status()
        self.waitJob()
        self.cleanJob()

    def waitJob(self):
        # cancel all not-running submitted subjobs
#        for k in range(self.nreplicas):
#            if self.status[k]['running_status'] == "R":
#                if self.cus[k].get_state() != "Running":
#                    self.cus[k].cancel()
#                    self.status[k]['running_status'] = "W"
        #update status
#        self.updateStatus()
#        self.print_status()
        #wait until running jobs complete
        self.cds.wait()

    def cleanJob(self):
        self.cds.cancel()
        self.pj.cancel()
        
    def launch_pilotjob(self):
	#pilotjob: PilotJob description
	#pilotjob: Variables defined in command.inp
        pcd = {'service_url': self.keywords.get('RESOURCE_URL'),
               'number_of_processes': self.keywords.get('TOTAL_CORES'),
               'working_directory': self.keywords.get('BJ_WORKING_DIR'),
               'queue': self.keywords.get('QUEUE'),
               'processes_per_node': self.ppn,
               'project': self.keywords.get('PROJECT'),
               'walltime': int(self.keywords.get('WALL_TIME'))}

        if self.keywords.get('SGE_WAYNESS') is not None:
                pcd['spmd_variation'] = self.keywords.get('SGE_WAYNESS')
         
	#pilotjob: Create pilot job with above description
        self.pj.create_pilot(pilot_compute_description=pcd)
        self.cds.add_pilot_compute_service(self.pj)
        self.pilotcompute = self.pj.list_pilots()[0]
        
    def _write_status(self):
        """
        Saves the current state of the RE job in the BASENAME.stat file using 
        pickle
        """
        status_file = '%s.stat'%self.basename
        f = self._openfile(status_file,'w')
        pickle.dump(self.status,f)
        #pickle.dump(self.node_status, f)
        f.close()

    def _read_status(self):
        """
        Loads the current state of the RE job from BASENAME.stat file using 
        pickle
        """
        status_file = '%s.stat'%self.basename
        f = self._openfile(status_file,'r')
        self.status = pickle.load(f)
        f.close()

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job. 
        It's fun to follow the progress in real time by doing:
        watch cat BASENAME_stat.txt
        """
        log = 'Replica  State  Status  Cycle \n'
        for k in range(self.nreplicas):
            log += ('%6d   %5d  %5s  %5d \n'%
                    (k,self.status[k]['stateid_current'], 
                     self.status[k]['running_status'],
                     self.status[k]['cycle_current']))                    
        log += 'Running = %d\n'%self.running
        log += 'Waiting = %d\n'%self.waiting

        logfile = '%s_stat.txt'%self.basename
        ofile = self._openfile(logfile,'w')
        ofile.write(log)
        ofile.close()

    def updateStatus(self, restart = False):
        """Scan the replicas and update their states."""
        for k in range(self.nreplicas):
            self._updateStatus_replica(k,restart)
        self._write_status()
        self._update_running_no()

    def _updateStatus_replica(self, replica, restart):
        """
        Update the status of the specified replica. If it has completed a cycle
        the input file for the next cycle is prepared and the replica is placed
        in the wait state.
        """
        this_cycle = self.status[replica]['cycle_current']
        if restart:
            if self.status[replica]['running_status'] == 'R':
                if self._hasCompleted(replica,this_cycle):
                    self.status[replica]['cycle_current'] += 1
                else:
                    print ('_updateStatus_replica(): Warning: restarting '
                           'replica %d (cycle %d)'%(replica,this_cycle))
            self._buildInpFile(replica)
            self.status[replica]['running_status'] = 'W'
        else:
            if self.status[replica]['running_status'] == 'R':
                if self._isDone(replica,this_cycle):
                    self.status[replica]['running_status'] = 'S'
                    if self._hasCompleted(replica,this_cycle):
                        self.status[replica]['cycle_current'] += 1
                    else:
                        print ('_updateStatus_replica(): Warning: restarting '
                               'replica %d (cycle %d)'%(replica,this_cycle))
                    self._buildInpFile(replica)
                    self.status[replica]['running_status'] = 'W'
                            
    def _update_running_no(self):
        """Update the number of running and waiting replicas."""
        self.running = 0
        self.waiting = 0
        for k in range(self.nreplicas):
            if self.status[k]['running_status'] == 'R':
                self.running += 1
            if self.status[k]['running_status'] == 'W':
                self.waiting += 1

    def _isDone(self,replica,cycle):
        """
        Generic function to check if a replica completed a cycle. 
        Calls in this case pilot-job version.
        """
        return self._isDone_PJ(replica,cycle)

    def _isDone_PJ(self,replica,cycle):
        """Return true if a replica has exited (done or failed)."""
       	#pilotjob: Get status of the compute unit
	#pilotjob: Query the replica to see if it is in the done state
        state = self.cus[replica].get_state()
        details = self.cus[replica].get_details()
        if state == 'Done' or state == 'Failed' or state == 'Canceled':
            if self.verbose:
                if details.has_key('start_time'):
                    if details.has_key('end_time'):
                        print '*'*80
                        print ('Replica: %d Start Time: %f End Time: %f'%
                               (replica,float(details['start_time']),
                                float(details['end_time'])))
                if details.has_key('end_queue_time'):
                    print ('End Queue Time: %f\n'%
                           float(details['end_queue_time']))
            return True
        else:
            return False
            
    def _hasCompleted(self,replica,cycle):
        """
        Attempts to check whether a replica has completed successfully from the
        bigjob compute unit. This is not expected to work during a restart when
        compute units are not available. In the latter case success is assumed. 
        MD engine modules are recommended to override this default routine with 
        one that implements a better test of success such as the existence of a 
        restart file or similar.
        """ 
        try:
            state = self.cus[replica].get_state()
        except:
            print ('_hasCompleted(): Warning: unable to query replica state. '
                   'Assuming success ...')
            return True

        if state == 'Done':
            return True
        else:
            return False

    def _njobs_to_run(self):
        # size of subjob buffer as a percentage of job slots 
        # (TOTAL_CORES/SUBJOB_CORES)
        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is None:
            subjobs_buffer_size = 0.5
        else:
            subjobs_buffer_size = float(subjobs_buffer_size)
        # find out how many replicas are waiting/(running/submitted)
        self._update_running_no()
        # launch new replicas if the number of submitted/running subjobs is 
        # less than the number of available slots 
        # (total_cores/subjob_cores) + 50%
        available_slots = (int(self.keywords.get('TOTAL_CORES')) / 
                           int(self.keywords.get('SUBJOB_CORES')))
        max_njobs_submitted = int((1.+subjobs_buffer_size)*available_slots)
        nlaunch = self.waiting - max(2,self.nreplicas - max_njobs_submitted)
        nlaunch = max(0,nlaunch)
        if self.verbose:
            print 'available_slots: %d'%available_slots
            print 'max_njobs_submitted: %d'%max_njobs_submitted
            print 'running/submitted subjobs: %d'%self.running
            print 'waiting replicas: %d'%self.waiting
            print 'replicas to launch: %d'%nlaunch
        return nlaunch

    def launchJobs(self):
        """
        Scans the replicas in wait state and randomly launches some of them
        if CPU's are available.
        """ 
        jobs_to_launch = self._njobs_to_run()
        if jobs_to_launch > 0:
            wait = [k for k in range(self.nreplicas) 
                    if self.status[k]['running_status'] == 'W']
            random.shuffle(wait)
            n = min(jobs_to_launch,len(wait))
            for k in wait[0:n]:
                if self.verbose:
                    print ('Launching replica %d cycle %d'
                           %(k,self.status[k]['cycle_current']))
                self.cus[k] = (
                    self._launchReplica(k,self.status[k]['cycle_current']))
                self.status[k]['running_status'] = 'R'

# gives random choice from a set with weight probabilities
    def _weighted_choice_sub(self,weights):
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
                
# _gibbs_re_j() produces a replica "j" to exchange with the given replica "i"
# based on independent sampling from the discrete Metropolis transition matrix
#
# T_rs = alpha_rs min[1,exp(-du_rs)] ; r not= s
# T_rr = 1 - sum_(s not= r) T_rs 
#
# where r and s are replica exchange permutations, r being the current
# permutation and s the new permutation. alpha_rs = 0 unless permutations
# r and s differ by a single replica swap and alpha_rs = 1/(n-1) otherwise,
# n being the number of replicas and (n-1) is the number of permutations s
# differing by permutation r by a single swap. du_rs is the change in
# reduced potential energy of the replica exchange ensemble in going from
# permutations r to permutation s (that is due to a replica swap).
# Based on the above we have
# du_rs = u_a(j)+u_b(i)-[u_a(i)+u_b(j)]
# where i and j are the replica being swapped and a and b, respectively, are the 
# states they occupy in the r permutations and b and a, respectively, those in
# the s permutations.
# 
# The energies u_a(i), i=1,n and a=1,n, are assumed stored in the input matrix U[a][i].
#
# In general, the set of replicas across which exchanges are considered is
# a subset of the n replicas. This list is passed in the 'replicas_waiting'
# list. Replica i ('repl_i') is assumed to be in this list.
#
    def _gibbs_re_j(self, repl_i, replicas_waiting, U):
        n = len(replicas_waiting)
        if n < 2:
            return repl_i
        #evaluate all i-j swap probabilities
        ps = zeros(n)
        du = zeros(n)
        eu = zeros(n)
        #
        sid_i = self.status[repl_i]['stateid_current'] 
        for j in range(n):
            repl_j = replicas_waiting[j]
            sid_j = self.status[repl_j]['stateid_current']
            du[j] = (U[sid_i][repl_j] + U[sid_j][repl_i] 
                  - U[sid_i][repl_i] - U[sid_j][repl_j])
        eu = exp(-du)
        #
        pii = 1.0
        i = -1
        f = 1./(float(n) - 1.)
        for j in range(n):
            repl_j = replicas_waiting[j]
            if repl_j == repl_i:
                i = j
            else:
                if eu[j] > 1.0:
                    ps[j] = f
                else:                    
                    ps[j] = f*eu[j]
                pii -= ps[j]
        try:
            ps[i] = pii
        except:
            self._exit('gibbs_re_j(): unrecoverable error: replica i not in the list of waiting replicas?')
        #index of swap replica within replicas_waiting list
        j = self._weighted_choice_sub(ps)
        #actual replica
        repl_j = replicas_waiting[j]
        return repl_j

    def doExchanges(self):
        """
        Perform n rounds of exchanges among waiting replicas using Gibbs 
        sampling.
        """
        print 'Entering doExchanges Method: %f'%time.time()

        # find out which replicas are waiting
        replicas_waiting = self.replicas_waiting_to_exchange
        states_waiting = self.states_waiting_to_exchange
        if len(replicas_waiting) > 1:
            # backtrack cycle
            for k in replicas_waiting:
                self.status[k]['cycle_current'] -= 1
                self.status[k]['running_status'] = 'E'
            # Matrix of replica energies in each state.
            # The computeSwapMatrix() function is defined by application 
            # classes (Amber/US, Impact/BEDAM, etc.)
            U = self._computeSwapMatrix(replicas_waiting,states_waiting)
            # Perform an exchange for each of the n replicas, m times
            # Perform an exchange for each of the n replicas, m times
            if self.nexchg_rounds >= 0:
                mreps = self.nexchg_rounds
            else:
                mreps = len(replicas_waiting)**(-self.nexchg_rounds)
            for reps in range(mreps):
                for repl_i in replicas_waiting:
                    repl_j = self._gibbs_re_j(repl_i,replicas_waiting,U)
                    if repl_j != repl_i:
                        #Swap state id's
                        #Note that the energy matrix does not change
                        sid_i = self.status[repl_i]['stateid_current'] 
                        sid_j = self.status[repl_j]['stateid_current']
                        self.status[repl_i]['stateid_current'] = sid_j
                        self.status[repl_j]['stateid_current'] = sid_i

# Uncomment to debug Gibbs sampling: actual and computed populations of 
# state permutations should match
# 
#                     self._debug_collect_state_populations(replicas_waiting,U)
#             self._debug_validate_state_populations(replicas_waiting,U)

            # write input files
            for k in replicas_waiting:
                # Creates new input file for the next cycle
                # Places replica back into "W" (wait) state 
                self.status[k]['cycle_current'] += 1
                self._buildInpFile(k)
                self.status[k]['running_status'] = 'W'

        print 'Exiting doExchanges Method: %f'%time.time()


    def _check_remote_resource(self, resource_url):
        """
        check if it's a remote resource. Basically see if 'ssh' is present
        """
        ssh_c = re.compile("(.+)\+ssh://(.*)")
        m = re.match(ssh_c, resource_url)
        if m:
            self.remote_protocol = m.group(1)
            self.remote_server = m.group(2)
            print resource_url + " : yes" + " " + remote_protocol + " " + remote_server
            return 1
        else:
            print resource_url + " : no"
            return 0

    def _setup_remote_workdir(self):
        """
        rsync local working directory with remote working directory
        """
        os.system("ssh %s mkdir -p %s" % (self.remote_server, self.keywords.get('REMOTE_WORKING_DIR')))
        extfiles = " "
        for efile in self.extfiles:
            extfiles = extfiles + " " + efile
        os.system("rsync -av %s %s/%s/" % (extfiles, self.remote_server, self.keywords.get('REMOTE_WORKING_DIR')))

        dirs = ""
        for k in range(self.nreplicas):
            dirs = dirs + " r%d" % k
        setup_script = """
cd %s ; \
for i in `seq 0 %d` ; do \
mkdir -p r$i ; \
 




"""
      
    def _debug_collect_state_populations(self, replicas, U):
        """
        Calculate the empirically observed distribution of state 
        permutations. Permutations not observed will NOT be counted.
        """
        try:
            self.npermt
        except (NameError,AttributeError):
            self.npermt = {}
  
        try:
            self.permt
        except (NameError,AttributeError):
            self.permt = {}

        curr_states = [self.status[i]['stateid_current'] for i in replicas]
        curr_perm = str(zip(replicas,curr_states))
        if self.npermt.has_key(curr_perm):
            self.npermt[curr_perm] += 1
        else:
            self.npermt[curr_perm] = 1
            self.permt[curr_perm] = copy(curr_states)
        
        #state id list
        # sids = [ self.status[k]['stateid_current'] for k in replicas]

        # try:
        #     self=perme0            
        # except:
        #     e = 0
        #     i = 0
        #     for repl in replicas:
        #         e += U[sids[i]][repl]
        #         i += 1
        #     self.perme0 = e

        # key = str(sids)
        # if self.npermt.has_key(key):
        #     self.npermt[key] += 1
        # else:
        #     self.npermt[key] = 1
        #     self.permt[key] = sids[:]

    def _debug_validate_state_populations(self, replicas, U):
        """
        Calculate the exact state distribution of all possible state
        permutations using the calculated energies. Compare this to the 
        empirical distribution using the Kullback-Liebler divergence:

        DKL = sum_k emp_k*ln(emp_k/exact_k)

        where k indexes the states and emp_k and exact_k are the empirical and
        exact densities of state k respectively. For numerical stability, if
        emp_k is 0 (the state is not observed) it is set to 1e-4. Although 
        this introduces an arbitrary bias, it makes DKL always well-defined.
        """
        from itertools import permutations

        # list of the currently occuplied states
        curr_states = [self.status[i]['stateid_current'] for i in replicas]
        # list of tuples of all possible state permutations
        curr_perm = str(zip(replicas,curr_states))
        print ('Swaps among replicas %s in states %s N! = %d permutations'%
               (str(replicas),str(curr_states),math.factorial(len(replicas))))
        emp  = []
        exact = []
        sumps = 0
        ss = 0
        for state_perm in permutations(curr_states):
            perm = str(zip(replicas,state_perm))
            # emperical distribution observed here
            if self.npermt.has_key(perm):
                emp.append(self.npermt[perm])
                ss += self.npermt[perm]
            else:
                emp.append(0.)
            # exact boltzmann weight of all permutations
            e = 0
            for i,j in zip(replicas,state_perm):
                e += U[j][i]
            p = math.exp(-e)
            exact.append(p)
            sumps += p
        exact = [p/sumps for p in exact]
        emp   = [float(p)/ss for p in emp]
        sum1 = 0.
        sum2 = 0.
        DKL = 0.
        print '%8s %9s %9s %s'%('','empirical','exact','state permutation')
        print '-'*80
        dP = 1.e-9 # this will show up as 0 but contribute a lot to DKL
        for k,state_perm in enumerate(permutations(curr_states)):
            perm = str(zip(replicas,state_perm))
            print '%8d %9.4f %9.4f %s'%(k+1,emp[k],exact[k],perm)
            sum1 += emp[k]
            sum2 += exact[k]
            if emp[k] > 0.: 
                empk = emp[k]
            else:           
                empk = dP
            if exact[k] > 0.: 
                exactk = exact[k]
            else:             
                exactk = dP
            DKL += empk*math.log(empk/exactk)
        print '-'*80
        print ('%8s %9.4f %9.4f (sum) Kullback-Liebler Divergence = %f'
               %('',sum1,sum2,DKL))
        print '='*80
        
        # ss = 0
        # for k in self.npermt.keys():
        #     ss += self.npermt[k]
        # ps = []
        # sumps = 0
        # for k in self.npermt.keys():
        #     sids = self.permt[k]
        #     e = 0
        #     i = 0
        #     for repl in replicas:
        #         e += U[sids[i]][repl]
        #         i += 1
        #     p = math.exp(-(e-self.perme0))
        #     ps.append(p)
        #     sumps += p

        # i = 0
        # for k in self.npermt.keys():
        #     print k, self.npermt[k]/float(ss), ps[i]/sumps
        #     i += 1

        # sys.exit(0)

            
            
