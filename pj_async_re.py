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

from configobj import ConfigObj

from gibbs_sampling import *
from pilot import PilotComputeService, ComputeDataService, State

__version__ = '0.2.1'

def _exit(message):
    """Print and flush a message to stdout and then exit."""
    print message
    sys.stdout.flush()
    print 'exiting...'
    sys.exit(1)

def _open(name, mode, max_attempts = 100, wait_time = 1):
    """
    Convenience function for opening files on an unstable filesystem.
    
    max_attempts : int
        maximum number of attempts to make at opening a file
    wait_time : int
        time (in seconds) to wait between attempts
    """
    attempts = 0
    f = None
    while f is None and attempts <= max_attempts:
        try:
            f = open(name,mode)
        except IOError:
            print ('Warning: unable to access file %s, re-trying in %d '
                   'second(s)...'%(name,wait_time))
            f = None
            attempts += 1
            time.sleep(wait_time)
    if attempts > max_attempts:
        _exit('Too many failures accessing file %s'%name)
    return f

class async_re_job(object):
    """
    Class to set up and run asynchronous file-based RE calculations
    """
    def __init__(self, command_file, options):
        self.command_file = command_file
        self.cus = {}
        self.jobname = os.path.splitext(os.path.basename(command_file))[0]
        self.keywords = ConfigObj(self.command_file)
        self._checkInput()
        self._printStatus()

    def _exit(self, message):
        _exit(message)

    def _openfile(self, name, mode, max_attempts = 100):
        f = _open(name,mode,max_attempts)
        return f

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
        elif name == 'waiting':
            return len(self.replicas_waiting)
        elif name == 'replicas_running':
            # Return a list of replica indices of replicas in a running state.
            self.updateStatus()
            return [k for k in range(self.nreplicas)
                    if self.status[k]['running_status'] == 'R']
        elif name == 'running':
            return len(self.replicas_running)
        else:
            return object.__getattribute__(self,name)
        
    def _printStatus(self):
        """Print a report of the input parameters."""
        print 'command_file =',self.command_file
        print 'jobname =',self.jobname
        for k,v in self.keywords.iteritems():
            print k,v

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
            self._exit('WALL_TIME (in minutes) needs to be specified')
        # variables required for PilotJob
        if self.keywords.get('COORDINATION_URL') is None:
            self._exit('COORDINATION_URL needs to be specified')
        if self.keywords.get('RESOURCE_URL') is None:
            self._exit('RESOURCE_URL needs to be specified')
        if self.keywords.get('QUEUE') is None:
            self._exit('QUEUE needs to be specified')
        if self.keywords.get('BJ_WORKING_DIR') is None:
            basedir = os.getcwd()
        else:
            basedir = self.keywords.get('BJ_WORKING_DIR')
        self.bj_working_dir = os.path.join(basedir,'agent') 
        if not os.path.exists(self.bj_working_dir):
            os.mkdir(self.bj_working_dir)

        if self.keywords.get('TOTAL_CORES') is None:
            self._exit('TOTAL_CORES needs to be specified')
        if self.keywords.get('SUBJOB_CORES') is None:
            self._exit('SUBJOB_CORES needs to be specified')

        # Optional variables
        #
        env = self.keywords.get('ENGINE_ENVIRONMENT')
        if env is not None and env != '':
            self.engine_environment = env.split(',')
        else:
            self.engine_environment = []
        # processors per node on this machine (can be auto-detected)
        if self.keywords.get('PPN') is not None: 
            self.ppn = int(self.keywords.get('PPN'))
        else:
            self.ppn = 1
        # spmd_variation for PilotJob (may override this later)
        if self.keywords.get('SPMD') is not None: 
            self.spmd = self.keywords.get('SPMD')
        else:
            self.spmd = 'single'
        # number of replicas (may be determined by other means)
        self.nreplicas = None
        
        if self.keywords.get('NEXCHG_ROUNDS') is not None:
            self.nexchg_rounds = int(self.keywords.get('NEXCHG_ROUNDS'))
        else:
            self.nexchg_rounds = 1

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
        if self.keywords.get('VERBOSE').lower() == 'yes': 
            self.verbose = True
        else:
            self.verbose = False


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
                    _exit('Replica directories already exist. Either turn off '
                          'RE_SETUP or remove the directories.')
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
                _exit('Internal error after restart. Not all jobs are in wait '
                      'state.')

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
               'working_directory': self.bj_working_dir,
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
        Pickle the current state of the RE job and write to in BASENAME.stat. 
        """
        status_file = '%s.stat'%self.basename
        f = _open(status_file,'w')
        pickle.dump(self.status,f)
        f.close()

    def _read_status(self):
        """
        Unpickle and load the current state of the RE job from BASENAME.stat.
        """
        status_file = '%s.stat'%self.basename
        f = _open(status_file,'r')
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
        ofile = _open(logfile,'w')
        ofile.write(log)
        ofile.close()

    def updateStatus(self, restart = False):
        """Scan the replicas and update their states."""
        for k in range(self.nreplicas):
            self._updateStatus_replica(k,restart)
        self._write_status()

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
        Scan the replicas in wait state and randomly launch some of them if 
        CPU's are available.
        """ 
        jobs_to_launch = self._njobs_to_run()
        if jobs_to_launch > 0:
            wait = self.replicas_waiting
            random.shuffle(wait)
            n = min(jobs_to_launch,len(wait))
            for k in wait[0:n]:
                if self.verbose:
                    print ('Launching replica %d cycle %d'
                           %(k,self.status[k]['cycle_current']))
                self.cus[k] = (
                    self._launchReplica(k,self.status[k]['cycle_current']))
                self.status[k]['running_status'] = 'R'

    def doExchanges(self):
        """Perform exchanges among waiting replicas using Gibbs sampling."""
        # NB: asking for self.replicas_waiting_to_exchange UPDATES the list,
        # therefore this must be kept static at each repetition.
        #
        replicas_to_exchange = self.replicas_waiting_to_exchange
        states_to_exchange = self.states_waiting_to_exchange
        nreplicas_to_exchange = len(replicas_to_exchange)
        if nreplicas_to_exchange < 2:
            return 0

        print 'Initiating exchanges amongst %d replicas:'%nreplicas_to_exchange
        exchange_start_time = time.time()
        # backtrack cycle of waiting replicas
        for k in replicas_to_exchange:
            self.status[k]['cycle_current'] -= 1
            self.status[k]['running_status'] = 'E'
        # Matrix of replica energies in each state.
        # The computeSwapMatrix() function is defined by application 
        # classes (Amber/US, Impact/BEDAM, etc.)
        matrix_start_time = time.time()
        swap_matrix = self._computeSwapMatrix(replicas_to_exchange,
                                              states_to_exchange)
        matrix_time = time.time() - matrix_start_time

        sampling_start_time = time.time()
        # Perform an exchange for each of the n replicas, m times
        if self.nexchg_rounds >= 0:
            mreps = self.nexchg_rounds
        else:
            mreps = nreplicas_to_exchange**(-self.nexchg_rounds)
        for reps in range(mreps):
            for repl_i in replicas_to_exchange:
                sid_i = self.status[repl_i]['stateid_current'] 
                curr_states = [self.status[repl_j]['stateid_current'] 
                               for repl_j in replicas_to_exchange]
                repl_j = pairwise_independence_sampling(repl_i,sid_i,
                                                        replicas_to_exchange,
                                                        curr_states,
                                                        swap_matrix)
                if repl_j != repl_i:
                    sid_i = self.status[repl_i]['stateid_current'] 
                    sid_j = self.status[repl_j]['stateid_current']
                    self.status[repl_i]['stateid_current'] = sid_j
                    self.status[repl_j]['stateid_current'] = sid_i

        # Uncomment to debug Gibbs sampling: 
        # Actual and observed populations of state permutations should match.
        # 
        #     self._debug_collect_state_populations(replicas_to_exchange)
        # self._debug_validate_state_populations(replicas_to_exchange,
        #                                        states_to_exchange,U)
        sampling_time = time.time() - sampling_start_time
        # Write new input files.
        for k in replicas_to_exchange:
            # Create new input files for the next cycle and place replicas back
            # into "W" (wait) state.
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = 'W'

        total_time = time.time() - exchange_start_time

        print '------------------------------------------'
        print 'Swap matrix computation time: %10.2f s'%matrix_time
        print 'Gibbs sampling time         : %10.2f s'%sampling_time
        print '------------------------------------------'
        print 'Total exchange time         : %10.2f s'%total_time



#     def _check_remote_resource(self, resource_url):
#         """
#         check if it's a remote resource. Basically see if 'ssh' is present
#         """
#         ssh_c = re.compile("(.+)\+ssh://(.*)")
#         m = re.match(ssh_c, resource_url)
#         if m:
#             self.remote_protocol = m.group(1)
#             self.remote_server = m.group(2)
#             print resource_url + " : yes" + " " + remote_protocol + " " + remote_server
#             return 1
#         else:
#             print resource_url + " : no"
#             return 0

#     def _setup_remote_workdir(self):
#         """
#         rsync local working directory with remote working directory
#         """
#         os.system("ssh %s mkdir -p %s" % (self.remote_server, self.keywords.get('REMOTE_WORKING_DIR')))
#         extfiles = " "
#         for efile in self.extfiles:
#             extfiles = extfiles + " " + efile
#         os.system("rsync -av %s %s/%s/" % (extfiles, self.remote_server, self.keywords.get('REMOTE_WORKING_DIR')))

#         dirs = ""
#         for k in range(self.nreplicas):
#             dirs = dirs + " r%d" % k
#         setup_script = """
# cd %s ; \
# for i in `seq 0 %d` ; do \
# mkdir -p r$i ; \
 




# """
      
    def _debug_collect_state_populations(self, replicas):
        """
        Calculate the empirically observed distribution of state permutations. 
        Permutations not observed will NOT be counted and will need to be
        added later for proper comparison to the exact distribution.
        """
        try:
            self.nperm
        except (NameError,AttributeError):
            self.nperm = {}
        curr_states = [self.status[i]['stateid_current'] for i in replicas]
        curr_perm = str(zip(replicas,curr_states))
        if self.nperm.has_key(curr_perm):
            self.nperm[curr_perm] += 1
        else:
            self.nperm[curr_perm] = 1

    def _debug_validate_state_populations(self, replicas, states, U):
        """
        Calculate the exact state permutation distribution and compare it to
        the observed distribution. The similarity of these distributions is
        measured via the Kullback-Liebler divergence.
        """
        empirical = sample_to_state_perm_distribution(self.nperm,replicas,
                                                      states)
        exact = state_perm_distribution(replicas,states,U)
        print '%8s %-9s %-9s %-s'%('','empirical','exact','state permutation')
        print '-'*80
        if len(empirical.keys()) > len(exact.keys()):
            perms = empirical.keys()
        else:
            perms = exact.keys()
        for k,perm in enumerate(perms):
            print '%8d %9.4f %9.4f %s'%(k+1,empirical[perm],exact[perm],perm)
        print '-'*80
        dkl = state_perm_divergence(empirical,exact)
        print 'Kullback-Liebler Divergence = %f'%dkl
        print '='*80
