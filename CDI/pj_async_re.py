# File Based Replica Exchange class
"""A module to prepare and run file-based asynchronous RE jobs
See documentation in doc/ directory.

Contributors: 

Emilio Gallicchio <emilio.gallicchio@rutgers.edu>
Brian Radak <radakb@biomaps.rutgers.edu>
Melissa Romanus <melissa.romanus@rutgers.edu>
Brian Radak <radakb@biomaps.rutgers.edu>

"""

import os, sys, time, re, pickle, random, math, copy
import pdb
import shutil, signal, glob
from configobj import ConfigObj
#pilotjob: Packages for Pilot API
from pilot import PilotComputeService, ComputeDataService, State

__version__ = '0.1.0'

class async_re_job:
    """
    Class to set up and run asynchronous file-based RE calculations
    """
    def __init__(self, command_file, options):
        self.command_file = command_file
        self.cus={}
        self.jobname = os.path.splitext( os.path.basename(command_file) )[0]

        self._parseInputFile()

        self._checkInput()

        self._printStatus()



    def _error(self, text):
        """ Print an error line to the log file """
        print text
        sys.stdout.flush()
    
    def _exit(self, text):
        """ Print an error and exit """
        self._error(text)
        print 'exiting...'
        sys.exit(1)

    def _parseInputFile(self):
        """
        Read keywords from control file. Requires the ConfigObj module
        """
        self.keywords = ConfigObj(self.command_file)
        #print self.keywords
        

    def _printStatus(self):
        """ Logs input parameters """
        print 'command_file=', self.command_file
        print 'jobname=', self.jobname
        for k, v in self.keywords.iteritems():
            print k, v

    def _openfile(self,name,mode):
        maxtries = 100
        tries = 0
        f = None
        while not f and tries <= maxtries:
            try:
                f = open(name,mode)
            except IOError:
                print "Warning: unable to access file %s, re-trying in 1 second ..." % name
                f = None
                tries += 1
                time.sleep(1)
        if tries > maxtries:
            self._exit("Too many failures accessing file %s: quitting." % name )
        return f

    def _checkInput(self):
        """ 
Checks that the required parameters are specified and parses
these and the other settings. 
"""
        #required options
        # the basename for the job
        self.basename = self.keywords.get('ENGINE_INPUT_BASENAME')
        if self.basename is None:
            self._exit("ENGINE_INPUT_BASENAME needs to be specified")
        # execution time in minutes
        self.walltime = float(self.keywords.get('WALL_TIME'))
        if self.walltime is None:
            self._exit("WALL_TIME needs to be specified")

        #pilotjob: Required variables for PilotJob
        if self.keywords.get('COORDINATION_URL') is None:
            self._exit("COORDINATION_URL needs to be specified")
        if self.keywords.get('RESOURCE_URL') is None:
            self._exit("RESOURCE_URL needs to be specified")
        if self.keywords.get('QUEUE') is None:
            self._exit("QUEUE needs to be specified")
        if self.keywords.get('BJ_WORKING_DIR') is None:
            self._exit("BJ_WORKING_DIR needs to be specified")
        if self.keywords.get('TOTAL_CORES') is None:
            self._exit("TOTAL_CORES needs to be specified")
        if self.keywords.get('SUBJOB_CORES') is None:
            self._exit("SUBJOB_CORES needs to be specified")

        #pilotjob: optional variables
        if self.keywords.get('PPN') is None: self.ppn = 1
        else: self.ppn = int(self.keywords.get('PPN'))

        if self.keywords.get('SPMD') is None: self.spmd = "single"
        else: self.spmd = self.keywords.get('SPMD')
        
        #number of replicas (may be determined by other means)
        if self.keywords.get('NREPLICAS') is not None:
            self.nreplicas = int(self.keywords.get('NREPLICAS'))
        else:
            self.nreplicas = None

        #initializes extfiles variable for 'setupJob'
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if self.extfiles is not None and self.extfiles != '':
            self.extfiles = self.extfiles.split(',')

    def setupJob(self):
        """
If RE_SETUP='yes' creates and populates subdirectories, one for each replica 
called r0, r1, ..., rN in the working directory. Otherwise reads saved state
from the ENGINE_BASENAME.stat file.

To populate each directory calls _buildInpFile(k) to prepare the MD engine
input file for replica k. Also creates soft links to the working directory 
for the accessory files specified in ENGINE_INPUT_EXTFILES.
"""
	#pilotjob: Initialize PilotJob at given COORDINATION_URL (CU)
        self.pj = PilotComputeService(self.keywords.get('COORDINATION_URL'))
	#pilotjob: Initialize PilotJob Data service (DU)
        self.cds=ComputeDataService()
	#pilotjob: Launch the PilotJob at the given COORDINATION_URL
        self.launch_pilotjob()

        if not (self.keywords.get('RE_SETUP') is None) and (self.keywords.get('RE_SETUP') == "yes"):
            # create replicas directories r1, r2, etc.
            for k in range(self.nreplicas):
                os.mkdir("r%d" % k)
            # create links for external files
            if self.extfiles != None:
                for k in range(self.nreplicas):
                    os.chdir("r%d" % k)
                    for file in self.extfiles:
                        if os.path.exists("../%s" % file) != True:
                            msg = "No such file: %s" % file
                            self._exit(msg)
                        os.symlink("../%s" % file, file)
                    os.chdir("..")
            # create status table
            self.status = []
            for k in range(self.nreplicas):
                st = {}
                st['stateid_current'] = k
                st['running_status'] = "W"
                st['cycle_current'] = 1
                self.status.append(st)
            # save status tables
            self._write_status()
            # create input files no. 1
            for k in range(self.nreplicas):
                self._buildInpFile(k)
            self.updateStatus()
        else:
            self._read_status()
            self.updateStatus(restart=True)

        self.print_status()
        #at this point all replicas should be in wait state
        for k in range(self.nreplicas):
            if self.status[k]['running_status'] != "W":
                self._exit("Internal error after restart. Not all jobs are in wait state")

    def scheduleJobs(self):
        # wait until bigjob enters executing
        while self.pilotcompute.get_state() != "Running":
            time.sleep(10)

        # Gets the wall clock time for a replica to complete a cycle
        # If unspecified it is estimated as 10% of job wall clock time
        # Note  
        replica_run_time = self.keywords.get('REPLICA_RUN_TIME')
        if self.keywords.get('REPLICA_RUN_TIME') is None:
            replica_run_time = self.walltime/10
            if replica_run_time < 1:
                replica_run_time = 1
        else:
            replica_run_time = int(self.keywords.get('REPLICA_RUN_TIME'))
        # double it to give time for current running processes 
        # and newly submitted processes to complete
        replica_run_time = 2*replica_run_time

        # Time in between cycles in seconds
        # If unspecified it is set as 30 secs
        if self.keywords.get('CYCLE_TIME')  is None:
            cycle_time = 30.0
        else:
            cycle_time = float(self.keywords.get('CYCLE_TIME'))

        start_time = time.time()
        end_time = start_time + 60*(self.walltime - replica_run_time) - cycle_time - 10
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
        pcd={"service_url":self.keywords.get('RESOURCE_URL'),
             "number_of_processes":self.keywords.get('TOTAL_CORES'),
             "working_directory": self.keywords.get('BJ_WORKING_DIR'),
             "queue":self.keywords.get('QUEUE'),
             "processes_per_node":self.ppn,
	     "project":self.keywords.get('PROJECT'),
             "walltime":int(self.keywords.get('WALL_TIME'))}
         
	#pilotjob: Create pilot job with above description
        self.pj.create_pilot(pilot_compute_description=pcd)
        self.cds.add_pilot_compute_service(self.pj)
        self.pilotcompute = self.pj.list_pilots()[0]
        
            
    def _write_status(self):
        """
Saves the current state of the RE job in the BASENAME.stat 
file using pickle
"""
        status_file = "%s.stat" % self.basename
        f = self._openfile(status_file, "w")
        pickle.dump(self.status, f)
        #pickle.dump(self.node_status, f)
        f.close()

    def _read_status(self):
        """
Loads the current state of the RE job from BASENAME.stat 
file using pickle
"""
        status_file = "%s.stat" % self.basename
        f = self._openfile(status_file, "r")
        self.status = pickle.load(f)
        f.close()

    def print_status(self):
        """
Writes to BASENAME_stat.txt a text version of the status of the RE job. 

It's fun to follow the progress in real time by doing:

watch cat BASENAME_stat.txt
"""
        logfile = "%s_stat.txt" % self.basename
        ofile = self._openfile(logfile,"w")
        log = "Replica  State  Status  Cycle \n"
        for k in range(self.nreplicas):
            log += "%6d   %5d  %5s  %5d \n" % (k, self.status[k]['stateid_current'], 
                self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()
        
    
    
#         
#    def _buildInpFile(self, replica):
#        """
#Generic function to prepare the input file for the MD engine for the specified
#replica for the current cycle and thermodynamic state. Calls specialized 
#functions depending on RE_TYPE and ENGINE.
#"""
#        if self.keywords.get('RE_TYPE') == 'BEDAM':
#            self._buildInpFile_BEDAM(self.basename, replica, self.status[replica]['stateid_current'], self.status[replica]['cycle_current'])


    def updateStatus(self,restart=False):
        """
Scans replica to update their state. 
"""
        for k in range(self.nreplicas):
            self._updateStatus_replica(k,restart)
        self._write_status()
        self._update_running_no()

    def _updateStatus_replica(self,replica,restart):
        """
Update the status of the specified replica. If it has completed a cycle the 
input file for the next cycle is prepared and the replica is placed in
the wait state.
"""
        this_cycle = self.status[replica]['cycle_current']


        if restart:

            if self.status[replica]['running_status'] == "R":
                if self._hasCompleted(replica,this_cycle):
                    self.status[replica]['cycle_current'] += 1
                else:
                    print "_updateStatus_replica(): Warning: restarting replica %d (cycle %d)" % (replica,this_cycle)
            self._buildInpFile(replica)
            self.status[replica]['running_status'] = "W"

        else:

            if self.status[replica]['running_status'] == "R":
                if self._isDone(replica,this_cycle):
                    self.status[replica]['running_status'] = "S"
                    if self._hasCompleted(replica,this_cycle):
                        self.status[replica]['cycle_current'] += 1
                    else:
                        print "_updateStatus_replica(): Warning: restarting replica %d (cycle %d)" % (replica,this_cycle)
                    self._buildInpFile(replica)
                    self.status[replica]['running_status'] = "W"
                            
    def _update_running_no(self):
        """
Updates the number of running replicas
"""
        self.running = 0
        self.waiting = 0
        for k in range(self.nreplicas):
            if self.status[k]['running_status'] == "R":
                self.running += 1
            if self.status[k]['running_status'] == "W":
                self.waiting += 1

    def _isDone(self,replica,cycle):
        """
Generic function to check if a replica completed a cycle. 
Calls in this case pilot-job version.
"""
        return self._isDone_PJ(replica,cycle)

    def _isDone_PJ(self,replica,cycle):
        """
Returns true if a replica has exited (done or failed)
"""
       	#pilotjob: Get status of the compute unit
	#pilotjob: Query the replica to see if it is in the done state
        state = self.cus[replica].get_state()
        details = self.cus[replica].get_details()

        if state == "Done" or state == "Failed" or state == "Canceled":
            if details.has_key("start_time"):
                if details.has_key("end_time"):
                    print "*********************************************************"
                    print "Replica: "+replica+"Start Time: "+float(details["start_time"])+"End Time: "+float(details["end_time"])
            if details.has_key("end_queue_time"):
                print "End Queue Time: "+float(details["end_queue_time"])+"\n"
            return True
        else:
            return False
            
      
    def _hasCompleted(self,replica,cycle):
        """
Attempts to check whether a replica has completed successfully from the bigjob
compute unit. This is not expected to work during a restart when compute units
are not available. In the latter case success is assumed. MD engine modules
are recommended to override this default routine with one that implements 
a better test of success such as the existence of a restart file or similar.
""" 
        try:
            state = self.cus[replica].get_state()
        except:
            print "_hasCompleted(): Warning: unable to query replica state. Assuming success ..."
            return True

        if state == "Done":
            return True
        else:
            return False


    def _njobs_to_run(self):
        # size of subjob buffer as a percentage of job slots (TOTAL_CORES/SUBJOB_CORES)
        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is None:
            subjobs_buffer_size = 0.5
        else:
            subjobs_buffer_size = float(subjobs_buffer_size)
        # find out how many replicas are waiting/(running/submitted)
        self._update_running_no()
        # launch new replicas if the number of submitted/running subjobs is less than
        # the number of available slots (total_cores/subjob_cores) + 50%
        available_slots = int(int(self.keywords.get('TOTAL_CORES'))/int(self.keywords.get('SUBJOB_CORES')))
        max_njobs_submitted = int((1.+subjobs_buffer_size)*available_slots)
        nlaunch = self.waiting - max(2,self.nreplicas - max_njobs_submitted)
        nlaunch = max(0,nlaunch)
        if self.keywords.get('VERBOSE') == "yes":
            print "available_slots: %d" % available_slots
            print "max_njobs_submitted: %d" % max_njobs_submitted
            print "running/submitted subjobs: %d" % self.running
            print "waiting replicas: %d" % self.waiting
            print "replicas to launch: %d" % nlaunch
        return nlaunch

    def launchJobs(self):
        """
Scans the replicas in wait state and randomly launches some of them
if CPU's are available.
""" 

        jobs_to_launch = self._njobs_to_run()
        if jobs_to_launch > 0:
            wait = []
            for k in range(self.nreplicas):
                if self.status[k]['running_status'] == "W":
                    wait.append(k)
            random.shuffle(wait)
            n = min(jobs_to_launch,len(wait))
            for k in wait[0:n]:
                if self.keywords.get('VERBOSE') == "yes":
                    print "Launching replica %d cycle %d" % (k,self.status[k]['cycle_current'])
                self.cus[k] = self._launchReplica(k,self.status[k]['cycle_current'])
                self.status[k]['running_status'] = "R"


    def doExchanges(self):
        """
Randomly selects a pair of replicas in wait state for exchange of 
thermodynamic parameters. 
"""
        print "Entering doExchanges Method: "+time.time()
        # find out which replicas are waiting
        self._update_running_no()
        if self.waiting > 1:
            replicas_waiting = []
            for k in range(self.nreplicas):
                if self.status[k]['running_status'] == "W" and self.status[k]['cycle_current'] > 1:
                    replicas_waiting.append(k)
            random.shuffle(replicas_waiting)
#           perform exchanges in pairs
            for k in range(0,len(replicas_waiting)-1,2):
                repl_a = replicas_waiting[k]
                repl_b = replicas_waiting[k+1]
                """
1. Places the two replicas in "E" (exchanging state)
2. rewinds the cycle to the previous completed cycle
3. performs the exchange by calling exchange routine from specialized 
   modules.
4. Creates new input file for the next cycle
5. Places replicas back into "W" (wait) state 
                """
                self.status[repl_a]['running_status'] = "E"
                self.status[repl_b]['running_status'] = "E"
                self.status[repl_a]['cycle_current'] -= 1
                self.status[repl_b]['cycle_current'] -= 1
                self._doExchange_pair(repl_a,repl_b)
                self.status[repl_a]['cycle_current'] += 1
                self.status[repl_b]['cycle_current'] += 1
                self._buildInpFile(repl_a)
                self._buildInpFile(repl_b)
                self.status[repl_a]['running_status'] = "W"
                self.status[repl_b]['running_status'] = "W"

        print "Exiting doExchanges Method: "+time.time()

