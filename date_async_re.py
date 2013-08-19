import os, sys, time
from pj_async_re import async_re_job

class pj_date_job(async_re_job):

    def _launchReplica(self,replica,cycle):
        """
Issues a command to launch /bin/date using PJ
"""
	#pilotjob: Compute Unit (i.e. Job) description
        compute_unit_description = {
            "executable": "/bin/date",
            "arguments": [""],
            "total_cpu_count": int(self.keywords.get('SUBJOB_CORES')),            
            "output": "sj-stdout-"+str(replica)+"-"+str(cycle)+".txt",
            "error": "sj-stderr-"+str(replica)+"-"+str(cycle)+".txt",   
            "working_directory":os.getcwd()+"/r"+str(replica),
            "spmd_variation":self.keywords.get('SPMD')
        }  
        if self.keywords.get('VERBOSE') == "yes":
            print "Launching %s in directory %s cycle %d" % ("/bin/date",os.getcwd()+"/r"+str(replica),cycle)
        compute_unit=self.cds.submit_compute_unit(compute_unit_description)

        #self.cus[replica]=compute_unit
        return compute_unit


class date_async_re_job(pj_date_job,async_re_job):
                            
    def _checkInput(self):
        async_re_job._checkInput(self)
        #make sure DATE is wanted
        if self.keywords.get('RE_TYPE') != 'DATE':
            self._exit("RE_TYPE is not DATE")
        #number of replicas
        if self.keywords.get('NREPLICAS') is None:
            self._exit("NREPLICAS needs to be specified")
        self.nreplicas = int(self.keywords.get('NREPLICAS'))

    def _buildInpFile(self, replica):
        pass

    def _doExchange_pair(self,repl_a,repl_b):
        pass

    def _computeSwapMatrix(self, replicas, states):
        U = [[ 0. for j in range(self.nreplicas)] 
             for i in range(self.nreplicas)]
        return U

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"
    
    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)
    
    commandFile = sys.argv[1]

    print ""
    print "===================================="
    print "DATE Asynchronous Replica Exchange "
    print "===================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    rx = date_async_re_job(commandFile, options=None)

    rx.setupJob()

    while rx.pilotcompute.get_state() != "Running":
        time.sleep(2)
    
# Gets the wall clock time for a replica to complete a cycle
# If unspecified it is estimated as 10% of job wall clock time  
    replica_run_time = rx.keywords.get('REPLICA_RUN_TIME')
    if replica_run_time is None:
        replica_run_time = rx.walltime/10

# Time in between cycles in seconds
# If unspecified it is set as 30 secs
    if rx.keywords.get('CYCLE_TIME')  is None:
        cycle_time = 30.0
    else:
        cycle_time = float(rx.keywords.get('CYCLE_TIME'))

    start_time = time.time()
    while time.time() < start_time + 60*(rx.walltime - replica_run_time):
        time.sleep(1)

        rx.updateStatus()
        rx.print_status()
        rx.launchJobs()
        
        time.sleep(cycle_time)

        rx.updateStatus()
        rx.print_status()        
        rx.doExchanges()
        
    rx.updateStatus()
    rx.print_status()

    rx.waitJob()
    rx.cleanJob()
