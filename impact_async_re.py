import os, re, random, math
from pj_async_re import async_re_job

class pj_impact_job(async_re_job):

    def _launchReplica(self,replica,cycle):
         """
Launches Impact sub-job using pilot-job
"""
         num_threads = os.getenv('OMP_NUM_THREADS')
         if num_threads == None:
             num_threads = 1
         else:
             num_threads = int(num_threads)

         input_file = "%s_%d.inp" % (self.basename, cycle)
         log_file = "%s_%d.log" % (self.basename, cycle)
         err_file = "%s_%d.err" % (self.basename, cycle)

         schrod_env = None

          # Parallelism
#           'number_of_processes': <Total number of processes to start>,
#           'processes_per_host':  <Nr of processes per host>,
#           'threads_per_process': <Nr of threads to start per process>,
#           'total_core_count':    <Total number of cores requested>,
#           'spmd_variation':      <Type and startup mechanism>,

	 #pilotjob: Compute Unit (i.e. Job) description
#         compute_unit_description = {
#            "executable": os.getcwd()+"/runimpact",
#            "environment": schrod_env,
#            "arguments": [input_file],
#            "number_of_processes": 1,
#            "threads_per_process": int(self.keywords.get('SUBJOB_CORES')),
#            "total_core_count": int(self.keywords.get('SUBJOB_CORES')),
#            "output": log_file,
#            "error": err_file,   
#            "working_directory":os.getcwd()+"/r"+str(replica),
#            "spmd_variation":self.keywords.get('SPMD')
#         }  


#            "environment": schrod_env,

	 #pilotjob: Compute Unit (i.e. Job) description
         compute_unit_description = {
            "executable": os.getcwd()+"/runimpact",
            "arguments": [input_file],
            "total_cpu_count": int(self.keywords.get('SUBJOB_CORES')),
            "output": log_file,
            "error": err_file,   
            "working_directory":os.getcwd()+"/r"+str(replica),
            "spmd_variation":self.keywords.get('SPMD')
         }  

         if self.keywords.get('VERBOSE') == "yes":
            print "Launching %s %s in directory %s cycle %d" % (os.getcwd()+"/runimpact",input_file,os.getcwd()+"/r"+str(replica),cycle)

#         compute_unit=self.cds.submit_compute_unit(compute_unit_description)
         compute_unit=self.pilotcompute.submit_compute_unit(compute_unit_description)
         return compute_unit

    def _getImpactData(self, file):
        """
Reads all of the Impact simulation data values temperature, energies, etc.
at each time step and puts into a big table
"""
        if not os.path.exists(file):
            msg = 'File does not exist: %s' % file
            self._exit(msg)
        step_line = re.compile("^ Step number:")
        number_line = re.compile("(\s+-*\d\.\d+E[\+-]\d+\s*)+")
        nsamples = 0
        data = []
        f = self._openfile(file ,"r")
        line = f.readline()
        while line:
            # fast forward until we get to the line: 
            # "Step number: ... "
            while line and not re.match(step_line, line): 
                line = f.readline()
            # read the step number
            if re.match(step_line, line):
                words = line.split()
                step = words[2]
                #now read up to 3 lines of numbers
                datablock = [int(step)]
                ln = 0
                while ln < 3:
                    line = f.readline()
                    if not line:
                        msg = "Unexpected end of file"
                        self._exit(msg)
                    if re.match(number_line, line):
                        for word in line.split():
                            datablock.append(float(word))
                        ln += 1
                data.append(datablock)
            line = f.readline()
        f.close()
        return data
        
    def _hasCompleted(self,replica,cycle):
        """
Returns true if an IMPACT replica has successfully completed a cycle.
"""
        try:
            #check existence of rst file
            rstfile = "r%d/%s_%d.rst" % (replica, self.basename,cycle)
            if not os.path.exists(rstfile):
                print "Warning: can not find file %s." % rstfile 
                return False
            #check that rst file is of the correct size
            if cycle > 1:
                rstfile_p = "r%d/%s_%d.rst" % (replica, self.basename,cycle-1)
                rstsize = os.path.getsize(rstfile)
                rstsize_p = os.path.getsize(rstfile_p)
                if not rstsize == rstsize_p:
                    print "Warning: files %s and %s have different size" % (rstfile,rstfile_p)
                    return False
            #check that we can read data from .out
            output_file = "r%s/%s_%d.out" % (replica,self.basename,cycle)
            datai = self._getImpactData(output_file)
            nf = len(datai[0])
            nr = len(datai)
        except:
            rstfile = "r%d/%s_%d.rst" % (replica, self.basename,cycle)
            rstfile_p = "r%d/%s_%d.rst" % (replica, self.basename,cycle-1)
            output_file = "r%s/%s_%d.out" % (replica,self.basename,cycle)
            print "Warning: unable to access some of these files: %s %s %s." % (rstfile,rstfile_p,output_file)
            return False
        return True

    #compute matrix of dimension-less energies: each column is a replica 
    #and each row is a state
    #so U[i][j] is the energy of replica j in state i. 
    #
    #Note that the matrix is sized to include all of the replicas and states 
    #but the energies of replicas not 
    #in waiting state, or those of waiting replicas for states not belonging to 
    #waiting replicas list are undefined.
    def _computeSwapMatrix(self, replicas, states):
        # U will be sparse matrix, but is convenient bc the indices of the
        # rows and columns will always be the same.
        U = [[ 0. for j in range(self.nreplicas)] 
             for i in range(self.nreplicas)]

        n = len(replicas)

        #collect replica parameters and potentials
        par = []
        pot = []
        for k in replicas:
            v = self._getPot(k,self.status[k]['cycle_current'])
            l = self._getPar(k)
            par.append(l)
            pot.append(v)
        print pot
        print par   

        for i in range(n):
            repl_i = replicas[i]
            for j in range(n):
                sid_j = states[j]
                # energy of replica i in state j
                U[sid_j][repl_i] = self._reduced_energy(par[j],pot[i])
        return U


