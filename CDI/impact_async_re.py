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

	 #pilotjob: Compute Unit (i.e. Job) description
         compute_unit_description = {
            "executable": os.getcwd()+"/runimpact",
            "environment": schrod_env,
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
        f = open(file ,"r")
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
        
    def _isDone(self,replica,cycle):
        """
Returns true if an IMPACT replica has completed a cycle. Basically checks
if the restart file exists.
This overrides the generic isDone using pilot-job, which cannot check if a
replica is done after a restart.
"""
        rstfile = "r%d/%s_%d.rst" % (replica, self.basename,cycle)
        if os.path.exists(rstfile):
            return True
        else:
            return False
