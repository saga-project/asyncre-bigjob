import os, re, random, math
from numpy import *
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

#
# Experimental Gibbs sampling for RE
#


# gives random choice from a set with weight probabilities
    def _weighted_choice_sub(self,weights):
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
                
    def _gibbs_re_j(self,i,par,pot):
        # produces a replica "j" to exchange with the given replica "i"
#        re_etot = 0
        n = len(par)
        ee = []
        for j in range(n):
            ej = self._reduced_energy(par[j],pot[j])
            ee.append(ej)
#            re_etot += ej
        ei = ee[i]
        ps = zeros(n)
        for j in range(n):
            # energy after (i,j) exchange
            eij = self._reduced_energy(par[i],pot[j]) + self._reduced_energy(par[j],pot[i])
            # note that total energy is not included, since it is the same for all of the
            # permutation states.
            ps[j] = -(eij - ei - ee[j])
#        ps.append(math.exp(-(eij - ei - ee[j])))
        ps = exp(ps)
        return self._weighted_choice_sub(ps)

    def doExchanges(self):
        """
Perform n rounds of exchanges among waiting replicas using Gibbs sampling.
"""
        # find out which replicas are waiting
        self._update_running_no()
        if self.waiting > 1:
            replicas_waiting = []
            for k in range(self.nreplicas):
                if self.status[k]['running_status'] == "W" and self.status[k]['cycle_current'] > 1:
                    replicas_waiting.append(k)

#        replicas_waiting = [replicas_waiting[4], replicas_waiting[8], replicas_waiting[12], replicas_waiting[16] ]

        # backtrack cycle
        for k in replicas_waiting:
            self.status[k]['cycle_current'] -= 1
            self.status[k]['running_status'] = "E"
        #collect replica parameters and potentials
        par = []
        pot = []
        for k in replicas_waiting:
            v = self._getPot(k,self.status[k]['cycle_current'])
            l = self._getPar(k)
            par.append(l)
            pot.append(v)
        # perform an exchange for each of the n replicas
        print pot
        print par
#        for reps in range(len(replicas_waiting)):
        for reps in range(1):


#        npermt = {}
#        permt = {}
#        for reps in range(1000):

            for i in range(len(replicas_waiting)):
                j = self._gibbs_re_j(i,par,pot)
                if i != j:
                    #swap state id's
                    ri = replicas_waiting[i]
                    rj = replicas_waiting[j]
                    sid_i = self.status[ri]['stateid_current'] 
                    sid_j = self.status[rj]['stateid_current']
                    self.status[ri]['stateid_current'] = sid_j
                    self.status[rj]['stateid_current'] = sid_i
                    #swap parameters
                    tmp = par[i]
                    par[i] = par[j]
                    par[j] = tmp

#            print par
#            key = str(par)
#            if npermt.has_key(key):
#                npermt[key] += 1
#            else:
#                npermt[key] = 1
#                permt[key] = copy.copy(par)

#        ss = 0
#        for k in npermt.keys():
#            ss += npermt[k]
#        ps = []
#        sumps = 0
#        for k in npermt.keys():
#            b = permt[k]
#            e = 0
#            for i in range(len(replicas_waiting)):
#                e += self._reduced_energy(b[i],pot[i])
#            p = math.exp(-e)
#            ps.append(p)
#            sumps += p
#    
#        print pot
#        i = 0
#        for k in npermt.keys():
#            print k, npermt[k]/float(ss), ps[i]/sumps
#            i += 1
        
        print par

        # write input files
        for k in replicas_waiting:
            # Creates new input file for the next cycle
            # Places replica back into "W" (wait) state 
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = "W"


