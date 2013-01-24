import os, re, random, math
from pj_async_re import async_re_job
import numpy as np

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        AMBERHOME = os.getenv('AMBERHOME')
        if AMBERHOME is None:
            self._exit('Cannot find AMBERHOME.')
            
        # Check which AMBER MD engine to use, default to sander
        engine = self.keywords.get('ENGINE').upper()
        sander_flags = [ 'AMBER', 'SANDER', 'AMBER-SANDER' ]
        pmemd_flags = [ 'PMEMD', 'AMBER-PMEMD' ]
        engine_name = ''
        if engine in sander_flags:  engine_name = 'sander'
        elif engine in pmemd_flags: engine_name = 'pmemd'
        else:                       self._exit('ENGINE is not AMBER')
        if self.spmd == 'mpi': engine_name += '.MPI'
        # else just assume that a serial executable is desired
        # TODO?: Cuda

        # Check that this executable exists, etc.
        self.exe = os.path.join(AMBERHOME,'bin',engine_name)
        if not os.path.exists(self.exe) or not os.access(self.exe,os.X_OK):
            self._exit('Could not find an executable called %s, expected it to'
                       ' be at %s'%(engine_name,self.exe))

        print 'AMBERHOME:',AMBERHOME
        print 'exe:',self.exe

    def _launchReplica(self,replica,cycle):
        """Launches Amber sub-job using pilot-job
        """
        input_file = "%s_%d.inp" % (self.basename, cycle)
        out_file = "%s_%d.out" % (self.basename, cycle)
        prm_file = "%s.parm7" % self.basename
        crd_file = "%s_%d.rst7" % (self.basename, cycle-1)
        rst_file = "%s_%d.rst7" % (self.basename, cycle)
        xyz_file = "%s_%d.nc" % (self.basename, cycle)
        info_file = "%s_%d.info" % (self.basename, cycle)
        log_file = "%s_%d.log" % (self.basename, cycle)
        err_file = "%s_%d.err" % (self.basename, cycle)

        arguments = ["-O",
                     "-i", input_file,
                     "-o", out_file, 
                     "-p", prm_file, 
                     "-c", crd_file, 
                     "-r", rst_file, 
                     "-x", xyz_file, 
                     "-inf", info_file]

        #pilotjob: Compute Unit (i.e. Job) description
        compute_unit_description = {
            "executable": self.exe,
            "environment": [],
            "arguments": arguments,
            "output": log_file,
            "error": err_file,   
            "working_directory":os.getcwd()+"/r"+str(replica),
            "number_of_processes": int(self.keywords.get('SUBJOB_CORES')),
            "spmd_variation": self.spmd,
            }

        if self.keywords.get('VERBOSE') == "yes":
            print ( "Launching %s in directory %s (cycle %d)" % 
                    (self.exe.split('/')[-1], os.getcwd()+"/r"+str(replica), 
                     cycle) )

#        compute_unit=self.cds.submit_compute_unit(compute_unit_description)
        compute_unit=self.pilotcompute.submit_compute_unit(compute_unit_description)
        return compute_unit
        
    def _isDone(self,replica,cycle):
        """
        Returns true if an Amber replica has completed a cycle. Basically 
        checks if the restart file exists.
        This overrides the generic isDone using pilot-job, which cannot check if
        a replica is done after a restart.
        """
        rstfile = "r%d/%s_%d.rst7" % (replica, self.basename,cycle)
        if os.path.exists(rstfile):
            return True
        else:
            return False

###########################################################################
#
# Work in Progress: Gibbs sampling style exchanges (see impact_async_re.py)
#
###########################################################################
# gives random choice from a set with weight probabilities
    def _weighted_choice_sub(self,weights):
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
                
    def _gibbs_re_j(self,i,nstates):
        # produces a replica "j" to exchange with the given replica "i"
        ee = [ self._reduced_energy(j,j) for j in range(nstates) ]

        ps = np.zeros(nstates)
        for j in range(nstates):
            # energy after (i,j) exchange
            eij = self._reduced_energy(i,j) + self._reduced_energy(j,i)
            ps[j] = -(eij - ee[i] - ee[j])
        ps = np.exp(ps)
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

        # backtrack cycle
        for k in replicas_waiting:
            self.status[k]['cycle_current'] -= 1
            self.status[k]['running_status'] = "E"

        for reps in range(1):
            for i in range(len(replicas_waiting)):
                j = self._gibbs_re_j(i,len(replicas_waiting))
                if i != j:
                    #swap state id's
                    ri = replicas_waiting[i]
                    rj = replicas_waiting[j]
                    sid_i = self.status[ri]['stateid_current'] 
                    sid_j = self.status[rj]['stateid_current']
                    self.status[ri]['stateid_current'] = sid_j
                    self.status[rj]['stateid_current'] = sid_i

        # write input files
        for k in replicas_waiting:
            # Creates new input file for the next cycle
            # Places replica back into "W" (wait) state 
            self.status[k]['cycle_current'] += 1
            self._buildInpFile(k)
            self.status[k]['running_status'] = "W"
