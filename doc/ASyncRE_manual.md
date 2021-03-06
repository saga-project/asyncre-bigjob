Asynchronous Replica Exchange (ASyncRE)
=======================================
v. 0.2.1


User Manual
-----------

by 

Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>

and

Brian Radak <radakb@biomaps.rutgers.edu>


Package Developers
------------------

Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>

Brian Radak <radakb@biomaps.rutgers.edu>

Melissa Romanus <melissa.romanus@rutgers.edu>

Introduction
------------

ASyncRE is a Python package to perform file-based asynchronous parallel replica exchange molecular simulations. 

The current implementation is aimed at computer clusters managed by a queuing system and supported by a shared filesystem. The [BigJob distributed computing infrastructure](https://github.com/saga-project/BigJob/wiki) is used for job launching and monitoring. Instructions on how to install BigJob on a cluster are available at the following [link](http://saga-project.github.io/BigJob/sphinxdoc/install/install.html). While primarily directed at NSF XSEDE clusters, BigJob supports most cluster configurations.

The ASyncRE package includes a core module which performs common tasks such as job staging through BigJob and exchanging of parameters among replicas. Support for arbitrary MD engines and RE schemes are introduced through user-provided modules. Currently, MD engine modules are available for the AMBER and IMPACT MD programs. A similar modular mechanism provides support for arbitrary RE schemes (temperature, Hamiltonian, etc.), including arbitrary multidimensional combinations of these (such as 2D RE temperature/Hamiltonian). The software is currently distributed with modules for multidimensional RE umbrella sampling with AMBER, and BEDAM lambda-RE alchemical binding free energy calculations with the Impact MD engine.

The algorithm implemented by ASyncRE can be summarized as follows:

1. Job files and executables for the replicas are set up as appropriate for the MD engine/RE scheme combination as specified by a user-provided module. Typically this is accomplished by parsing a set of template input files according to the thermodynamic and potential energy settings which distinguish one replica from another. Each replica lives in a separate sub-directory of the working directory. These replica sub-directories are named r0, r1, ..., r<M-1>, where M is the number of replicas.

2. Periodically, a subset of the replicas are submitted to BigJob for execution and enter a "R" (running) state. When a replica completes a run (of length specified in the MD engine input file), referred to in what follows as a "cycle", it enters a "W" (waiting) state which makes it eligible for exchange of thermodynamic and other parameters with other replicas and for the initiation of a new run cycle.

3. Periodically, exchanges of thermodynamic parameters are performed between the replicas in the waiting state based on the appropriate replica exchange rules (as specified in user modules) based on their thermodynamic states (temperature for example), or potential energy settings, together with the necessary structural and energetic information obtained from the MD engine output files.

Internally a dictionary named `status` is used to keep track of the status of the replicas (the current cycle, thermodynamic state, running status, etc.). The main items are:

<dl>
 <dt>status[repl]['stateid_current']: </dt> 
     <dd>The id of the current thermodynamic state held by replica `repl`. The state id is a unique integer id from 0 to N-1 assigned to each thermodynamic state. During exchanges, replicas swap state id's. The mapping between thermodynamic parameters and id's is user-defined. See for example the `_buildBEDAMStates()` routine in the `bedamtempt_async_re.py` module.</dd>

 <dt>status[repl]['running_status']: </dt>
     <dd>The running status of replica `repl`. It is either "R" for running (either in the buffer area waiting to execute or actually consuming CPU cycles), or "W" for waiting (to run). When in "W" state the replica undergoes parameter exchanges with other replicas in the "W" state. When a running replica finishes a cycle (whether successfully or unsuccessfully) it transitions from "R" to "W". Vice versa when a replica is submitted to BigJob it enters the "R" state.</dd>

 <dt>status[repl]['cycle_current']: </dt>
     <dd>The current cycle of replica `repl`. A cycle of n means that the replica has completed n-1 runs and it is either running or waiting to execute the nth run. </dd>
</dl>

The `status` data structure is check-pointed periodically to a pickle file called `<basename>.stat` in the working directory. When restarting, the `status` data structure is restored from this file. 

Installation
------------

See [README.md](https://github.com/saga-project/asyncre-bigjob/blob/master/README.md) in top-level directory.

Execution
---------

A typical command to initiate an ASyncRE simulations is as follows:

    ssh cluster_head_node
    source ~/.bigjob/bin/activate
    cd working_directory
    python amber_us.py control_file.cntl > LOG 2>& &

The second command above activates the virtualenv python environment (see BigJob documentation). The 'amber_us.py' is a short user-provided python script that loads the appropriate modules and launches the asynchronous RE simulation. For AMBER/umbrella sampling it is something like:

    import sys
    from amberus_async_re import amberus_async_re_job
    # Parse control file
    usage = "%prog <ConfigFile>"
    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)    
    commandFile = sys.argv[1]
    # initializes asynchronus RE simulation
    rx = amberus_async_re_job(commandFile, options=None)
    rx.setupJob()
    # starts submitting replicas to BigJob
    rx.scheduleJobs()

For RE modalities not natively supported by the current ASyncRE package (see "Writing extension modules" below) the execution script is usually preceded by a custom class definition. For example (again, see below for the details):

    import sys, math, os, ...
    from amber_async_re import pj_amber_job
     
    class mycrazyREscheme_async_re_job(pj_amber_job):
        def _checkInput(self):
           ...
        def _buildInpFile(self, replica):
           ...
        def ...
           ...
     
    if __name__ == '__main__':
        # Parse control file
        usage = "%prog <ConfigFile>"
        if len(sys.argv) != 2:
          print "Please specify ONE input file"
          sys.exit(1)    
        commandFile = sys.argv[1]
        # initializes asynchronus RE simulation
        rx = mycrazyREscheme_async_re_job(commandFile, options=None)
        rx.setupJob()
        # starts submitting replicas to BigJob
        rx.scheduleJobs()

The `control_file.cntl` above refers to a file contains keyword=value pairs, such as:

    #-Main setting------------------------
    ENGINE = 'AMBER'
    RE_TYPE = 'AMBERUS'
    RE_SETUP = 'yes'
    VERBOSE="yes"
    ENGINE_INPUT_BASENAME = 'DMP_US'
    ENGINE_INPUT_EXTFILES = 'DMP_US.parm7,DMP_US_0.rst7'
    #-RE/simulation settings---------------
    FORCE_CONSTANTS = '5.0,5.0:5.0,5.0:5.0,5.0:5.0,5.0:5.0,5.0:5.0,5.0'
    BIAS_POSITIONS = '275.,275.:275.,280.:280.,275.:280.,280.:280.,285.:285.,280'
    #-BigJob settings----------------------
    WALL_TIME=200
    COORDINATION_URL="redis://ILikeBigJob_wITH-REdIS@gw68.quarry.iu.teragrid.org:6379"
    RESOURCE_URL="pbs://localhost"
    QUEUE="batch"
    BJ_WORKING_DIR='/home/radakb/amber_us_2d/agent'
    TOTAL_CORES=16
    SUBJOB_CORES=8
    #---------------------------------------

The command above causes, among other things, the submission of a job to the local queuing system named "bliss_job" (a fixed name assigned by BigJob). Execution terminates after a specified amount of wall-clock time. The internal state of the simulation is check-pointed periodically and at the end of execution so that it can be restarted at a later time (see below). Failed runs are automatically detected and the relevant replicas are reset and restarted. 

ASyncRE has proven to be quite robust. Jobs with ~100's of replicas and 1000's of CPU cores running continuously for 1-2 days are routinely conducted by our groups on XSEDE.

Control settings
----------------

As shown above the control file is a simple keyword/value pair listing. The syntax follows the ConfigObj python package. See the [ConfigObj documentation](http://www.voidspace.org.uk/python/configobj.html) for details. The keywords recognized by ASyncRE are as follows:

**Main settings:**

<dl>
<dt>RE_TYPE</dt>
<dd>Identifies the replica exchange scheme. Required. "AMBERUS" and "BEDAMTEMPT" are currently the two natively recognized types (corresponding to AMBER Umbrella Sampling, and BEDAM+Temperature respectively). However feasible values vary depending on the available extension application modules in your installation.</dd>

<dt>ENGINE</dt>
<dd>Sets the MD engine. Required. "AMBER" and "IMPACT" are currently the two natively recognized values. Feasible values vary depending on the available extension application modules in your installation.</dd>

<dt>ENGINE_INPUT_BASENAME</dt>
<dd>Basename of the job. Required. Used, depending on the application, to locate/create input files and associated files, and to write the check-pointing files "ENGINE_INPUT_BASENAME.stat" and "ENGINE_INPUT_BASENAME_stat.txt". The latter lists the current status of the replicas (cycle number, state, running/waiting, etc.).</dd>

<dt>RE_SETUP</dt>
<dd>Whether to setup a new RE simulation (create replica directories, etc.). 'no' is used to restart a previously interrupted RE job. Defaults to 'no'. </dd>

<dt>ENGINE_INPUT_EXTFILES</dt>
<dd>List of structure files etc. that are copied from working directory to the replicas directories to start each replica. Default to the null value.</dd>

<dt>VERBOSE</dt>
<dd>If set to 'yes' prints detailed information on the progress of the simulation, exchanges, etc. Defaults to 'no'.</dd>
</dl>

**BigJob-related settings:**

<dl>
<dt>TOTAL_CORES</dt>
<dd>The number of CPU cores requested from the queuing system. On most cluster configurations the corresponding request of compute nodes is determined automatically. See BigJob documentation. TOTAL_CORES should be smaller than the number of replicas otherwise few replicas will be found in the waiting state at any one time and as a result exchanges will occur with insufficient frequency. A good value for TOTAL_CORES is such so as to accommodate roughly half of the replicas. Defaults to "1".</dd>

<dt>PPN</dt>
<dd>Processes per node. Required by BigJob on some architectures. Defaults to "1".</dd>

<dt>SUBJOB_CORES</dt>
<dd>The number of CPU cores utilized by each replica. Set as needed based on parallelism. Defaults to "1".</dd>

<dt>SPMD</dt>
Type of replica parallel execution. Could be either "single" or "mpi". See BigJob documentation. Defaults to "single".</dd>

<dt>SUBJOBS_BUFFER_SIZE</dt>
<dd>The size of the job buffer area expressed as a Fraction of TOTAL_CORES. When a replica completes execution BigJob immediately launches a new one taken from this buffer instead of waiting for a replica to be submitted. Defaults to 0.5.</dd>

<dt>WALL_TIME</dt>
<dd>Requested execution time in minutes. Time during which ASyncRE is waiting for the queued BigJob to begin execution is not counted towards this limit. This value is also passed to the queuing system as a job attribute. ASyncRE stops submitting replicas shortly before WALL_TIME is exceeded (see REPLICA_RUN_TIME below) to give time replicas to complete execution. No default, required setting.</dd>

<dt>REPLICA_RUN_TIME</dt>
<dd>Estimated wall clock time in minutes for a replica to complete a cycle. Used to determine when to stop submitting jobs to BigJob.  See WALL_TIME above. If unspecified it is estimated as 10% of job wall clock time. </dd>

<dt>CYCLE_TIME</dt>
<dd>Period in seconds between exchanges. This also sets the frequency with which the status of running replicas is updated. Defaults to 30 seconds. Note that setting it to a too small value can easily overwhelm the cluster head node and the filesystem, especially when dealing with many replicas and file/reading writing and computations related to exchanges are expensive.</dd>

<dt>QUEUE</dt>
<dd>The name of the queue where to submit the BigJob. Consult the cluster documentation for the appropriate queue. When not set the default queue may be selected.</dd>

<dt>PROJECT</dt>
<dd>Accounting string for the computing resource. Something like "5674209". Defaults to the null value.</dd>

<dt>BJ_WORKING_DIR</dt>
<dd>The directory where BigJob stores log files etc. Required setting.</dd>

<dt>COORDINATION_URL</dt>
<dd>The address of a suitable redis server. See the BigJob documentation. Required setting.</dd>

<dt>RESOURCE_URL</dt>
<dd>The address of the computing resource where to submit the BigJob. See BigJob documentation. Required setting.</dd>
</dl>

**Application and MD-engine specific settings:**

These are parsed and interpreted by application extension modules. See below for native module or documentation provided with the extension modules.


Writing extension modules
-------------------------

The ASyncRE package is structured around a core python class named `async_re_job`. Extension modules are user-provided derivative classes of `async_re_job` which define routines not available in the parent class or override generic routines of the parent class. For example the definition of the `bedamtemp_async_re_job` class which implements BEDAM RE applications starts with:

    class bedamtempt_async_re_job(async_re_job):
       
      def _checkInput(self):
          async_re_job._checkInput(self)
          #make sure BEDAM + TEMPERATURE is wanted
          if self.keywords.get('RE_TYPE') != 'BEDAMTEMPT':
              self._exit("RE_TYPE is not BEDAMTEMPT")
          #BEDAM runs with IMPACT
          if self.keywords.get('ENGINE') != 'IMPACT':
              self._exit("ENGINE is not IMPACT")
          ....
      ....

The `bedamtempt_async_re_job` class inherits all of the methods of the `async_re_job` class. However the new class can define new methods or override existing ones. For example the `_checkInput()` method above, used to parse the control file, overrides the generic one present in the parent `async_re_job` class. The new `_checkInput()` method calls the parent one to obtain from the control file the ASyncRE core settings described above. This is followed by commands to parse the settings specific to the application in question (starting with RE_TYPE and ENGINE in this case). 

This basic class inheritance strategy is used throughout to write extension modules. An extension class is required to provide at a minimum the following methods called by the core class:

<dl>
<dt>_checkInput(self):</dt>
<dd> Parses application and MD engine settings from the control file. For example here's how to read a list of temperatures:</dd>
</dl>

      def _checkInput(self):
          async_re_job._checkInput(self)
          ...
          #list of temperatures
          if self.keywords.get('TEMPERATURES') is None:
              self._exit("TEMPERATURES needs to be specified")
          self.temperatures = self.keywords.get('TEMPERATURES').split(',')
          ....
<dl>
<dt>_buildInpFile(self, repl):</dt> 
<dd>Creates the input files for replica replica 'repl' to prepare it for execution at the current cycle. For example the current temperature setting for the replica may be used to replace a '@temperature@' placeholder from a template input file:</dd>
</dl>

    def _buildInpFile(self, replica):
        # gets job basename
        basename = self.basename
        # thermodynamic state id
        stateid = self.status[replica]['stateid_current']
        # current cycle for this replica
        cycle = self.status[replica]['cycle_current']
        # path of template input file (in working directory)
        template = "%s.inp" % basename
        # path of replica input file (in its r{}/ directory)
        inpfile = "r%d/%s_%d.inp" % (replica, basename, cycle)
        # map state id to actual temperature value
        temperature = self.temperatures[stateid]
        # read template file
        tfile = self._openfile(template, "r")
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace("@n@",str(cycle))
        tbuffer = tbuffer.replace("@nm1@",str(cycle-1))
        tbuffer = tbuffer.replace("@temperature@",temperature)
        # write out
        ofile = self._openfile(inpfile, "w")
        ofile.write(tbuffer)
        ofile.close()
        # now the replica is ready for execution

<dl>
<dt>_doExchange_pair(self,repl_a,repl_b):</dt>
<dd>Performs an exchange attempt between replica 'repl_a' and 'repl_b'. This usually involves obtaining energies, structural information, thermodynamic and potential parameters as needed from output files or internal data structures. The required side effect of this function is the update of the thermodynamic state id tags of the two replicas. For example:</dd>
</dl>

    sid_a = self.status[repl_a]['stateid_current']
    sid_b = self.status[repl_b]['stateid_current']
    if <exchange successful>:
       self.status[repl_a]['stateid_current'] = sid_b
       self.status[repl_a]['stateid_current'] = sid_a

**Note** The `_doExchange_pair()` interface above based on the default Metropolis exchange algorithm will be soon replaced by a more efficient Gibbs sampling exchange algorithm already implemented for the Impact and AMBERUS extension modules once incompatibilities between the two implementations are resolved.

<dl>
<dt>_launchReplica(self,replica,cycle):</dt>
<dd>Instructs BigJob on how to launch a replica, which is typically a process common for all applications using the same MD engine. An example for AMBER (sander) is illustrated below. The routine is required to return the BigJob compute unit id of the replica being launched. Also, note that, technically, the submission of the replica to BigJob does not necessarily imply immediate execution; rather, the replica job is typically placed in a buffer area (see above) and will begin execution on when sufficient CPU resources on a compute node become available. For example:</dd>
</dl>

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
            "environment": [],
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

<dl>
<dt>_hasCompleted(self,repl,cy):</dt>
<dd>Optional. Returns 'True' if replica 'repl' has successfully completed cycle 'cy', otherwise 'False'. This is used to automatically relaunch replica runs that, for whatever reason, have failed to complete. Specification of this routine is optional because a default version, based on querying BigJob, exists in the core ASyncRE module. However the default routine is unable to detect the return status of a replica during a restart (see RE_SETUP above) because BigJob has no knowledge of it. In this case the default routine assumes success and returns 'True'. By providing a customizing routine ASyncRE will be able to detect failed replicas at restart and relaunch them automatically. For AMBER jobs, detecting the presence of a restart file is usually sufficient. For example:</dd>
</dl>

    def _hasCompleted(self,replica,cycle):
        rstfile = "r%d/%s_%d.rst7" % (replica, self.basename, cycle)
        if os.path.exists(rstfile):
            return True
        else:
            return False    

AMBER specifics:
----------------

Using different AMBER engines:

Currently, the user must make a choice between AMBER's two main MD engines, SANDER and PMEMD, in the ASyncRE input file. Invoking 'AMBER' will default to the SANDER engine, due to its broader capabilities ('AMBER-SANDER' and 'SANDER' are also recognized). The more performance tuned PMEMD can be requested via 'AMBER-PMEMD' or just 'PMEMD'. 

Executables compiled for use with MPI are automatically used when the BigJob setting 'SUBJOB_CORES' is greater than one (this also flags the proper 'SPMD' setting, see "Control settings" above). This choice is intentionally limited because AMBER MPI executables all exit with an error if mpirun (or equivalent commands) is called with fewer than two processors and there is currently no way (or reason) for AsyncRE to detect this error.

Neither pmemd.CUDA nor pmemd.CUDA.MPI are currently supported.

**AMBER control settings:**

<dl>
<dt>AMBER_GROUPFILE</dt>
<dd>Indicates an AMBER groupfile specifying input files (output files and other flags are ignored) defining each state and the starting structure of the replica initially in that state (ignored if RE_SETUP is 'no'). In this way, ASyncRE is perfectly compatible with input formats for traditional synchronous RE as implemented in AMBER. This is also currently the only way to specify different starting coordinates for each replica, a capability that can be especially important for replica exchange umbrella sampling.</dd>

<dt>AMBER_RESTRAINT_TEMPLATE</dt>
<dd>In replica exchange umbrella sampling simulations, this indicates which file to look in for determining restraints. If it is not specified, then a file using the ENGINE_INPUT_BASENAME with the extension '.RST' is expected. Note that it does not matter at all whether or not a restraint file is specified in the mdin file. ASyncRE will overwrite such information as needed.<dd>
</dl>

Multidimensional umbrella sampling control settings:
----------------------------------------------------

<dl>
<dt>FORCE_CONSTANTS and BIAS_POSITIONS</dt>
<dd>These settings should be strings of values indicating the respective quantities for harmonic biasing potentials. In one dimension, different states are delimited by a comma. In higher dimensions, each coordinate is delimited by a comma and each state is delimited by a colon. Examples:</dd>
</dl>

four states in one dimension:

    FORCE_CONSTANTS = '50.0,50.0,50.0,50.0'
    BIAS_POSITIONS = '1.0,1.5,2.0,2.5'

six states in two dimensions:

    FORCE_CONSTANTS = '5.0,5.0:5.0,5.0:5.0,5.0:5.0,5.0:5.0,5.0:5.0,5.0'
    BIAS_POSITIONS = '275.,275.:275.,280.:280.,275.:280.,280.:280.,285.:285.,280'

The order parameter units and bias form (i.e. Ubias = k(x-x0)^2 or (k/2)(x-x0)^2, FORCE_CONSTANTS --> k, BIAS_POSITIONS --> x0) will depend on the specific MD engine.
