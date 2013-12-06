import os

import amberio.ambertools as at
from amberio.amberrun import read_amber_groupfile, amberrun_from_files
from pj_async_re import async_re_job, _exit

__all__ = ['pj_amber_job', 'amber_states_from_configobj',
           'extract_amber_coordinates', 'SUPPORTED_AMBER_ENGINES',
           'DISANG_NAME', 'DUMPAVE_EXT'
           ]

# TODO?: Cuda
SUPPORTED_AMBER_ENGINES = {'AMBER': 'sander', 'SANDER': 'sander', 
                           'AMBER-SANDER': 'sander', 'PMEMD': 'pmemd', 
                           'AMBER-PMEMD': 'pmemd'
                           }
DISANG_NAME = 'restraint.RST' # hardcoded AMBER restraint file name
DUMPAVE_EXT = 'TRACE' # hardcoded file extension for restraint coordinates

def which(program):
    def is_exe(filename):
        return os.path.exists(filename) and os.access(filename,os.X_OK)

    path,name = os.path.split(program)
    if path:
        if is_exe(program):
            return program
    for path in os.environ['PATH'].split(os.pathsep):
        exe = os.path.join(path,program)
        if is_exe(exe):
            return exe
    return None

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        try:
            engine = self.keywords.get('ENGINE').upper()
            engine = SUPPORTED_AMBER_ENGINES[engine]
        except KeyError:
            _exit('Requested ENGINE (%s) is either invalid or not '
                  'currently supported.'%self.keywords.get('ENGINE'))

        if int(self.keywords.get('SUBJOB_CORES')) > 1:
            self.spmd = 'mpi'
            engine = '%s.MPI'%engine
            if not at.AMBER_MPI_EXES:
                _exit('Cannot find AMBER MPI executables. Are these compiled '
                      'and in AMBERHOME/bin?')
        else:
            self.spmd = 'single'
            if not at.AMBER_SERIAL_EXES:
                _exit('Cannot find AMBER serial executables. Are these '
                      'compiled and in AMBERHOME/bin?')
      
        self.exe = os.path.join(at.AMBERHOME,'bin',engine)
        self.states = amber_states_from_configobj(self.keywords,self.verbose)
        self.nreplicas = len(self.states)

    def _buildInpFile(self, repl, state = None):
        """
        For a given replica:
        1) determine the current state 
        2) write a new mdin file (change to a restart input if cycle > 1)
        3) link to a new prmtop 
        4) link to a new ref file (as needed)
        5) link to the inpcrd from cycle = 0 if cycle = 1
        """
        if state is None:
            sid = self.status[repl]['stateid_current']
        else:
            sid = state
        cyc = self.status[repl]['cycle_current']
        # Make a copy of one of the existing AmberRun state templates.
        title = ' replica %d : state %d : cycle %d'%(repl,sid,cyc)
        self.states[sid].mdin.title = title
        # Modify the template as appropriate.
        if cyc > 1: 
            self.states[sid].restart()
        if self.states[sid].has_restraints:
            rstr_title = title
            rstr_file = 'r%d/%s'%(repl,DISANG_NAME)
            self.states[sid].rstr.write_amber_restraint_file(rstr_file,title)
            trace_file = '%s_%d.%s'%(self.basename,cyc,DUMPAVE_EXT)
            self.states[sid].mdin.nmr_vars['DUMPAVE'] = trace_file
        self.states[sid].mdin.write_amber_mdin('r%d/mdin'%repl)
        # Links
        prmtop = self.states[sid].filenames['prmtop']
        self._linkReplicaFile('prmtop',prmtop,repl)
        if self.states[sid].has_refc:
            refc = self.states[sid].filenames['ref']
            self._linkReplicaFile('refc',refc,repl)
        if cyc == 1:
            inpcrd = self.states[sid].filenames['inpcrd']
            self._linkReplicaFile('%s_0.rst7'%self.basename,inpcrd,repl) 

    def _launchReplica(self, repl, cyc):
        """Launch an AMBER sub-job using pilot-job. 

        The input files for AMBER that define a state are assumed to be 
        the default names mdin, prmtop, and refc. These files are always
        re-written or symlinked to in _buildInpFile().
        """
        # Working directory for this replica
        wdir = '%s/r%d'%(os.getcwd(),repl)

        # Cycle dependent input and output file names
        inpcrd = '%s_%d.rst7'%(self.basename,cyc-1)
        mdout  = '%s_%d.out'%(self.basename,cyc)
        mdcrd  = '%s_%d.nc'%(self.basename,cyc)
        restrt = '%s_%d.rst7'%(self.basename,cyc)
        stdout = '%s_%d.log'%(self.basename,cyc)
        stderr = '%s_%d.err'%(self.basename,cyc)

        args = ['-O','-c',inpcrd,'-o',mdout,'-x',mdcrd,'-r',restrt]

        amber_env = ['AMBERHOME=%s'%at.AMBERHOME, 'MKL_HOME=%s'%at.MKL_HOME]
        amber_env.extend(self.engine_environment)

        script_name = 'run'
        run_script = open('%s/%s'%(wdir,script_name),'w')
        for env in amber_env:
            run_script.write('export %s\n'%env)
        run_script.write('EXE=%s\n\n'%self.exe)
        run_script.write('cd %s\n'%wdir)
        run_script.write('$EXE %s\n\n'%(' '.join(args)))
        # run_script.write('cd ..\n')
        # run_script.write('python calc_all_us_state_energies.py %s %s %d\n'
        #                  %(self.command_file,'%s/%s'%(wdir,restrt),repl))
        run_script.close()

        # Compute Unit (i.e. Job) description
        bash = which('bash')
        cpt_unit_desc = {
            'executable': '%s %s'%(bash,script_name),
            'environment': [],
            'arguments': [],
            'output': stdout,
            'error': stderr,   
            'working_directory': wdir,
            'number_of_processes': int(self.keywords.get('SUBJOB_CORES')),
            'spmd_variation': 'single',
            }

        compute_unit = self.pilotcompute.submit_compute_unit(cpt_unit_desc)
        return compute_unit
        
    def _hasCompleted(self, repl, cyc):
        """
        Return true if an AMBER replica has completed a cycle.

        Basically checks if the restart file exists.
        """
        # TODO: Parse the output file and look for more sure signs of 
        #       completion?
        rst = 'r%d/%s_%d.rst7'%(repl,self.basename,cyc)
        if os.path.exists(rst):
            return async_re_job._hasCompleted(self,repl,cyc)
        else:
            return False

    def _extractLastCoordinates(self, repl):
        """
        Return a 3N list of coordinates from the last restart (rst7) file 
        of a given replica.
        """
        cyc = self.status[repl]['cycle_current']
        rst = 'r%d/%s_%d.rst7'%(repl,self.basename,cyc)
        return at.rst7(rst).coords

def extract_amber_coordinates(replica, cycle, basename):
    restrt_name = 'r%d/%s_%d.rst7'%(replica,basename,cycle)
    return at.rst7(restrt_name).coords

def amber_states_from_configobj(keywords, verbose=False):
    """Return an AmberRunCollection from an ASyncRE command file."""
    # keywords = ConfigObj(command_file)
    try:
        engine = keywords.get('ENGINE').upper()
        engine = SUPPORTED_AMBER_ENGINES[engine]
    except KeyError:
        _exit('Requested ENGINE (%s) is either invalid or not currently '
              'supported.'%keywords.get('ENGINE'))        
    
    # Set up the general state/replica information - 2 methods
    # 
    # (1) If present, read the AMBER groupfile and define the states,
    if keywords.get('AMBER_GROUPFILE') is not None:
        groupfile = keywords.get('AMBER_GROUPFILE')
        states = read_amber_groupfile(groupfile,engine)
        if verbose:
            print ('Created %d replicas from AMBER groupfile: %s'
                   %(len(states),groupfile))
    # (2) otherwise assume that the states can be inferred from the extfiles 
    # and input from a specific application (e.g. umbrella sampling).
    else:
        basename = keywords.get('ENGINE_INPUT_BASENAME')
        extfiles = keywords.get('ENGINE_INPUT_EXTFILES')
        if extfiles is not None and extfiles != '':
            extfiles = extfiles.split(',')
        else:
            extfiles = None
        nreplicas = int(keywords.get('NREPLICAS'))
        if nreplicas is None:
            _exit('Could not determine the replica count from the input '
                  'provided (set NREPLICAS directly or provide an AMBER '
                  'groupfile)')
        try:
            states = amberrun_from_files(basename,extfiles,nreplicas,'-O',
                                         engine)
        except IOError:
            _exit('Problem creating replicas, not enough information?')
        if verbose:
            print ('Created %d replicas using the provided '
                   'ENGINE_INPUT_EXTFILES and ENGINE_INPUT_BASENAME'
                   %nreplicas)
    return states
        
