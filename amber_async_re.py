import os
import sys
import random
import math
import copy

import numpy as np

from pj_async_re import async_re_job
from amberio.ambertools import AMBERHOME,KB,rst7
from amberio.amberrun import *

__all__ = ['pj_amber_job',
           'AMBERHOME','KB',
           'SUPPORTED_AMBER_ENGINES','DISANG_NAME','DUMPAVE_EXT']

# TODO?: Cuda
SUPPORTED_AMBER_ENGINES = {'AMBER': 'sander', 'SANDER': 'sander', 
                            'AMBER-SANDER': 'sander', 'PMEMD': 'pmemd', 
                            'AMBER-PMEMD': 'pmemd'}
DISANG_NAME = 'restraint.RST' # hardcoded AMBER restraint file name
DUMPAVE_EXT = 'TRACE' # hardcoded file extension for restraint coordinates

class pj_amber_job(async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        
        # ===========================
        # Set up the AMBER executable
        # ===========================
        engine = self.keywords.get('ENGINE').upper()
        if SUPPORTED_AMBER_ENGINES.has_key(engine):
            engine = SUPPORTED_AMBER_ENGINES[engine]
        else:
            self._exit('Requested ENGINE (%s) is either invalid or not '
                       'currently supported.'%engine)
        if self.spmd == 'mpi' or int(self.keywords.get('SUBJOB_CORES')) > 1:
            self.spmd = 'mpi' # Always use MPI if SUBJOB_CORES > 1.
            engine += '.MPI'
        # Check that this file exists and is exectuable.
        self.exe = os.path.join(AMBERHOME,'bin',engine)
        if not os.path.exists(self.exe) or not os.access(self.exe,os.X_OK):
            self._exit('Could not find an executable: %s\nExpected it to'
                       ' be at %s'%(engine,self.exe))

        # ========================================================
        # Set up the general state/replica information - 2 methods
        # ========================================================
        # (1) If present, read the AMBER groupfile and define the states,
        if self.keywords.get('AMBER_GROUPFILE') is not None:
            groupfile = self.keywords.get('AMBER_GROUPFILE')
            self.states = read_amber_groupfile(groupfile,engine)
            self.nreplicas = len(self.states)
            if self.verbose:
                print ('Creating %d replicas from AMBER groupfile: %s'
                       %(self.nreplicas,groupfile))
        # (2) otherwise assume that the states can be inferred from the
        # extfiles and input from a specific application (e.g. umbrella 
        # sampling).
        else:
            if self.nreplicas is None:
                self._exit('Could not determine the replica count from the'
                           ' input provided (set NREPLICAS directly or provide'
                           ' an AMBER groupfile)')
            print 'basename',self.basename
            print 'extfiles',self.extfiles
            print 'nreplicas',self.nreplicas
            self.states = amberrun_from_files(self.basename,self.extfiles,
                                              self.nreplicas,'-O',engine)
            if self.verbose:
                print ('Creating %d replicas using the provided'
                       ' ENGINE_INPUT_EXTFILES and ENGINE_INPUT_BASENAME'
                       %self.nreplicas)

    def _buildInpFile(self, repl):
        """
        For a given replica:
        1) determine the current state 
        2) write a new mdin file (change to a restart input if cycle > 1)
        3) link to a new prmtop 
        4) link to a new ref file (as needed)
        5) link to the inpcrd from cycle = 0 if cycle = 1
        """
        sid = self.status[repl]['stateid_current']
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
            self.states[sid].mdin.set_namelist_value('DUMPAVE',trace_file,None)
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

        # Compute Unit (i.e. Job) description
        cpt_unit_desc = {
            'executable': self.exe,
            'environment': ['AMBERHOME=%s'%AMBERHOME],
            'arguments': args,
            'output': stdout,
            'error': stderr,   
            'working_directory': wdir,
            'number_of_processes': int(self.keywords.get('SUBJOB_CORES')),
            'spmd_variation': self.spmd,
            }

        compute_unit = self.pilotcompute.submit_compute_unit(cpt_unit_desc)
        return compute_unit
        
    def _hasCompleted(self, repl, cyc):
        """Returns True if an Amber replica has completed a cycle.

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
        return rst7(rst).coords

    def _state_params_are_same(self, variable, namelist):
        """
        Return False if any two states have different values of a 
        variable in the specified namelist. If all states have the same 
        value, then return that value.

        This routine can be useful if a particular exchange protocol 
        assumes that certain state parameters (e.g. temperature) are the
        same in all states.
        """
        value = self.states[0].mdin.namelist_value(variable,namelist)
        for state in self.states[1:]:
            this_value = state.mdin.namelist_value(variable,namelist)
            if this_value != value: 
                return False
        return value

    #
    # CURRENTLY BROKEN!
    #
    def _computeSwapMatrix(self, replicas, states):
        # U will be a sparse matrix, but is convenient bc the indices of the
        # rows and columns will always be the same.
        U = [[ 0. for j in range(self.nreplicas)] 
             for i in range(self.nreplicas)]
        for repl_i in replicas:
            cyc = self.status[repl_i]['cycle_current']  
            inpcrd_i = '%s/r%d/%s_%d.rst7'%(os.getcwd(),repl_i,self.basename,
                                            cyc)
            for sid_j in states:
                state_dir = '%s/r%d'%(os.getcwd(),sid_j)
                # energy of replica i in state j
                os.chdir('r%d'%sid_j)
                u_ji = self.states[sid_j].reduced_energy(inpcrd_i)
                U[sid_j][repl_i] = u_ji
                os.chdir('..')
        return U
