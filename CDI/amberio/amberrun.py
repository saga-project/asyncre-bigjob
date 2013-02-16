"""
FILE: amberrun.py - plugin for interacting with AMBER MD runs

DESCRIPTION: 

AUTHOR: Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>

REFERENCES: AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
"""
def ReadAmberGroupfile(groupfile):
    """Read an AMBER groupfile and return an AmberRunCollection.
    """
    runList = AmberRunCollection()
    for line in open(groupfile,'r'):
        # ignore everything after #'s if they are present
        try:
            line = line[:line.index('#')]
        except ValueError:
            pass
        if len(line) > 0:
            mode = None
            filenames = {}
            tokens = line.split()
            for i,token in enumerate(tokens):
                if token in ['-O','-A']: mode = token
                if token == '-i': filenames['mdin'] = tokens[i+1]
                if token == '-o': filenames['mdout'] = tokens[i+1]
                if token == '-p': filenames['prmtop'] = tokens[i+1]
                if token == '-c': filenames['inpcrd'] = tokens[i+1]
                if token == '-r': filenames['restrt'] = tokens[i+1]
                if token == '-ref': filenames['ref'] = tokens[i+1]
                if token == '-x': filenames['mdcrd'] = tokens[i+1]
                if token == '-inf': filenames['mdinfo'] = tokens[i+1]
            runList.append(AmberRun(mode=mode,**filenames))
    return runList

class AmberRunCollection(list):
    """A list of AmberRun objects.
    """
    def __init__(self,*runs):
        list.__init__(self)
        for run in runs: self.extend(run)

    def append(self, item):
        if not isinstance(item, AmberRun):
            raise TypeError('AmberRunCollections must contain AmberRuns!')
        list.append(self, item)

    def extend(self, items):
        if hasattr(items, '__iter__'):
            for item in items:
                self.append(item)

    def SetEngine(self, engine):
        for run in self: run.SetEngine(engine)

class AmberRun(object):
    """
    An AMBER run is mostly defined by its input files (e.g. mdin, prmtop) but 
    also its output files (e.g. mdout, mdcrd). This object contains a dict
    (called filenames) of these. It also allows full manipulation of the mdin
    and restraint files as objects (called mdin and rstr). If no NMR restraints
    are present than rstr = None.

    (see mdin.py and rstr.py for details of these object)
    """
    def __init__(self,mode=None,basename=None,engine='sander',**filenames):
        import copy
        from mdin import ReadAmberMdinFile
        # ==================================
        # Dictionaries for referencing files 
        # ==================================
        # Input flags (used to write arguments)
        self.file_flags = {'mdin' : '-i', 'mdout' : '-o',
                           'prmtop' : '-p', 'inpcrd' : '-c', 
                           'restrt' : '-r', 'ref' : '-ref', 
                           'mdcrd' : '-x', 'mdinfo' : '-inf'}
        # Set the default filenames of required files.
        self.default_filenames = {'mdin' : 'mdin', 'mdout' : 'mdout', 
                                  'prmtop' : 'prmtop', 'inpcrd' : 'inpcrd', 
                                  'restrt' : 'restrt', 'ref' : 'refc', 
                                  'mdcrd' : 'mdcrd', 'mdinfo' : 'mdinfo'}
        # ================================================
        # Determine filenames and build file based objects
        # ================================================
        self.filenames = copy.deepcopy(self.default_filenames)
        # Set new filenames as specified by the input.
        for file in filenames.keys():
            if file in self.filenames.keys():
                self.filenames[file] = filenames[file]
        # Read mdin and determine certain attributes of this run.
        self.mdin = ReadAmberMdinFile(self.filenames['mdin'])
        # Trajectory output format (ASCII or NetCDF)
        self.useBinTraj = False
        if self.mdin.GetVariableValue('ioutfm','cntrl') == 1: 
            self.useBinTraj = True
        # Start/Restart
        self.isRestart = False
        if self.mdin.GetVariableValue('irest','cntrl') == 1: 
            self.isRestart = True
        # NMRopt restraints
        self.AddRestraints(self.mdin.GetVariableValue('DISANG',None))
        # ====================================================
        # Additional parameters that are not filenames/objects
        # ====================================================
        # Set the file writing mode O(verwrite), A(ppend), or None
        self.mode = mode
        # If requested, set the basename of all output files.
        if basename is not None: self.SetBasename(basename)
        # Set the MD engine (also effects defaults in mdin)
        self.SetEngine(engine)

    def SetEngine(self, engine):
        """Set the MD engine (and all of the inherent defaults).
        """
        self.engine = engine
        self.mdin.SetDefaults(engine)

    def SetBasename(self, basename):
        """Set the basename of all output files.
        """
        self.filenames['mdout'] = basename + '.out'
        self.filenames['restrt'] = basename + '.rst7'
        if self.useBinTraj: 
            self.filenames['mdcrd'] = basename + '.nc'
        else:
            self.filenames['mdcrd'] = basename + '.crd'
        self.filenames['mdinfo'] = basename + '.info'

    def GetArgs(self):
        """
        Return a list of the command line arguments needed to launch a run.
        Arguments equal to the defaults will be omitted.
        """
        args = []
        if self.mode is not None: args.append(self.mode)
        for file in self.filenames.keys():
            if self.filenames[file] != self.default_filenames[file]:
                args += [self.file_flags[file],self.filenames[file]]
        return args

    def AddRestraints(self, rstr_file, trace_file=None, print_step=None):
        try:
            from rstr import ReadAmberRestraintFile
            if trace_file is None:
                trace_file = self.mdin.GetVariableValue('DUMPAVE',None)
            if print_step is None:
                print_step = self.mdin.GetVariableValue('istep1','wt',
                                                        "'DUMPFREQ'")
            self.mdin.AddRestraints(rstr_file,trace_file,print_step)
            self.rstr = ReadAmberRestraintFile(rstr_file)
        except IOError:
            self.rstr = None
            # Let this slide so that AMBER ultimately reports the error.
            # (It might be a dummy file to establish the DISANG variable.)
            print ('WARNING! Unable to read AMBER restraint file: %s\n'
                   '         Continuing anyway...'%rstr_file)
            
    # def Energy(self, crd_file):
    #     import tempfile,shutil,os,copy,commands
    #     # - Make a temporary directory
    #     # - Make a modified mdin object
    #     # - Move to the tmp dir and get the exe and args
    #     # - Call the MD engine
    #     # - Parse the output and return the energy
    #     cwd    = os.getcwd()
    #     tmpdir = tempfile.mkdtemp(prefix='amber-snglpnt-')
    #     os.chdir(tmpdir)
    #     
    #     snglpnt_mdin = copy.deepcopy(self.mdin)
    #     snglpnt_mdin.
    #   
    #     os.chdir(cwd)
    #     shutil.rmtree(tmpdir)

        
    #     return 0.
