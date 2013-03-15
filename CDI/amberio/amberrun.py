"""
FILE: amberrun.py - plugin for interacting with AMBER MD runs

DESCRIPTION: 

AUTHOR: Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>

REFERENCES: AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
"""
def ReadAmberGroupfile(groupfile, engine='sander'):
    """
    Read an AMBER groupfile and return an AmberRunCollection object.

    REQUIRED ARGUMENTS:
    groupfile - AMBER groupfile, may contain comments, etc.

    OPTIONAL ARGUMENTS:
    engine - MD engine expected to run this groupfile (default = sander). This 
    is used to make the runs aware of default values.
    """
    runList = AmberRunCollection()
    for line in open(groupfile,'r'):
        # ignore everything after #'s if they are present
        try:
            line = line[:line.index('#')]
        except ValueError:
            pass
        if len(line) > 0:
            runList.append(ParseAmberArguments(line,engine))
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

def ParseAmberArguments(args, engine='sander'):
    """
    Parse AMBER command line arguments and return an AmberRun object.
    
    REQUIRED ARGUMENTS:
    args - a string of AMBER command line arguments

    OPTIONAL ARGUMENTS:
    engine - MD engine expected to run this groupfile (default = sander). This 
    is used to make the runs aware of default values.
    """
    # write modes
    modes = ['-O', '-A']
    # input flags 
    file_flags = {'mdin' : '-i', 'mdout' : '-o', 'prmtop' : '-p',
                  'inpcrd' : '-c', 'restrt' : '-r', 'ref' : '-ref', 
                  'mdcrd' : '-x', 'mdvel' : '-v', 'mden' : '-e', 
                  'mdinfo' : '-inf', 'cpin' : '-cpin', 'cpout' : '-cpout',
                  'cprestrt' : '-cprestrt'}
    # start with default arguments
    parsed_mode = None
    parsed_filenames = { 'mdin' : 'mdin', 'mdout' : 'mdout', 
                         'prmtop' : 'prmtop', 'inpcrd' : 'inpcrd', 
                         'restrt' : 'restrt', 'ref' : 'refc', 
                         'mdcrd' : 'mdcrd', 'mdvel' : 'mdvel', 
                         'mdinfo' : 'mdinfo', 'cpin' : 'cpin', 
                         'cpout' : 'cpout', 'cprestrt' : 'cprestrt' }
    args = args.split()
    for mode in modes:
        if mode in args:
            parsed_mode = mode
    for file,flag in file_flags.iteritems():
        if flag in args:
            parsed_filenames[file] = args[args.index(flag)+1] 
    return AmberRun(mode=parsed_mode,engine=engine,**parsed_filenames)

class AmberRun(object):
    """
    An AMBER run is mostly defined by its input files (e.g. mdin, prmtop) but 
    also its output files (e.g. mdout, mdcrd). This object contains a dict
    (called filenames) of these. It also allows full manipulation of the mdin
    and restraint files as objects (called mdin and rstr). If no NMR restraints
    are present than rstr = None.

    (see mdin.py and rstr.py for details of these object)
    """
    def __init__(self, mode, engine, basename=None,**filenames):
        import os
        from mdin import ReadAmberMdinFile
        # ================================================
        # Determine filenames and build file based objects
        # ================================================
        # Set the default filenames, then rename based on given input
        self.filenames = {'mdin' : 'mdin', 'mdout' : 'mdout', 
                          'prmtop' : 'prmtop', 'inpcrd' : 'inpcrd', 
                          'restrt' : 'restrt', 'ref' : 'refc', 
                          'mdcrd' : 'mdcrd', 'mdvel' : 'mdvel', 
                          'mdinfo' : 'mdinfo', 'cpin' : 'cpin', 
                          'cpout' : 'cpout', 'cprestrt' : 'cprestrt' }
        for file,filename in filenames.iteritems():
            self.filenames[file] = filename
        # Read mdin and construct an mdin object.
        self.mdin = ReadAmberMdinFile(self.filenames['mdin'],engine)
        # Check that prmtop and inpcrd files actually exist.
        if not os.path.exists(self.filenames['prmtop']):
            raise Exception('Could not find AMBER prmtop file %s.'
                            %self.filenames['prmtop'])
        if not os.path.exists(self.filenames['inpcrd']):
            raise Exception('Could not find AMBER inpcrd file %s.'
                            %self.filenames['inpcrd'])
        # Define an AmberRestraint object if restraints are present
        self.rstr = None
        if self.HasRestraints(): # Uses info from mdin object
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
        # Trajectory output format (ASCII or NetCDF)
        if self.mdin.GetVariableValue('ioutfm','cntrl') == 1:
            self.filenames['mdcrd'] = basename + '.nc'
        else:
            self.filenames['mdcrd'] = basename + '.crd'
        self.filenames['mdinfo'] = basename + '.info'

    def GetArgs(self):
        """
        Return a list of the command line arguments needed to launch a run.
        Arguments equal to the defaults will be omitted.
        """
        file_flags = {'mdin' : '-i', 'mdout' : '-o', 'prmtop' : '-p',
                      'inpcrd' : '-c', 'restrt' : '-r', 'ref' : '-ref', 
                      'mdcrd' : '-x', 'mdvel' : '-v', 'mden' : '-e', 
                      'mdinfo' : '-inf', 'cpin' : '-cpin', 'cpout' : '-cpout',
                      'cprestrt' : '-cprestrt'}
        default_filenames = {'mdin' : 'mdin', 'mdout' : 'mdout', 
                             'prmtop' : 'prmtop', 'inpcrd' : 'inpcrd', 
                             'restrt' : 'restrt', 'ref' : 'refc', 
                             'mdcrd' : 'mdcrd', 'mdvel' : 'mdvel', 
                             'mdinfo' : 'mdinfo', 'cpin' : 'cpin', 
                             'cpout' : 'cpout', 'cprestrt' : 'cprestrt' }
        args = []
        if self.mode is not None:
            args.append(self.mode)
        for file,filename in self.filenames.iteritems():
            if filename != default_filenames[file]:
                args += [file_flags[file],filename]
        return args

    def SetRestart(self, isRestart=True):
        """Change whether or not this run is a restart.
        """
        self.mdin.SetRestart(isRestart)

    def AddRestraints(self, rstr_file, trace_file=None, print_step=None):
        from rstr import ReadAmberRestraintFile
        if trace_file is None:
            trace_file = self.mdin.GetVariableValue('DUMPAVE',None)
        if print_step is None:
            print_step = self.mdin.GetVariableValue('istep1','wt',
                                                        "'DUMPFREQ'")
        self.mdin.AddRestraints(rstr_file,trace_file,print_step)
        self.rstr = ReadAmberRestraintFile(rstr_file)
            
    def HasReferenceCoordinates(self):
        return (self.mdin.GetVariableValue('ntr','cntrl') == 1)

    def HasRestraints(self):
        return (self.mdin.GetVariableValue('nmropt','cntrl') != 0)

    def RunSinglePoint(self, inpcrd=None):
        """
        Run a single AMBER energy calculation on a coordinate file. If none is
        given, use the latest coordinate for this run.
        """
        import os,tempfile,copy,commands,shutil
        from ambertools import AmberMdout
        # 1) Make a temporary directory (tmpdir)
        cwd    = os.getcwd()
        tmpdir = tempfile.mkdtemp(prefix='amber-snglpnt-')
        os.chdir(tmpdir)
        # 2) Copy and modify the mdin object to run a single point energy
        snglpnt_mdin = copy.deepcopy(self.mdin)
        snglpnt_mdin.SetSinglePoint()
        # 3) Move to tmpdir, write a new mdin, get new path to other inputs
        snglpnt_mdin.WriteAmberMdinFile('mdin')
        prmtop = '%s/%s'%(cwd,self.filenames['prmtop'])
        if inpcrd is None: 
            inpcrd = '%s/%s'%(cwd,self.filenames['inpcrd'])
        else:
            inpcrd = '%s/%s'%(cwd,inpcrd)
        ref = '%s/%s'%(cwd,self.filenames['ref'])
        if self.HasRestraints():
            rstr = snglpnt_mdin.GetVariableValue('DISANG',None)
            snglpnt_mdin.SetVariableValue('DISANG',cwd + '/' + rstr,None)
        # 4) Call the MD engine and create an AmberMdout object from mdout
        cmd = '%s -O -p %s -c %s -ref %s'%(self.engine,prmtop,inpcrd,ref)
        commands.getstatusoutput(cmd)
        snglpnt_mdout = AmberMdout('mdout')
        # 5) Get out and clean up!
        os.chdir(cwd)
        shutil.rmtree(tmpdir)
        return snglpnt_mdout

    def Energy(self, inpcrd=None):
        snglpnt_mdout = self.RunSinglePoint(inpcrd)
        return snglpnt_mdout.data['EPtot'][0]
        
    def ReducedEnergy(self, inpcrd=None):
        """
        Use the current settings to calculate the "reduced" energy of an AMBER 
        coordinate file. If none is given, use the latest input coordinates for
        this run.
        """
        from ambertools import AMBERHOME,KB
        snglpnt_mdout = self.RunSinglePoint(inpcrd)

        temp0 = self.mdin.GetVariableValue('temp0','cntrl')
        beta = 1./(KB*temp0)
        energy = snglpnt_mdout.data['EPtot'][0]

        # TODO: NpT and muVT reduced energies
        # pV = 0.
        # if snglpnt_mdout.data.has_key('VOLUME'):
        #     V = snglpnt_mdout.data['VOLUME'][0]
        #     pres0 = self.mdin.GetVariableValue('pres0','cntrl')
        #     pV = pres0*V
        return beta*energy
