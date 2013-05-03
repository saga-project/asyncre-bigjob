"""
FILE: amberrun.py - plugin for interacting with AMBER MD runs

DESCRIPTION: 

AUTHOR: Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>

REFERENCES: AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
"""
import os
import copy
import tempfile
import commands
import shutil

from amberio.ambertools import AMBERHOME,KB,AmberMdout
from amberio.mdin import read_amber_mdin
from amberio.rstr import read_amber_restraint

__all__ = ['read_amber_groupfile','parse_amber_args','amberrun_from_files',
           'AmberRunCollection','AmberRun','AMBER_EXTS']

# Recognized extensions for AMBER file formats
AMBER_EXTS = {'mdin': ['.in','.inp','.mdin'],
              'prmtop': ['.prmtop','.parm7','.parm'],
              'inpcrd': ['.inpcrd','.rst7','.mdcrd','.crd','rst']}
# file writing modes
_write_modes = ['-O', '-A']
# default filenames
_default_filenames = {'mdin': 'mdin', 'mdout': 'mdout', 'prmtop': 'prmtop', 
                      'inpcrd': 'inpcrd', 'restrt': 'restrt', 'ref': 'refc', 
                      'mdcrd': 'mdcrd', 'mdvel': 'mdvel', 'mden': 'mden',
                      'mdinfo': 'mdinfo', 'cpin': 'cpin', 'cpout': 'cpout', 
                      'cprestrt': 'cprestrt'}
# command line flags for each file type
_file_flags = {'mdin': '-i', 'mdout': '-o', 'prmtop': '-p', 'inpcrd': '-c', 
               'restrt': '-r', 'ref': '-ref', 'mdcrd': '-x', 'mdvel': '-v', 
               'mden': '-e', 'mdinfo': '-inf', 'cpin': '-cpin', 
               'cpout': '-cpout', 'cprestrt': '-cprestrt'}

def read_amber_groupfile(groupfile, engine='sander'):
    """
    Read an AMBER groupfile and return an AmberRunCollection object.

    REQUIRED ARGUMENTS:
    groupfile - AMBER groupfile, may contain comments, etc.

    OPTIONAL ARGUMENTS:
    engine - MD engine expected to run this groupfile (default = sander). This 
    is used to make the runs aware of default values.
    """
    run_list = AmberRunCollection()
    for line in open(groupfile,'r'):
        # ignore everything after #'s if they are present
        try:
            line = line[:line.index('#')]
        except ValueError:
            pass
        if len(line) > 0:
            run_list.append(parse_amber_args(line,engine))
    return run_list

def amberrun_from_files(basename, filenames=None, nruns=1, mode='', 
                        engine='sander'):
    """
    Create an AmberRunCollection of size 'nruns' by 1) searching the current
    directory for files with the prefix 'basename' and/or 2) a sequence of 
    filenames. In both cases the file types are determined by recognizing
    common file extensions.
    """
    # These are the bare minimum files that can define an AmberRun.
    files = {'mdin': None, 'prmtop': None, 'inpcrd': None}
    # First, try to match against the known filenames based off of their file 
    # extensions.
    if filenames is not None:
        for file in filenames:
            ext = os.path.splitext(file)[1]
            if ext in AMBER_EXTS['mdin']:
                files['mdin'] = file
            elif ext in AMBER_EXTS['prmtop']:
                files['prmtop'] = file
            elif ext in AMBER_EXTS['inpcrd']:
                files['inpcrd'] = file
    # Next, check in the current directory using the basename.
    for input in AMBER_EXTS.keys():
        for ext in AMBER_EXTS[input]:
            file = '%s%s'%(basename,ext)
            if os.path.exists(file) and files[input] is None:
                files[input] = file
    # There's no way to delineate reference and input coordinates w/o using a 
    # groupfile. Best guess is to assume they are the same.
    files['ref'] = files['inpcrd']
    # Give an error if no match of required files could be made.
    for input in AMBER_EXTS.keys():
        if files[input] is None:
            raise IOError('Insufficient specifications for AmberRunCollection.')
    # Define states with these file names
    run_list = AmberRunCollection()
    for n in range(nruns):
        run_list.append(AmberRun(mode=mode,engine=engine,**files))
    return run_list


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

    def set_engine(self, engine):
        for run in self: run.set_engine(engine)


def parse_amber_args(args, engine='sander'):
    """
    Parse AMBER command line arguments and return an AmberRun object.
    
    REQUIRED ARGUMENTS:
    args - a string of AMBER command line arguments

    OPTIONAL ARGUMENTS:
    engine - MD engine expected to run this groupfile (default = sander). This 
    is used to make the runs aware of default values.
    """
    parsed_mode = ''
    parsed_filenames = copy.copy(_default_filenames)
    args = args.split()
    for mode in _write_modes:
        if mode in args:
            parsed_mode = args.pop(args.index(mode))
    for file,flag in _file_flags.iteritems():
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
    def __init__(self, mode, engine, basename=None, **filenames):
        # ================================================
        # Determine filenames and build file based objects
        # ================================================
        # Set the default filenames, then rename based on given input
        self.filenames = copy.copy(_default_filenames)
        for file,filename in filenames.iteritems():
            self.filenames[file] = filename
        # Read mdin and construct an mdin object.
        self.mdin = read_amber_mdin(self.filenames['mdin'],engine)
        # Check that prmtop and inpcrd files actually exist.
        if not os.path.exists(self.filenames['prmtop']):
            raise Exception('Could not find AMBER prmtop file %s.'
                            %self.filenames['prmtop'])
        if not os.path.exists(self.filenames['inpcrd']):
            raise Exception('Could not find AMBER inpcrd file %s.'
                            %self.filenames['inpcrd'])
        # Define an AmberRestraint object if restraints are present
        self.rstr = None
        if self.has_restraints(): # Uses info from mdin object
            try:
                self.add_restraints(self.mdin.namelist_value('DISANG',None))
            except TypeError:
                print ('WARNING! nmropt > 0, but no DISANG input provided. A '
                       'restraint object was not initialized.')

        # ====================================================
        # Additional parameters that are not filenames/objects
        # ====================================================
        # Set the file writing mode
        self.mode = mode
        # If requested, set the basename of all output files.
        if basename is not None: self.set_basename(basename)
        # Set the MD engine (also effects defaults in mdin)
        self.set_engine(engine)

    def set_engine(self, engine):
        """Set the MD engine (and all of the inherent defaults).
        """
        self.engine = engine
        self.mdin.set_defaults(engine)

    def set_basename(self, basename):
        """Set the basename of all output files.
        """
        self.filenames['mdout'] = basename + '.out'
        self.filenames['restrt'] = basename + '.rst7'
        # Trajectory output format (ASCII or NetCDF)
        if self.mdin.namelist_value('ioutfm','cntrl') == 1:
            self.filenames['mdcrd'] = basename + '.nc'
        else:
            self.filenames['mdcrd'] = basename + '.crd'
        self.filenames['mdinfo'] = basename + '.info'

    def arguments(self):
        """
        Return a list of the command line arguments needed to launch a run.
        Arguments equal to the defaults will be omitted.
        """
        args = [self.mode]
        for file,filename in self.filenames.iteritems():
            if filename != _default_filenames[file]:
                args.extend([_file_flags[file],filename])
        return args

    def restart(self, is_restart=True):
        """Change whether or not this run is a restart.
        """
        if is_restart:
            self.mdin.set_namelist_value('irest',1,'cntrl')
            self.mdin.set_namelist_value('ntx',5,'cntrl')
        else:
            self.mdin.set_namelist_value('irest',0,'cntrl')

    def add_restraints(self, rstr_file, trace_file=None, print_step=None):
        # Set the &cntrl namelist to read nmr options.
        self.mdin.set_namelist_value('nmropt',1,'cntrl')
        # Add a &wt section with restraint output options.
        self.mdin.set_namelist_value('type',"'DUMPFREQ'",'wt',"'type'")
        if print_step is not None:
            self.mdin.set_namelist_value('istep1',print_step,'wt',"'DUMPFREQ'")
        self.mdin.set_namelist_value('type',"'END'",'wt',"'type'")
        # Add nmr variables that don't belong to a namelist.
        self.mdin.set_namelist_value('DISANG',rstr_file,None)
        if trace_file is not None:
            self.mdin.set_namelist_value('DUMPAVE',trace_file,None)
        self.mdin.set_namelist_value('LISTIN','POUT',None)
        # Construct an AmberRestraint object containing the restraint info.
        self.rstr = read_amber_restraint(rstr_file)
            
    def has_refc(self):
        return (self.mdin.namelist_value('ntr','cntrl') != 0)

    def has_restraints(self):
        return (self.mdin.namelist_value('nmropt','cntrl') != 0)

    def restrt_is_binary(self):
        return (self.mdin.namelist_value('ntxo','cntrl') != 1)

    def mdcrd_is_binary(self):
        return (self.mdin.namelist_value('ioutfm','cntrl') != 0)

    def write_amber_mdin(self, outfile=None, mode='w'):
        if outfile is None:
            outfile = self.filenames['mdin']
        self.mdin.write_amber_mdin(outfile,mode)
        
    def write_amber_restraint_file(self, outfile=None, title='', mode='w'):
        if outfile is None:
            outfile = self.mdin.namelist_value('DISANG',None)
        self.rstr.write_amber_restraint_file(outfile,title,mode)

    def snglpnt(self, inpcrd=None, disang='tmp.RST'):
        """
        Run a single AMBER energy calculation on a coordinate file. If none is
        given, use the last set of input coordinates for this run.
        """
        # Create a copy of the current AmberRun and set it to a snglpnt. 
        snglpnt_run = copy.deepcopy(self)
        snglpnt_run.mdin.set_namelist_value('imin',0,'cntrl')
        snglpnt_run.mdin.set_namelist_value('nstlim',0,'cntrl')
        snglpnt_run.mdin.set_namelist_value('irest',0,'cntrl')
        snglpnt_run.mdin.set_namelist_value('ntx',5,'cntrl')
        snglpnt_run.mdin.set_namelist_value('ntpr',1,'cntrl')
        snglpnt_run.mdin.set_namelist_value('ntwx',0,'cntrl')
        if inpcrd is not None:
            snglpnt_run.filenames['inpcrd'] = inpcrd
        rundir = os.getcwd()
        try:
            new_prmtop = os.path.abspath(snglpnt_run.filenames['prmtop'])
            new_inpcrd = os.path.abspath(inpcrd)
            new_refc = os.path.abspath(snglpnt_run.filenames['ref'])
            print 'prmtop:',new_prmtop
            print 'inpcrd:',new_inpcrd
            print 'refc:',new_refc

            tmpdir = tempfile.mkdtemp(prefix='amber-snglpnt-')
            os.chdir(tmpdir)
            snglpnt_run.filenames['prmtop'] = new_prmtop
            snglpnt_run.filenames['inpcrd'] = new_inpcrd
            snglpnt_run.filenames['ref'] = new_refc
            snglpnt_run.filenames['mdin'] = 'mdin'
            snglpnt_run.filenames['mdout'] = 'mdout'
            if snglpnt_run.has_restraints():
                snglpnt_run.mdin.set_namelist_value('DISANG',disang,None)
                snglpnt_run.write_amber_restraint_file(disang)
            snglpnt_run.write_amber_mdin()

            cmd = '%s %s'%(snglpnt_run.engine,' '.join(snglpnt_run.arguments()))
            print snglpnt_run.filenames
            print cmd
            commands.getstatusoutput(cmd)
            snglpnt_mdout = AmberMdout('mdout')
        finally:
            # Change back to the working directory and clean up output files.
            os.chdir(rundir)
            # shutil.rmtree(tmpdir)
        return snglpnt_mdout

    def energy(self, inpcrd=None, energy_component='EPtot'):
        """Calculate the energy of an AMBER coordinate file.
        """
        snglpnt_mdout = self.snglpnt(inpcrd)
        return snglpnt_mdout.data[energy_component][0]
        
    def reduced_energy(self, inpcrd=None, energy_component='EPtot'):
        """Calculate the "reduced" energy of an AMBER coordinate file.
        """
        snglpnt_mdout = self.snglpnt(inpcrd)
        if self.mdin.namelist_value('ntt','cntrl') != 0:
            temp0 = self.mdin.namelist_value('temp0','cntrl')
            beta = 1./(KB*temp0)
        else:
            raise Exception('The reduced energy is undefined for constant '
                            'energy ensembles.')
        energy = snglpnt_mdout.data[energy_component][0]
        # TODO: NpT and muVT reduced energies
        # pV = 0.
        # if snglpnt_mdout.data.has_key('VOLUME'):
        #     V = snglpnt_mdout.data['VOLUME'][0]
        #     pres0 = self.mdin.namelist_value('pres0','cntrl')
        #     pV = pres0*V
        return beta*energy
