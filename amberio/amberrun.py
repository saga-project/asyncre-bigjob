"""
Plugin for interacting with AMBER MD runs via Python

References: 
    AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
"""
import os
import copy
import tempfile
import commands
import shutil

import amberio.ambertools as at
from amberio.mdin import read_amber_mdin
from amberio.rstr import read_amber_restraint

__author__ = 'Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>'

__all__ = ['read_amber_groupfile','parse_amber_args','amberrun_from_files',
           'AmberRunCollection','AmberRun',
           'AMBER_EXTS','WRITE_MODES','DEFAULT_FILENAMES','FILE_FLAGS']

# Recognized extensions for AMBER file formats
AMBER_EXTS = {'mdin': ['.in','.inp','.mdin'],
              'prmtop': ['.prmtop','.parm7','.parm'],
              'inpcrd': ['.inpcrd','.rst7','.mdcrd','.crd','rst']}
# file writing modes
WRITE_MODES = ['-O', '-A']
# default filenames
DEFAULT_FILENAMES = {'mdin': 'mdin', 'mdout': 'mdout', 'prmtop': 'prmtop', 
                     'inpcrd': 'inpcrd', 'restrt': 'restrt', 'ref': 'refc', 
                     'mdcrd': 'mdcrd', 'mdvel': 'mdvel', 'mden': 'mden',
                     'mdinfo': 'mdinfo', 'cpin': 'cpin', 'cpout': 'cpout', 
                     'cprestrt': 'cprestrt'}
# command line flags for each file type
FILE_FLAGS = {'mdin': '-i', 'mdout': '-o', 'prmtop': '-p', 'inpcrd': '-c', 
              'restrt': '-r', 'ref': '-ref', 'mdcrd': '-x', 'mdvel': '-v', 
              'mden': '-e', 'mdinfo': '-inf', 'cpin': '-cpin', 
              'cpout': '-cpout', 'cprestrt': '-cprestrt'}

def read_amber_groupfile(groupfile, engine='sander'):
    """Read an AMBER groupfile and return an AmberRunCollection object."""
    run_list = AmberRunCollection()
    for line in open(groupfile,'r'):
        # Ignore everything after #'s if they are present.
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
    Create an AmberRunCollection of size 'nruns' by 1) searching the 
    current directory for files with the prefix 'basename' and/or 2) a 
    sequence of filenames. In both cases the file types are determined
    by recognizing common file extensions.
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
        run_list.append(AmberRun(engine,mode,**files))
    return run_list


class AmberRunCollection(list):
    """A list of AmberRun objects."""
    def __init__(self, *runs):
        list.__init__(self)
        for run in runs: self.extend(run)

    def append(self, item):
        if not isinstance(item,AmberRun):
            raise TypeError('AmberRunCollections must contain AmberRuns!')
        if hasattr(self,'engine'):
            item.engine = self.engine
        list.append(self,item)

    def extend(self, items):
        if hasattr(items, '__iter__'):
            for item in items:
                self.append(item)

    def __setattr__(self, name, value):
        list.__setattr__(self,name,value)
        if name == 'engine':
            for run in self:
                run.engine = value
    
    def write_amber_groupfile(self, outfile, mode='w'):
        closeafter = True
        if hasattr(outfile,'write'):
            closeafter = False
        elif isinstance(outfile,str):
            outfile = open(outfile,mode)
        else:
            raise TypeError("'outfile' must be either a string or file object.")
        for run in self:
            outfile.write('%s\n'%' '.join(run.arguments))
        if closeafter:
            outfile.close()

def parse_amber_args(args, engine='sander'):
    """Parse AMBER command line arguments and return an AmberRun object."""
    arg_list = args.split()
    parsed_mode = ''
    for mode in WRITE_MODES:
        try:
            parsed_mode = arg_list.pop(args.index(mode))
        except ValueError:
            pass

    parsed_filenames = copy.copy(DEFAULT_FILENAMES)
    for file,flag in FILE_FLAGS.iteritems():
        try:
            parsed_filenames[file] = arg_list[arg_list.index(flag)+1]
        except ValueError:
            pass
    return AmberRun(engine,parsed_mode,**parsed_filenames)


class AmberRun(object):
    """
    An AMBER run is mostly defined by its input files (e.g. mdin, 
    prmtop) but also its output files (e.g. mdout, mdcrd). This object 
    contains a dict (called filenames) of these. It also allows full 
    manipulation of the mdin and restraint files as objects (called mdin
    and rstr). If no NMR restraints are present than rstr = None.

    (see mdin.py and rstr.py for details of these object)
    """
    def __init__(self, engine, mode='', basename=None, **filenames):
        # Determine filenames and build file based objects
        # 
        # Set the default filenames, then rename based on the input.
        self.filenames = copy.copy(DEFAULT_FILENAMES)
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
        # Define an AmberRestraint object if restraints are present.
        self.rstr = None
        if self.has_restraints: # Use info from the mdin file.
            try:
                self.add_restraints(self.mdin.namelist_value('DISANG',None))
            except TypeError:
                print ('WARNING! nmropt > 0, but no DISANG input provided. A '
                       'restraint object was not initialized.')

        # Additional parameters that are not filenames/objects
        #
        self.mode = mode # file writing mode
        self.basename = basename # output file basename
        self.engine = engine # MD engine (also effects defaults in mdin)

    def __setattr__(self, name, value):
        object.__setattr__(self,name,value)
        if name == 'engine':
            self.mdin.engine = value
        elif name == 'basename' and value is not None:
            self.filenames['mdout'] = basename + '.out'
            self.filenames['restrt'] = basename + '.rst7'
            # Trajectory output format (ASCII or NetCDF)
            if self.mdin.namelist_value('ioutfm','cntrl') == 1:
                self.filenames['mdcrd'] = basename + '.nc'
            else:
                self.filenames['mdcrd'] = basename + '.crd'
            self.filenames['mdinfo'] = basename + '.info'

    def __getattr__(self, name):
        if name == 'has_refc':
            return (self.mdin.namelist_value('ntr','cntrl') != 0)
        elif name == 'has_restraints':
            return (self.mdin.namelist_value('nmropt','cntrl') != 0)
        elif name == 'restrt_is_binary':
            return (self.mdin.namelist_value('ntxo','cntrl') != 1)
        elif name == 'mdcrd_is_binary':
            return (self.mdin.namelist_value('ioutfm','cntrl') != 0)
        elif name == 'arguments':
            args = [self.mode]
            for file,filename in self.filenames.iteritems():
                if filename != DEFAULT_FILENAMES[file]:
                    args.extend([FILE_FLAGS[file],filename])
            return args
        else:
            object.__get__attr(self,name)

    def restart(self, is_restart=True):
        """Change whether or not this run is a restart."""
        if is_restart:
            self.mdin.set_namelist_value('irest',1,'cntrl')
            self.mdin.set_namelist_value('ntx',5,'cntrl')
        else:
            self.mdin.set_namelist_value('irest',0,'cntrl')

    def binary_formatting(self, is_binary=True):
        """Change whether or not this run uses netCDF formats."""
        if is_binary:
            self.mdin.set_namelist_value('ntxo',2,'cntrl')
            self.mdin.set_namelist_value('ioutfm',1,'cntrl')
        else:
            self.mdin.set_namelist_value('ntxo',1,'cntrl')
            self.mdin.set_namelist_value('ioutfm',0,'cntrl')

    def add_restraints(self, rstr_file, trace_file=None, print_step=None):
        """Add or modify restraints."""
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

    def write_amber_mdin(self, outfile=None, mode='w'):
        """Write a new mdin file. Default overwrite currents file."""
        if outfile is None:
            outfile = self.filenames['mdin']
        self.mdin.write_amber_mdin(outfile,mode)
        
    def write_amber_restraint_file(self, outfile=None, title='', mode='w'):
        """Write a new restraint file. Default overwrites current file."""
        if outfile is None:
            outfile = self.mdin.namelist_value('DISANG',None)
        self.rstr.write_amber_restraint_file(outfile,title,mode)

    def snglpnt(self, inpcrd=None, disang='tmp.RST'):
        """
        Run a single AMBER energy calculation on a coordinate file. If 
        none is given, use the last set of input coordinates for this 
        run.
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
            if snglpnt_run.has_restraints:
                snglpnt_run.mdin.set_namelist_value('DISANG',disang,None)
                snglpnt_run.write_amber_restraint_file(disang)
            snglpnt_run.write_amber_mdin()

            cmd = '%s %s'%(snglpnt_run.engine,' '.join(snglpnt_run.arguments()))
            print snglpnt_run.filenames
            print cmd
            commands.getstatusoutput(cmd)
            snglpnt_mdout = at.AmberMdout('mdout')
        finally:
            # Change back to the working directory and clean up output files.
            os.chdir(rundir)
            # shutil.rmtree(tmpdir)
        return snglpnt_mdout

    def energy(self, inpcrd=None, energy_component='EPtot'):
        """Calculate the energy of an AMBER coordinate file."""
        snglpnt_mdout = self.snglpnt(inpcrd)
        return snglpnt_mdout.data[energy_component][0]
        
    def reduced_energy(self, inpcrd=None, energy_component='EPtot'):
        """Calculate the "reduced" energy of an AMBER coordinate file."""
        snglpnt_mdout = self.snglpnt(inpcrd)
        if self.mdin.namelist_value('ntt','cntrl') != 0:
            temp0 = self.mdin.namelist_value('temp0','cntrl')
            beta = 1./(at.KB*temp0)
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
