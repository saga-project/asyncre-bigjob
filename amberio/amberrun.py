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
from amberio.mdin import AmberMdin 
from amberio.rstr import AmberRestraint

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
        run_list.append(AmberRun(engine,mode,None,**files))
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

    def state_params_are_same(self, namelist, variable):
        """
        Return false if any two states have different values of a variable in 
        the specified namelist. If all states have the same value, then return 
        that value.

        This routine can be useful if, for example, a particular replica 
        exchange protocol assumes that certain state parameters (e.g. 
        temperature) are the same in all states.
        """
        value = self[0].mdin.__getattribute__(namelist)[variable]
        for state in self[1:]:
            this_value = state.mdin.__getattribute__(namelist)[variable]
            if this_value != value:
                return False
        return value

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
    return AmberRun(engine,parsed_mode,None,**parsed_filenames)


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
        self.mdin = AmberMdin.from_mdin(self.filenames['mdin'],engine)
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
                self.add_restraints(self.mdin.nmr_vars['DISANG'])
            except (TypeError,KeyError):
                self.rstr = None
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
            if self.mdin.cntrl['ioutfm'] == 1:
                self.filenames['mdcrd'] = basename + '.nc'
            else:
                self.filenames['mdcrd'] = basename + '.crd'
            self.filenames['mdinfo'] = basename + '.info'

    def __getattribute__(self, name):
        if name == 'has_refc':
            return (self.mdin.cntrl['ntr'] != 0)
        elif name == 'has_restraints':
            return (self.mdin.cntrl['nmropt'] != 0)
        elif name == 'restrt_is_binary':
            return (self.mdin.cntrl['ntxo'] == 2)
        elif name == 'mdcrd_is_binary':
            return (self.mdin.cntrl['ioutfm'] == 1)
        elif name == 'arguments':
            args = [object.__getattribute__(self,'mode')]
            for file,filename in \
                    object.__getattribute__(self,'filenames').iteritems():
                if filename != DEFAULT_FILENAMES[file]:
                    args.extend([FILE_FLAGS[file],filename])
            return args
        else:
            return object.__getattribute__(self,name)

    def restart(self, is_restart=True):
        """Change whether or not this run is a restart."""
        if is_restart:
            self.mdin.cntrl['irest'] = 1
            self.mdin.cntrl['ntx'] = 5
        else:
            self.mdin.cntrl['irest'] = 0

    def binary_formatting(self, is_binary=True):
        """Change whether or not this run uses netCDF formats."""
        if is_binary:
            self.mdin.cntrl['ntxo'] = 2
            self.mdin.cntrl['ioutfm'] = 1
        else:
            self.mdin.cntrl['ntxo'] = 1
            self.mdin.cntrl['ioutfm'] = 0
           
    def add_restraints(self, rstr_file, trace_file=None, print_step=None):
        """
        Add or modify restraints:
        
        (1) Set the &cntrl namelist to read nmr options.
        (2) Set the appropiate nmr input files.
        (3) Add or modify wt namelists for output (optional).
        (4) Set the appropriate nmr output files (optional).
        """
        self.mdin.cntrl['nmropt'] = 1
        self.mdin.nmr_vars['DISANG'] = rstr_file
        self.mdin.nmr_vars['LISTIN'] = 'POUT'
        self.rstr = AmberRestraint.from_disang(rstr_file)

        if print_step is not None:
            self.mdin.modify_or_add_wt("'DUMPFREQ'",0,**{'istep1': print_step})
        if trace_file is not None:
            self.mdin.nmr_vars['DUMPAVE'] = trace_file


    def write_amber_mdin(self, outfile=None, mode='w'):
        """Write a new mdin file. Default overwrite currents file."""
        if outfile is None:
            outfile = self.filenames['mdin']
        self.mdin.write_amber_mdin(outfile,mode)
        
    def write_restraint_file(self, outfile=None, title='', mode='w'):
        """Write a new restraint file. Default overwrites current file."""
        if outfile is None:
            outfile = self.mdin.nmr_vars['DISANG']
        self.rstr.title = title
        self.rstr.write(outfile,mode)

    def snglpnt(self, inpcrd=None, disang='tmp.RST'):
        """
        Run a single AMBER energy calculation on a coordinate file. If 
        none is given, use the last set of input coordinates for this 
        run.
        """
        # Create a copy of the current AmberRun and set it to a snglpnt. 
        snglpnt_run = copy.deepcopy(self)
        snglpnt_run.mdin.cntrl['imin'] = 1
        snglpnt_run.mdin.cntrl['nstlim'] = 0
        snglpnt_run.mdin.cntrl['irest'] = 0
        snglpnt_run.mdin.cntrl['ntx'] = 5
        snglpnt_run.mdin.cntrl['ntpr'] = 1
        snglpnt_run.mdin.cntrl['ntwx'] = 0
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
                snglpnt_run.mdin.nmr_vars['DISANG'] = disang
                snglpnt_run.write_restraint_file(disang)
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
        if self.mdin.cntrl['ntt'] != 0:
            beta = 1./(at.KB*self.mdin.cntrl['temp0'])
        else:
            raise Exception('The reduced energy is undefined for constant '
                            'energy ensembles.')
        energy = snglpnt_mdout.data[energy_component][0]
        # TODO: NpT and muVT reduced energies
        # pV = 0.
        # if snglpnt_mdout.data.has_key('VOLUME'):
        #     V = snglpnt_mdout.data['VOLUME'][0]
        #     pV = self.mdin.cntrl['pres0']*V
        return beta*energy
