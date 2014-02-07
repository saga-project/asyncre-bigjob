"""
Load Python libraries included with AmberTools (12 and later) and ascertain
the install status of AMBER.
"""
import os
import sys

__all__ = ['AMBERHOME', 'MKL_HOME',
           'AmberMdout', 'AmberParm', 'Rst7',
           'AMBER_SERIAL_EXES', 'AMBER_MPI_EXES', 'AMBER_CUDA_EXES',
           'KB'
           ]

# Check for environmental variables (required and optional).
#
AMBERHOME = os.getenv('AMBERHOME')
MKL_HOME = os.getenv('MKL_HOME')
if AMBERHOME == '':
    raise Exception('AMBERHOME is not set.')
else:
    try:
        sys.path.append(os.path.join(AMBERHOME,'bin'))
    except:
        raise Exception('Problem adding $AMBERHOME/bin to PYTHONPATH.')

# Load python libraries from Jason Swail. These are useful for file I/O of
# various AMBER file formats, particularly .parm7/.prmtop, .inpcrd/.rst7 
# (including netCDF), and mdout files.
#
try:
    from mdoutanalyzer.mdout import AmberMdout
    from chemistry.amber.readparm import AmberParm,Rst7
except ImportError:
    raise Exception('Could not load AMBER python libraries. These are'
                    ' only available in AmberTools12 and later.')

# Check the status of the AMBER executables.
#
def check_engine_status(engines):
    """Return false if any engine does not exist or is not executable."""
    for engine in engines:
        exe = os.path.join(AMBERHOME,'bin',engine)
        if not os.path.exists(exe) or not os.access(exe,os.X_OK):
            return False
    return True

AMBER_SERIAL_EXES = check_engine_status(['sander','pmemd'])
AMBER_MPI_EXES = check_engine_status(['sander.MPI','pmemd.MPI'])
AMBER_CUDA_EXES = check_engine_status(['pmemd.cuda','pmemd.cuda.MPI'])

# Constants 
#
# from src/sander/constants.F90:
BOLTZMANN = 1.380658e-23     # Boltzmann's constant in J/K
AVOGADRO = 6.0221367e23      # Avogadro's number in particles/mol
JPKC = 4184                  # Joules per kcal
KB = BOLTZMANN*AVOGADRO/JPKC # in kcal/mol-K
# from src/sander/runmd.F90:
BARKC = 1.6604345e4*4.184    # converts kcal/mol to bar?

