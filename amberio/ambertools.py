"""
Load Python libraries included with AmberTools (12 and later)

"""
import os
import sys

try:
    AMBERHOME = os.getenv('AMBERHOME')
except:
    raise Exception('AMBERHOME is not set.')
sys.path.append(os.path.join(AMBERHOME,'bin'))
try:
    from mdoutanalyzer.mdout import AmberMdout
    from chemistry.amber.readparm import rst7
except ImportError:
    raise Exception('Could not load AMBER python libraries. These are'
                    ' only available in AmberTools12 and later.')

# Constants 
# from src/sander/constants.F90:
BOLTZMANN = 1.380658e-23             # Boltzmann's constant in J/K
AVOGADRO  = 6.0221367e23             # Avogadro's number in particles/mol
JPKC      = 4184                     # Joules per kcal
KB        = BOLTZMANN*AVOGADRO/JPKC  # in kcal/mol-K
# from src/sander/runmd.F90:
BARKC     = 1.6604345e4*4.184        # converts kcal/mol to bar?

__author__ = 'Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>'

__all__ = ['AMBERHOME','KB','AmberMdout','rst7']
