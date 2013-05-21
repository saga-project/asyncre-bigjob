"""
Plugin for I/O of AMBER coordinate files

"""
try:
    from netCDF4 import Dataset
except ImportError:
    print ('WARNING: Unable to find netCDF4. Reading of restart files in '
           'netCDF format will not be available.')

from amberio.ambertools import rst7 as _rst7_parmed

__author__ = 'Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>'

__all__ = ['rst7']


class rst7(object):
    """Create an rst7 file object from an AMBER restart file.

    Both ASCII and netCDF formats are supported.
    """
    def __init__(self, inpcrd_name, is_binary=True):
        self.filename = inpcrd_name
        self.valid = False
        self.title = ''
        self.box = []
        self.hasbox = False
        self.vels = []
        self.hasvels = False

        if is_binary:
            # Try to read the restart file in netCDF format. If that fails,
            # fall back to reading a formatted restart file.
            try:
                tmp = Dataset(inpcrd_name,'r')
                self.valid = True
                self.time = tmp.variables['time'][:]
                self.coords = tmp.variables['coordinates'][:].ravel().tolist()
                self.natom = len(self.coords)/3
                if tmp.variables.has_key('velocities'):
                    self.hasvels = True
                    self.vels = tmp.variables['velocities'][:].ravel().tolist()
                if (tmp.variables.has_key('cell_lengths') and
                    tmp.variables.has_key('cell_angles')):
                    self.hasbox = True
                    self.box = tmp.variables['cell_lengths'][:].tolist()
                    self.box.extend(tmp.variables['cell_angles'][:].tolist())
            except RuntimeError:
                print ('WARNING: Unable to read %s as a netCDF restart. '
                       'Attempting to read as a formatted AMBER7 restart...'
                       %self.filename)
                is_binary = False

        if not is_binary:
            # Use the rst7 object from parmed in AmberTools.
            tmp = _rst7_parmed(inpcrd_name)
            self.valid = tmp.valid
            self.title = tmp.title
            self.time = tmp.time
            self.coords = tmp.coords
            self.natom = tmp.natom
            if tmp.hasbox:
                self.hasbox = tmp.hasbox
                self.box = tmp.box
            if tmp.hasvels:
                self.hasvels = tmp.hasvels
                self.vels = tmp.vels
