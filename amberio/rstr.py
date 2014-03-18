import re

from numpy import pi, asarray, inf

from namelist import Namelist, NamelistCollection
from coordinates import Bond, Angle, Dihedral

               
class AmberRestraint(NamelistCollection):
    def __init__(self, title='', *rsts):
        NamelistCollection.__init__(self,*rsts)
        self.title = str(title).rstrip()

    @classmethod
    def from_disang(cls, filename):
        """Create an AmberRestraint from a file with &rst namelists."""
        nls,lines = NamelistCollection.separate_nls(filename)
        rst = AmberRestraint(title=''.join(lines))
        for nl in nls.matches('rst'):
            try:
                iat = [int(i) for i in nl.pop('iat').split(',')]
            except KeyError:
                raise KeyError("'iat' is required for nmropt restraints!")
            try:
                rstwt = [float(w) for w in nl.pop('rstwt').split(',')]
                rst.append(GenDistCoordRestraint(iat=iat,rstwt=rstwt,**nl))
            except KeyError:
                if len(iat) == 2:
                    rst.append(BondRestraint(iat=iat,**nl))
                elif len(iat) == 3:
                    rst.append(AngleRestraint(iat=iat,**nl))
                elif len(iat) == 4:
                    rst.append(TorsionRestraint(iat=iat,**nl))
                else:
                    raise ValueError("Bad 'iat' specification in %s"%filename)
        if len(rst) < 1:
             print 'WARNING! No &rst namelists found in %s.'%filename
        return rst

    @classmethod
    def from_mdout(cls, filename, rstwt=None):
        params = _extract_restraint_params_from_mdout(filename)
        rst = AmberRestraint()
        for i,iat in enumerate(params['iat']):
            #    Since the restraint report in mdout files does not contain 
            # rstwt data, the only way to identify generalized distance 
            # coordinates is to provide these a priori. If the 'rstwt' argument
            # is not 'None', it is expected to be a sequence of rstwts. The 
            # elements in rstwt are then assigned to the restraints in order 
            # until the sequence is exhausted.
            try:
                rstwt[i][:]
                rst.append(GenDistCoordRestraint(iat=iat,rstwt=rstwt[i]))
                is_gendistcoord = True
            except (TypeError,IndexError):
                is_gendistcoord = False
            if not is_gendistcoord:
                if len(iat) == 2:
                    rst.append(BondRestraint(iat=iat))
                elif len(iat) == 3:
                    rst.append(AngleRestraint(iat=iat))
                elif len(iat) == 4:
                    rst.append(TorsionRestraint(iat=iat))
                else:
                    raise Exception('Bad iat specification')
            for param in ['R1','R2','R3','R4','RK2','RK3']:
                rst[i][param.lower()] = params[param][i]
        return rst

    @property
    def r1(self):
        return asarray([rst['r1'] for rst in self])

    @property
    def r2(self):
        return asarray([rst['r2'] for rst in self])

    @property
    def r3(self):
        return asarray([rst['r3'] for rst in self])

    @property
    def r4(self):
        return asarray([rst['r4'] for rst in self])

    @property
    def rk2(self):
        return asarray([rst['rk2'] for rst in self])

    @property
    def rk3(self):
        return asarray([rst['rk3'] for rst in self])

    @property
    def image_dist(self):
        return asarray([rst.image_dist for rst in self])

    @property
    def kfac(self):
        return asarray([rst.kfac for rst in self])

    def __setattr__(self, name, value):
        if name in ['r0','r1','r2','r3','r4','k0','rk2','rk3']:
            try:
                for rst,v in zip(self,value):
                    rst[name] = float(v)
            except TypeError:
                for rst in self:
                    rst[name] = float(value)
        else:
            NamelistCollection.__setattr__(self,name,value)

    def write(self, outfile, mode='w'):
        """Write an AMBER restraint file with the current restraints."""
        closeafter = False
        if hasattr(outfile,'write'):
            pass
        else:
            closeafter = True
            outfile = open(outfile,mode)
        outfile.write('%s\n%s'%(self.title,str(self)))
        if closeafter: 
            outfile.close()

    def coordinates(self, x, n=None):
        """Return the restraint coordinate values from a 3N coordinate list."""
        return [rst.coordinate(x) for rst in self[:n]]


class NmroptRestraint(Namelist):
    def __init__(self, *args, **kwargs):
        Namelist.__init__(self,'rst',' ','=',' ',72,72,*args,**kwargs)

    def __setitem__(self, key, value):
        if key == 'k0':
            # convenience variable for harmonic potential
            self['rk2'] = float(value)
            self['rk3'] = float(value)
        elif key == 'r0':
            # convenience variable for harmonic potential
            self['r2'] = float(value)
            self['r3'] = float(value)
            self['r1'] = self._r1_harmonic
            self['r4'] = self._r4_harmonic
        elif key in ['rk2','rk3']:
            # some restraints use internal units (e.g. angles)
            Namelist.__setitem__(self,key,float(value)) 
        else:
            Namelist.__setitem__(self,key,value)

    def __getitem__(self, key):
        try:
            return Namelist.__getitem__(self,key)
        except KeyError:
            self[key] = value = self.default(key)
            return value

    @property
    def image_dist(self):
        return inf

    @property
    def kfac(self):
        return 1.0

    def __str__(self):
        # Convert list values to strings.
        self['iat'] = ','.join([str(i) for i in self['iat']])
        try:
            self['rstwt'] = ','.join([str(i) for i in self['rstwt']])
        except KeyError:
            pass
        txt = Namelist.__str__(self)
        # Convert back to lists.
        self['iat'] = [int(i) for i in self['iat'].split(',')]
        try:
            self['rstwt'] = [float(i) for i in self['rstwt'].split(',')]
        except KeyError:
            pass
        return txt

    def default(self, key):
        if key in ['r2','r3','rk2','rk3']:
            self[key] = 0.0
            return 0.0
        elif key == 'r1':
            self[key] = self._r1_harmonic
            return self._r1_harmonic 
        elif key == 'r4':
            self[key] = self._r4_harmonic
            return self._r4_harmonic 
        else:
            raise KeyError('No default for key %s!'%key)

class BondRestraint(NmroptRestraint):
    """A restraint on the "bond" between two atoms."""
    def coordinate(self, x):
        i,j = self['iat'][0]-1,self['iat'][1]-1
        return Bond(x,i,j)

    @property
    def _r1_harmonic(self):
        return 0.0
    
    @property
    def _r4_harmonic(self):
        return self['rk2'] + 500.0
    

class AngleRestraint(NmroptRestraint):
    """A restraint on the angle between three atoms."""
    def coordinate(self, x):
        i,j,k = self['iat'][0]-1,self['iat'][1]-1,self['iat'][2]-1
        return Angle(x,i,j,k)*(180/pi) # in degrees

    @property
    def kfac(self):
        return (pi/180)**2

    @property
    def _r1_harmonic(self):
        return 0.0

    @property
    def _r4_harmonic(self):
        return 180.0


class TorsionRestraint(NmroptRestraint):
    """A restraint on the dihedral defined by four atoms."""
    def coordinate(self, x):
        i,j = self['iat'][0]-1,self['iat'][1]-1,
        k,l = self['iat'][2]-1,self['iat'][3]-1
        return Dihedral(x,i,j,k,l)*(180/pi) # in degrees

    @property
    def kfac(self):
        return (pi/180)**2

    @property
    def image_dist(self):
        return 360

    @property
    def _r1_harmonic(self):
        return self['r2'] - 180.0

    @property
    def _r4_harmonic(self):
        return self['r2'] + 180.0


class GenDistCoordRestraint(NmroptRestraint):
    """A restraint on a linear combination of distances."""
    def coordinate(self, x):
        i_s = [i-1 for i in self['iat'][0::2]]
        j_s = [j-1 for j in self['iat'][1::2]]
        return sum([w*Bond(x,i,j) for i,j,w in zip(i_s,j_s,self['rstwt'])])
            
    @property
    def _r1_harmonic(self):
        return 0.0

    @property
    def _r4_harmonic(self):
        return self['r2'] + 500.0

def _extract_restraint_params_from_mdout(mdout_name):
    """Extract the restraint parameters from the report in an mdout file."""
    iat1_pattern = "(\([ 0-9]+?\))" # pattern for atom 1
    iat1p_pattern = "(\([ 0-9]+?\))*?" # pattern for atoms 1+
    iat_pattern = ("^.*?%s.*?%s[^\)\(]*?%s[^\)\(]*?%s[^\)\(]*$"
                   %(iat1_pattern,iat1p_pattern,iat1p_pattern,iat1p_pattern))
    iat_pattern_obj = re.compile(r"%s"%iat_pattern)

    params = {'iat': [], 'R1': [], 'R2': [], 'R3': [], 'R4': [],
              'RK2': [], 'RK3': []}
    pattern = ("^R1 =(?P<R1>[ 0-9\-\+\.]+) "
               "+R2 =(?P<R2>[ 0-9\-\+\.]+) "
               "+R3 =(?P<R3>[ 0-9\-\+\.]+) "
               "+R4 =(?P<R4>[ 0-9\-\+\.]+) "
               "+RK2 =(?P<RK2>[ 0-9\-\+\.]+) "
               "+RK3 =(?P<RK3>[ 0-9\-\+\.]+)")
    pattern_obj = re.compile(r'%s'%pattern)

    at_toadd = []
    for line in open(mdout_name,'r'):
        # Accumulate a list of integers between parantheses. These may or may
        # not correspond to the atom selection of a restraint.
        iat_match = iat_pattern_obj.search(line)
        if iat_match is not None:
            for iat_str in iat_match.groups():
                if iat_str is not None:
                    at_toadd.append(int(iat_str.strip(')').strip('(').strip()))
        # Look for a restraint report. If one is found, the atom list DOES
        # belong to a restraint. Append the atom list.
        param_match = pattern_obj.search(line)
        if param_match is not None:
            params['iat'].append(at_toadd)
            at_toadd = []
            for param in params.keys():
                if param != 'iat':
                    params[param].append(float(param_match.group(param)))
    return params

if __name__ == '__main__':
    import sys
    test = AmberRestraint.from_disang(sys.argv[1])
    test.write(sys.stdout)
