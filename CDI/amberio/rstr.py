"""                                                                             
FILE: rstr.py - plugin for AMBER nmropt style restraint files

DESCRIPTION: This module provides a convenient implementation of restraints
in AMBER for use in Python. This is useful, for example, if one wants to
determine just the restraint energy of a set of coordinates without
running a full energy evaluation in AMBER.

AUTHOR: Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>

REFERENCES: AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
"""
from math import pi
import sys

__all__ = ['ReadAmberRestraintFile','AmberRestraint','NmroptRestraint']

def ReadAmberRestraintFile(rstFilename):
    """
    Read an AMBER restraint file (usually file extension RST) and return an
    AmberRestraint object containing each restraint found therein.

    REQUIRED ARGUMENTS:
    rstFilename - file containing AMBER "&rst" namelists

    RETURN VALUES:
    amberRst - AmberRestraint object (list of NmroptRestraint)
    """
    from namelist import ReadNamelists
    amberRstr = AmberRestraint()
    # Read the namelists from the restraint file and look for &rst namelists
    namelists = ReadNamelists(rstFilename)
    for nl in namelists:
        if nl.name == 'rst':
            # only the iat keyword is required to define a restraint
            if not nl.has_key('iat'):
                msg = "'iat' must be specified to define an nmropt restraint."
                raise Exception(msg)
            amberRstr.append( NmroptRestraint(nl.pop('iat'),**nl) )
    if len(amberRstr) < 1:
        raise Exception('No &rst namelists found in file: %s'%rstFilename)
    return amberRstr

class AmberRestraint(list):
    """
    A collection (list) of NmroptRestraints defining an AMBER restraint energy.

    (See NmroptRestraint for further details.)
    """
    def __init__(self,*rstrs):
        list.__init__(self)
        for rstr in rstrs: self.extend(rstr)
        
    def append(self, item):
        if not isinstance(item, NmroptRestraint):
            raise TypeError('AmberRestraints must contain NmroptRestraints!')
        list.append(self, item)
    
    def extend(self, items):
        if hasattr(items, '__iter__'):
            for item in items:
                self.append(item)
            
    def SetRestraintParameters(self,**rstr_params):
        """
        Convenience function for setting any of r0, r1, r2, r3, r4, k0, rk2, and
        rk3 by list assignment. 

        Parameters will be assigned in order until either no parameters or 
        restraints are left. Thus, paramaters can be set for restraints 1 and 2,
        but not 1 and 3. For the latter case the individual restraints must
        be modified directly.
        """
        nrstrs = len(self)
        for params,values in rstr_params.iteritems():
            nvalues = len(values)
            # case 1: more restraints than parameter values
            if nrstrs >= nvalues:
                for i,value in enumerate(values):
                    self[i].SetRestraintParameters(**{params:float(value)})
            # case 1: more parameter values than restraints 
            #         (better to ask for forgiveness than permission)
            else:
                print ('Warning: %d restraint parameters were specified, but'
                       ' only %s restraints are defined.'%(nvalues,nrstrs))
                for i,rstr in enumerate(self):
                    value = float(rstr_params[params][i])
                    rstr.SetRestraintParameters(**{params:value})

    def Energy(self,crds):
        """Calculate the total energy from all restraints.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        energy - total restraint energy (in kcal/mol)
        """
        energy = 0.
        for rst in self: energy += rst.Energy(crds)
        return energy

    def EnergyAndGradients(self,crds):
        """
        THIS IS CURRENTLY BROKEN AND RETURNS ALL GRADIENTS AS ZERO!

        Calculate the total energy and gradients from all restraints.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        energy - restraint energy (in kcal/mol)
        gradients - 3N list of cartesian gradients (in kcal/mol-x)

        x - either Angstroms or radians
        """
        energy = 0.
        gradients = [ 0. for n in range(len(crds)) ]
        for rst in self:
            e,g = rst.EnergyAndGradients(crds)
            energy += e
            for i in range(len(g)): gradients[i] += g[i]
        return energy,gradients

    def DecomposeEnergy(self,crds):
        """Decompose (by restraint type) the total energy from all restraints.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        energy - dict containing decomposed energies (in kcal/mol); the 
        'Restraint' key corresponds to the total energy
        """
        energy = { 'Restraint':0., 'Bond':0., 'Angle':0., 'Torsion':0.,
                  'Gen. Dist. Coord.':0. }
        for rst in self: energy[rst.rstType] += rst.Energy(crds)
        eRestraint = 0.
        for key in energy.keys(): eRestraint += energy[key]
        energy['Restraint'] = eRestraint
        return energy

    def PrintRestraintReport(self,crds=None,anames=None):
        """
        Print a restraint report in the same format as AMBER output.
        Information not available from a restraint file (e.g. atom names) will
        be omitted unless additional information is provided.

        OPTIONAL ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)
        anames - N list of atom names
        """
        for rst in self:
            print '******' 
            rst.PrintRestraintReport(crds,anames)
        print '%23sNumber of restraints read = %5d'%('',len(self))

    def PrintRestraintEnergyReport(self,crds):
        """Print a restraint energy report in the same format as AMBER output.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)
        """
        e = self.DecomposeEnergy(crds)
        print (' NMR restraints: Bond =%9.3f   Angle = %9.3f   Torsion = %9.3f'
               %(e['Bond'],e['Angle'],e['Torsion']))
        if e['Gen. Dist. Coord.'] > 0.:
            print ('               : Gen. Dist. Coord. = %9.3f'
                   %e['Gen. Dist. Coord.'])

    def WriteAmberRestraintFile(self,outfile,title=''):
        """ Write an AMBER restraint file with all of the current restraints.

        REQUIRED ARGUMENTS:
        outfile - file object or name for writing
        """
        if hasattr(outfile,'write'):
            pass
        elif isinstance(outfile,str):
            outfile = open(outfile,'w')
        else:
            raise TypeError("'outfile' must be either a string or file object.")
        outfile.write('%s\n'%title)
        for rst in self: rst.WriteRstNamelist(outfile,closeafter=False)
        if outfile is not sys.__stdout__: outfile.close()
      
class NmroptRestraint(object):
    """Create a restraint object like that in the AMBER nmropt module. The
    restraint form is a harmonic flat-bottom well with six parameters: four
    positions (r1-r4) and two force constants (rk2 and rk3).

    From the AMBER 12 Manual (section 6.1.1 p. 204):
    ---
    the restraint is a well with a square bottom with parabolic sides out to a 
    defined distance, and then linear sides beyond that. If R is the value of 
    the restraint in question:

    - R < r1 Linear, with the slope of the "left-hand" parabola at the point R=r1
    - r1 <= R < r2 Parabolic, with restraint energy k2(R-r2)^2
    - r2 <= R <r3 E = 0
    - r3 <= R <r4 Parabolic with restraint energy k3(R-r3)^2
    - r4 <= R Linear, with the slope of the "right-hand" parabola at the point R=r4
    ----

    Frequently one only desires a purely harmonic well, in which case r1 << r2,
    r2=r3, r3 << r4, and rk2 = rk3. These defaults can be obtained by only
    specifying r0 and k0 (as in sander from AMBER 10 onward).

    REQUIRED ARGUMENTS:
    iat - list of atom indices defining the restraint

    OPTIONAL ARGUMENTS:
    rstr_params - any of r0, r1, r2, r3, r4, k0, rk2, and rk3 can be set by
    direct assignment. r0 and k0 will override all other specifications.

    NB: As in AMBER, angle positions are in degrees while angle force constants
    are in radians. Distances are always in Angstroms.
    """
    def __init__(self,iat,**rstr_params):
        # Determine the restraint type from atom and rstwt info
        if isinstance(iat,(tuple,list)):
            self.iat = iat
        else:
            self.iat = tuple([ int(i) for i in iat.split(',') ])
        self.nAtoms = len(self.iat)
        if 'rstwt' in rstr_params.keys(): 
            self.rstwt = [ float(w) for w in rstr_params['rstwt'].split(',') ]
            self.rstType = 'Gen. Dist. Coord.'
            if len(self.rstwt) != self.nAtoms/2:
                msg = ('Wrong number of rstwt values provided for %d atom Gen.' 
                       ' Dist. Coord. Expected %d, but got %d.'
                       %(self.nAtoms,self.nAtoms/2,len(self.rstwt)))
                raise Exception(msg)
        else:                             
            self.rstwt = None
            if self.nAtoms == 2:   self.rstType = 'Bond'
            elif self.nAtoms == 3: self.rstType = 'Angle'
            elif self.nAtoms == 4: self.rstType = 'Torsion'
            else:
                msg = 'Invalid restraint specification (or not supported yet).'
                raise Exception(msg)
        # Set the restraint parameters
        self.r  = [0., 0., 0., 0.]
        self.rk = [0., 0.]
        self.SetRestraintParameters(**rstr_params)
        
    def SetRestraintParameters(self,**rstr_params):
        """Set any of r0, r1, r2, r3, r4, k0, rk2, and rk3 by assignment.
        """
        # A pure harmonic restraint can be set with just r0, 
        if rstr_params.has_key('r0'):
            r0 = float(rstr_params['r0'])
            self.r[1:3] = [r0,r0]
            if self.rstType == 'Torsion':
                 self.r[0] = r0 - 180.
                 self.r[3] = r0 + 180.
            elif self.rstType == 'Angle':
                self.r[0] = 0.
                self.r[3] = r0 + 180.
            else: # all distance restraints
                self.r[0] = r0 - 500.
                self.r[3] = r0 + 500.
            # Use radians for angle and torsion calculations
            if self.rstType == 'Angle' or self.rstType == 'Torsion':
                self.r = [ r*pi/180. for r in self.r ]
        # otherwise the four positions need to be set individually.
        else:
            rs = {'r1':0,'r2':1,'r3':2,'r4':3}
            for r,i in rs.iteritems():
                if rstr_params.has_key(r):
                    self.r[i] = float(rstr_params[r])
                    # Use radians for angle and torsion calculations
                    if self.rstType == 'Angle' or self.rstType == 'Torsion':
                        self.r[i] = self.r[i]*pi/180.

        # Check the relative restraint positions.
        if not self.r[0] <= self.r[1] <= self.r[2] <= self.r[3]:
            msg = ('Restraint positions must be monotonically increasing'
                   ' (r1 <= r2 <= r3 <= r4).')
            raise ValueError(msg)

        #  A pure harmonic restraint can be set with just k0,
        if rstr_params.has_key('k0'):
            k0 = float(rstr_params['k0'])
            self.rk = [k0,k0]
        # otherwise the two force constants need to be set individually.
        else:
            if rstr_params.has_key('rk2'):
                self.rk[0] = float(rstr_params['rk2'])
            if rstr_params.has_key('rk3'):
                self.rk[1] = float(rstr_params['rk3'])

    def Coord(self,crds):
        """Calculate the restraint coordinate.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        r - the restraint coordinate (in Angstroms or radians)
        """
        r,drdx = self.CoordAndGradients(crds)
        return r

    def CoordAndGradients(self,crds):
        """
        THIS IS CURRENTLY BROKEN AND RETURNS ALL GRADIENTS AS ZERO!

        Calculate the restraint coordinate and the gradient along the atom 
        cartesian coordinates (in au).

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        r - the restraint coordinate (in Angstroms or radians)
        drdx - 3N list of cartesian gradients (unitless or radians/Angstrom)
        """
        import coordinates
        r = 0.
        drdx = [ 0. for n in range(len(crds)) ]
        if self.rstType == 'Bond': 
            i = self.iat[0] - 1
            j = self.iat[1] - 1
            r,drdxi,drdxj = coordinates.BondAndGradients(crds,i,j)
            drdx[3*i:3*(i+1)] = drdxi
            drdx[3*j:3*(j+1)] = drdxj
        elif self.rstType == 'Angle':
            i = self.iat[0] - 1
            j = self.iat[1] - 1
            k = self.iat[2] - 1
            r,drdxi,drdxj,drdxk = coordinates.AngleAndGradients(crds,i,j,k)
            drdx[3*i:3*(i+1)] = drdxi
            drdx[3*j:3*(j+1)] = drdxj
            drdx[3*k:3*(k+1)] = drdxk
        elif self.rstType == 'Torsion':
            i = self.iat[0] - 1
            j = self.iat[1] - 1
            k = self.iat[2] - 1
            l = self.iat[3] - 1
            r,drdxi,drdxj,drdxk,drdxl = (
                coordinates.DihedralAndGradients(crds,i,j,k,l) )
            drdx[3*i:3*(i+1)] = drdxi
            drdx[3*j:3*(j+1)] = drdxj
            drdx[3*k:3*(k+1)] = drdxk
            drdx[3*l:3*(l+1)] = drdxl
            # Get the closest periodic image/phase for torsions.
            rmean = (self.r[1] + self.r[2])/2.
            isNearestImage = False
            while not isNearestImage:
                if r - rmean > pi:
                    r -= 2*pi
                elif rmean - r > pi:
                    r += 2*pi
                else:
                    isNearestImage = True
        elif self.rstType == 'Gen. Dist. Coord.':
            for k in range(len(self.rstwt)):
                i = self.iat[2*k+0] - 1
                j = self.iat[2*k+1] - 1
                x,drdxi,drdxj = coordinates.BondAndGradients(crds,i,j)
                r += self.rstwt[k]*x
                for m in range(3):
                    drdx[3*i+m] += self.rstwt[k]*drdxi[m]
                    drdx[3*j+m] += self.rstwt[k]*drdxj[m]
        else:
            msg = 'ERROR! Problem with restraint type (%s)'%self.rstType
            raise Exception(msg)
        return r,drdx

    def Energy(self,crds):
        """Calculate the restraint energy.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        energy - restraint energy (in kcal/mol)
        """
        energy,gradients = self.EnergyAndGradients(crds)
        return energy

    def EnergyAndGradients(self,crds):
        """
        THIS IS CURRENTLY BROKEN AND RETURNS ALL GRADIENTS AS ZERO!
        
        Calculate the restraint energy and gradients.

        REQUIRED ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)

        RETURN VALUES:
        energy - restraint energy (in kcal/mol)
        gradients - 3N list of cartesian gradients (in kcal/mol-x)

        x - either Angstroms or radians
        """
        r,drdx = self.CoordAndGradients(crds)
        dedr = 0.
        energy = 0.
        # The following is the harmonic, flat-bottomed well described in the
        # AMBER manual. See the class documentation for more details.
        if r < self.r[0]:
            dr = self.r[0] - self.r[1]
            dedr = 2.*self.rk[0]*dr 
            energy = dedr*(r-self.r[0]) + self.rk[0]*dr**2
        elif self.r[0] <= r < self.r[1]:
            dr = r - self.r[1]
            dedr = 2.*self.rk[0]*dr
            energy = self.rk[0]*dr**2
        elif self.r[1] <= r < self.r[2]:
            dedr = 0.
            energy = 0.
        elif self.r[2] <= r < self.r[3]:
            dr = r - self.r[2]
            dedr = 2.*self.rk[1]*dr
            energy = self.rk[1]*dr**2
        else:
            dr = self.r[3] - self.r[2]
            dedr = 2.*self.rk[1]*dr
            energy = dedr*(r-self.r[3]) + self.rk[1]*dr**2
        # Use the chain rule to get the gradient (dE/dx) along the atomic 
        # cartesian coordinates:
        #
        # dE/dx = (dE/dr)(dr/dx)
        #
        # where E is the energy, r is the restraint coordinate, and x is some 
        # atomic cartesian coordinate.
        gradients = [ dedr*drdxi for drdxi in drdx ]
        return energy,gradients

    def PrintRestraintReport(self,crds=None,anames=None):
        """
        Print a restraint report in the same format as AMBER output.
        Information not available from a restraint file (e.g. atom names) will
        be omitted unless additional information is provided.

        OPTIONAL ARGUMENTS:
        crds - 3N list of coordinates (in Angstroms)
        anames - N list of atom names
        """
        # Use blanks if atom names are not available
        if anames is None: anames = [ '' for i in range(self.nAtoms) ]
        else:              anames = [ anames[i-1] for i in self.iat ]
        # The format for restraints on 4 or fewer atoms is rather simple:
        if self.nAtoms <= 4:
            line = ' '
            for i in range(self.nAtoms):
                line += '%-4s(%5d)-'%(anames[i],self.iat[i])
            print line[:-1]
        # The format for more complicated restraints is a little muckier:
        else:
            line1 = ' '
            for i in range(4): 
                line1 += '%-4s(%5d)-'%(anames[i],self.iat[i])
            line2 = ' '
            for i in range(4,self.nAtoms): 
                line2 += '%-4s(%5d)-'%(anames[i],self.iat[i])
            print line1
            print line2
        # Report angles/torsions in degrees (force constants stay in radians)
        r = self.r
        if self.rstType == 'Angle' or self.rstType == 'Torsion':
            r = [ ri*180./pi for ri in self.r ]
        print ('R1 =%8.3f R2 =%8.3f R3 =%8.3f R4 =%8.3f RK2 =%8.3f RK3 = %8.3f'
               %(r[0],r[1],r[2],r[3],self.rk[0],self.rk[1]))
    
        # Current coordinate values can only be gotten from a crd file.
        if crds is not None: 
            Rcurr = self.Coord(crds)
            delta_avg = abs(Rcurr - (self.r[1]+self.r[2])/2.)
            min_delta = min(abs(Rcurr - self.r[1]),abs(Rcurr - self.r[2]))
            # Report angles and torsions in degrees
            if self.rstType == 'Angle' or self.rstType == 'Torsion':
                Rcurr *= 180./pi
                delta_avg *= 180./pi
                min_delta *= 180./pi
            print (' Rcurr: %8.3f  Rcurr-(R2+R3)/2: %8.3f  '
                   'MIN(Rcurr-R2,Rcurr-R3): %8.3f'%(Rcurr,delta_avg,min_delta))

    def WriteRstNamelist(self,outfile=sys.stdout,closeafter=True):
        """Write an &rst namelist with the current restraint information.

        OPTIONAL ARGUMENTS:
        outfile - file object or name to write data to, default=stdout
        closeafter (bool) - flag to close outfile after writing, default=True
        """
        if hasattr(outfile,'write'):
            pass
        elif isinstance(outfile,str):
            outfile = open(outfile,'w')
        else:
            raise TypeError("'outfile' must be either a string or file object.")
        iat = ''
        for i in self.iat: iat += '%d,'%i
        # Report angles/torsions in degrees (force constants stay in radians)
        r = self.r
        if self.rstType == 'Angle' or self.rstType == 'Torsion':
            r = [ ri*180./pi for ri in self.r ]
        outfile.write(' &rst iat=%s r1=%f r2=%f r3=%f r4=%f rk2=%f rk3=%f'
                      %(iat,r[0],r[1],r[2],r[3],self.rk[0],self.rk[1]))
        if self.rstwt is not None:
            rstwt = ''
            for w in self.rstwt: rstwt += '%f,'%w
            outfile.write(' rstwt=%s'%rstwt)
        outfile.write(' / \n')
        if closeafter and outfile is not sys.__stdout__: outfile.close()

if __name__ == '__main__':
    import sys
    argc = len(sys.argv)
    print '=== AmberRestraint Test Suite ==='
    if 5 > argc < 2:
        print 'usage: AmberRestraint.py RST [[inpcrd] prmtop]'
        print 
        print 'RST - AMBER restraint file'
        print 'inpcrd - AMBER crd file (optional test on coordinates)'
        print 'prmtop - AMBER parm file (optional test for reporting)'
        print
        sys.exit()

    # Read restraints from a file
    rstFile = sys.argv[1]
    print 'reading AMBER restraint file: %s'%rstFile
    print '>>> rstTest =  ReadAmberRestraintFile(%s)'%rstFile
    rstTest = ReadAmberRestraintFile(rstFile)

    # Create an nmropt restraint in pure python
    print 'testing creation of restraints in pure python:'
    print '>>> rstTest.append(NmroptRestraint((1,2)))'
    print '(makes a dummy restraint between atoms 1 and 2)'
    rstTest.append(NmroptRestraint((1,2)))
    print '>>> rstTest[-1].SetRestraintParameters(r0=1.0,k0=10.)'
    print '(sets a pure harmonic potential with only 2 parameters)'
    rstTest[-1].SetRestraintParameters(r0=1.0,k0=10.)

    # Read crd and parm files if present
    crds = None
    anames = None
    if argc > 2:
        from chemistry.amber.readparm import rst7
        crdFile = sys.argv[2]
        print 'Reading coordinates from file: %s'%crdFile
        crds = rst7(crdFile).coords
    if argc > 3:
        from chemistry.amber.readparm import AmberParm
        prmFile = sys.argv[3]
        print 'reading parm info from file: %s'%prmFile
        anames = AmberParm(prmFile).parm_data['ATOM_NAME']

    # Print a report of the restraints so far
    print 'printing an AMBER-style report:'
    print '>>> rstTest.PrintRestraintReport(crds,anames)'
    rstTest.PrintRestraintReport(crds,anames)
    
    if crds is not None:
        # Print a report like those found after MD steps
        print
        print 'printing an AMBER-style restraint energy decomposition:'
        print '>>> rstTest.PrintRestraintEnergyReport(crds)'
        rstTest.PrintRestraintEnergyReport(crds)

        # Test the total energy and forces against an AMBER forcedump.dat
        print
        print 'calculating total restraint energy and forces:'
        print '>>> energy,gradients = rstTest.EnergyAndGradients(crds)'
        e,g = rstTest.EnergyAndGradients(crds)
        print 'RESTRAINT  = %12.4f'%e
        print 'Forces (same format as forcedump.dat)'
        for i in range(0,len(g),3):
            print ' % 18.16e % 18.16e % 18.16e'%(-g[i+0],-g[i+1],-g[i+2])           
    print
    print 'writing a new restraint file to stdout:'
    print '>>> rstTest.WriteAmberRestraintFile(sys.stdout)'
    rstTest.WriteAmberRestraintFile(sys.stdout)
