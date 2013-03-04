from amberrun import AmberRun 
import os,shutil

if os.path.exists('snglptn'): shutil.rmtree('snglptn')
foo = AmberRun('-O',engine='sander',inpcrd='DMP_US_1.rst7',ref='DMP_US_1.rst7')

print "energy = %f kcal/mol"%foo.Energy()
print "       = %f kT"%foo.ReducedEnergy()
