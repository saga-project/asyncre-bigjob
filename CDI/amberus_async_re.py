import os, sys, time, random, math
from pj_async_re import async_re_job
from amber_async_re import pj_amber_job
# AMBER plugins
from amberio.AmberRestraint import ReadAmberRestraintFile
from chemistry.amber.readparm import rst7

BOLTZMANN_CONSTANT = 1.38*6.022/4184

class amberus_async_re_job(pj_amber_job,async_re_job):

    def _checkInput(self):
        pj_amber_job._checkInput(self)
        #make sure AMBER umbrella sampling is wanted
        if self.keywords.get('RE_TYPE') != 'AMBERUS':
            self._exit("RE_TYPE is not AMBERUS")
        #input files
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if not (self.extfiles is None):
            if self.extfiles != '':
                self.extfiles = self.extfiles.split(',')
        #flag for turning off exchange
        if ( self.keywords.get('DO_EXCHANGES') == 'False' 
             or self.keywords.get('DO_EXCHANGES') == 'No'):
            self.do_exchanges = False
        else:
            self.do_exchanges = True
        #list of force constants
        if self.keywords.get('FORCE_CONSTANTS') is None:
            self._exit("FORCE_CONSTANTS needs to be specified")
        kbiasline = self.keywords.get('FORCE_CONSTANTS')
        if [kbiasline] == kbiasline.split(':'):
            self.kbias = [[item] for item in kbiasline.split(',')]
        else:
            self.kbias = [item.split(',') for item in kbiasline.split(':')]
        self.nreplicas = len(self.kbias)

        #list of bias positions
        if self.keywords.get('BIAS_POSITIONS') is None:
            self._exit("BIAS_POSITIONS needs to be specified")
        posbiasline = self.keywords.get('BIAS_POSITIONS')
        if [posbiasline] == posbiasline.split(':'):
            self.posbias = [[item] for item in posbiasline.split(',')]
        else:
            self.posbias = [item.split(',') for item in posbiasline.split(':')] 
        if len(self.posbias) != self.nreplicas:
            msg = ('Number of FORCE_CONSTANTS not equal to number of'
                   ' BIAS_POSITIONS')
            self._exit(msg)

        #simulation temperature
        if self.keywords.get('TEMPERATURE') is None:
            self._exit("TEMPERATURE is a required parameter")
        temperature = float(self.keywords.get('TEMPERATURE'))
        self.beta = 1./(BOLTZMANN_CONSTANT*temperature)

    def _buildInpFile(self, replica):
        """
        Builds input file for a AMBER umbrella sampling replica based on a
        template input file, BASENAME.inp, for the specified replica and cycle.
        """
        basename = self.basename
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']
        restraint_template = '%s.RST'%basename
        restraint_file = 'r%d/%s_%d.RST'%(replica,basename,cycle)
        
        # 1) Make an AmberRestraint object from the template
        # 2) Modify it to match the current state
        # 3) Write a new restraint file for the replica
        umbrella_potential = ReadAmberRestraintFile(restraint_template)
        for m in range(len(self.kbias[:][0])): # iterate bias dimensions
            k  = float(self.kbias[stateid][m])
            r0 = float(self.posbias[stateid][m])
            umbrella_potential[m].rk[0] = k
            umbrella_potential[m].rk[1] = k
            umbrella_potential[m].r[0] = r0 - 100.
            umbrella_potential[m].r[1] = r0
            umbrella_potential[m].r[2] = r0
            umbrella_potential[m].r[3] = r0 + 100.
        umbrella_potential.WriteAmberRestraintFile(restraint_file)

        input_template = '%s.inp'%basename
        input_file = 'r%d/%s_%d.inp'%(replica,basename,cycle)
        # read template buffer
        tfile = open(input_template,'r')
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace('@n@',str(cycle))
        # write out
        ofile = open(input_file,'w')
        ofile.write(tbuffer)
        ofile.close()
      
    def _doExchange_pair(self,repl_a,repl_b):
        """Perform exchange of bias parameters.        
        """
        if self.do_exchanges:
            # cycle count and state id of the replicas
            cycle_a = self.status[repl_a]['cycle_current']
            #sid_a = self.status[repl_a]['stateid_current']
            cycle_b = self.status[repl_b]['cycle_current'] 
            #sid_b = self.status[repl_b]['stateid_current']

            # final coordinates from the last cycle of each replica 
            rst_file_a = 'r%d/%s_%d.rst7'%(repl_a,self.basename,cycle_a)
            crds_a = rst7(rst_file_a).coords
            rst_file_b = 'r%d/%s_%d.rst7'%(repl_b,self.basename,cycle_b)
            crds_b = rst7(rst_file_b).coords

            # current biasing potential of each replica
            RST_a = 'r%d/%s_%d.RST'%(repl_a,self.basename,cycle_a)
            umbrella_a = ReadAmberRestraintFile(RST_a)
            RST_b = 'r%d/%s_%d.RST'%(repl_b,self.basename,cycle_b)
            umbrella_b = ReadAmberRestraintFile(RST_b)

            # do the energy evaluations
            u_aa = umbrella_a.Energy(crds_a)
            u_ab = umbrella_a.Energy(crds_b)
            u_ba = umbrella_b.Energy(crds_a)
            u_bb = umbrella_b.Energy(crds_b)
            delta = (u_ab + u_ba) - (u_aa + u_bb)
            u = self.beta*delta
       
            Exchange = True
            P_ab = 1.
            csi = 0.
            if u > 0.:
                csi = random.random()
                P_ab = math.exp(-u)
                if csi > P_ab: Exchange = False
            if Exchange:
                sid_a = self.status[repl_a]['stateid_current']
                sid_b = self.status[repl_b]['stateid_current']
                self.status[repl_a]['stateid_current'] = sid_b
                self.status[repl_b]['stateid_current'] = sid_a

            if self.keywords.get('VERBOSE') == 'yes':
                # extract the actual coordinates for reporting purposes
                print ('======================================================='
                       '=========================')
                print 'Exchange Attempt : Replicas (%d,%d)'%(repl_a,repl_b)
                print 'Replica %d Coordinate/Umbrella Info:'%repl_a
                umbrella_a.PrintRestraintReport(crds_a)
                umbrella_a.PrintRestraintEnergyReport(crds_a)
                print 'Replica %d Coordinate/Umbrella Info:'%repl_b
                umbrella_b.PrintRestraintReport(crds_b)
                umbrella_b.PrintRestraintEnergyReport(crds_b)
                print ('======================================================='
                       '=========================')
                print 'U_a(x_b) - U_a(x_a) + U_b(x_a) - U_b(x_b) = %f kcal/mol'%delta
                if Exchange:
                    print 'Accepted! P(a<->b) = %f >= %f'%(P_ab,csi)
                    print ('New states for replicas (%s,%s) are (%s,%s)'
                           %(repl_a,repl_b,
                             self.status[repl_a]['stateid_current'], 
                             self.status[repl_b]['stateid_current']))
                else:
                    print 'Rejected! P(a<->b) = %f < %f'%(P_ab,csi)
                print ('======================================================='
                       '=========================')
                
if __name__ == '__main__':

    BIGJOB_VERBOSE=100

    # Parse arguments:
    usage = "%prog <ConfigFile>"
    
    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)
    
    commandFile = sys.argv[1]

    print ""
    print "===================================="
    print "AMBER Asynchronous Replica Exchange "
    print "===================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    rx = amberus_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
