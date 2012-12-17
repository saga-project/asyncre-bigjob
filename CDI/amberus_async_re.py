import sys, time, random, math
from pj_async_re import async_re_job
from amber_async_re import pj_amber_job

BOLTZMANN_CONSTANT = 1.3806*6.022/4148

class amberus_async_re_job(pj_amber_job,async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        #make sure AMBER umbrella sampling is wanted
        if self.keywords.get('RE_TYPE') != 'AMBERUS':
            self._exit("RE_TYPE is not AMBERUS")
        #AMBERUS runs with Amber
        if self.keywords.get('ENGINE') != 'AMBER':
            self._exit("ENGINE is not AMBER")
        #input files
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if not (self.extfiles is None):
            if self.extfiles != '':
                self.extfiles = self.extfiles.split(',')
        #list of force constants
        if self.keywords.get('FORCE_CONSTANTS') is None:
            self._exit("FORCE_CONSTANTS needs to be specified")
        kbiasline = self.keywords.get('FORCE_CONSTANTS')
        if [kbiasline] == kbiasline.split(':'):
            self.kbias = [[item] for item in kbiasline.split(',')]
        else:
            self.kbias = [item.split(',') for item in kbiasline.split(':')]
        self.nreplicas = len(self.kbias)
        #conversion for angle force constants (AMBER uses kcal/mol-rad^2!)
        nbias = len(self.kbias[0])
        if self.keywords.get('BIAS_IS_ANGLE') is None:
            self.bias_is_angle = [ False for n in range(nbias) ]
        else:
            self.bias_is_angle = self.keywords.get('BIAS_IS_ANGLE').split(',')
            if len(self.bias_is_angle) != nbias:
                self._exit("# of biases and BIAS_IS_ANGLE flags don't match")
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

        restraint_template = "%s.RST" % basename
        restraint_file = "r%d/%s_%d.RST" % (replica, basename, cycle)
        # read template buffer
        tfile = open(restraint_template, "r")
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        for i in range(len(self.kbias[:][0])):
            rk = self.kbias[stateid][i]
            r0 = self.posbias[stateid][i]
            tbuffer = tbuffer.replace('@rk%d@'%i,rk)
            tbuffer = tbuffer.replace('@r0%d@'%i,r0)
        # write out
        ofile = open(restraint_file, "w")
        ofile.write(tbuffer)
        ofile.close()

        input_template = "%s.inp" % basename
        input_file = "r%d/%s_%d.inp" % (replica, basename, cycle)
        # read template buffer
        tfile = open(input_template, "r")
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace("@n@",str(cycle))
        # write out
        ofile = open(input_file, "w")
        ofile.write(tbuffer)
        ofile.close()
      
    @staticmethod
    def bias_energy(bias_coords, force_constants, bias_positions, isAngle=None):
        """Calculate the (harmonic) bias energy of coordinates in a given state.
        """
        if isAngle == None: isAngle = [ False for i in range(len(bias_coords)) ]

        dr2 = [ (float(r) - float(r0))**2 
                for r,r0 in zip(bias_coords,bias_positions) ]
        uBias = 0.
        for i in range(len(dr2)):
            if isAngle[i]:
                uBias += float(force_constants[i])*(math.pi/180.)**2*dr2[i]
            else:     
                uBias += float(force_constants[i])*dr2[i]
        return uBias

    def _doExchange_pair(self,repl_a,repl_b):
        """Perform exchange of bias parameters.        
        """
        cycle_a = self.status[repl_a]['cycle_current']
        sid_a = self.status[repl_a]['stateid_current']
        rk_a = self.kbias[sid_a]
        r0_a = self.posbias[sid_a]
        bias_coord_a = self._extractLastRCs(repl_a,cycle_a)

        cycle_b = self.status[repl_b]['cycle_current'] 
        sid_b = self.status[repl_b]['stateid_current']
        rk_b = self.kbias[sid_b]
        r0_b = self.posbias[sid_b]
        bias_coord_b = self._extractLastRCs(repl_b,cycle_b)

        isAngle = self.bias_is_angle
        u_aa = amberus_async_re_job.bias_energy(bias_coord_a,rk_a,r0_a,isAngle)
        u_ab = amberus_async_re_job.bias_energy(bias_coord_b,rk_a,r0_a,isAngle)
        u_ba = amberus_async_re_job.bias_energy(bias_coord_a,rk_b,r0_b,isAngle)
        u_bb = amberus_async_re_job.bias_energy(bias_coord_b,rk_b,r0_b,isAngle)
        delta = (u_ab + u_ba) - (u_aa + u_bb)
        
        if self.keywords.get('VERBOSE') == "yes":
            print 'Pair Info:'
            
            print 'replica = %d'%repl_a
            print ' bias coordinate : bias position'
            for r,r0 in zip(bias_coord_a,r0_a): print ' %15s : %13s'%(r,r0)
                
            print 'replica = %d'%repl_b
            print ' bias coordinate : bias position'
            for r,r0 in zip(bias_coord_b,r0_b): print ' %15s : %13s'%(r,r0)

            print "delta = %f kcal/mol"%delta

        csi = random.random()
        P_ab = math.exp(-self.beta*delta)
        if P_ab > csi:
            self.status[repl_a]['stateid_current'] = sid_b
            self.status[repl_b]['stateid_current'] = sid_a

            if self.keywords.get('VERBOSE') == "yes":
                print "Accepted %f %f" % (P_ab,csi)
                print (self.status[repl_a]['stateid_current'], 
                       self.status[repl_b]['stateid_current'])
        else:
            if self.keywords.get('VERBOSE') == "yes":
                print "Rejected %f %f" % (P_ab,csi)

    def _extractLastRCs(self,repl,cycle):
        """Extracts the last set of reaction coordinates from NMRopt output
        """
        trace_file = "r%s/%s_%d.TRACE" % (repl,self.basename,cycle)
        data = self._getAmberUSData(trace_file)
        return data[-1][1:]

if __name__ == '__main__':

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
