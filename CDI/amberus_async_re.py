import sys, time, random, math
from pj_async_re import async_re_job
from amber_async_re import pj_amber_job


class amberus_async_re_job(pj_amber_job,async_re_job):

    def _checkInput(self):
        async_re_job._checkInput(self)
        #make sure BEDAM is wanted
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
        self.kbias = self.keywords.get('FORCE_CONSTANTS').split(',')
        self.nreplicas = len(self.kbias)
        #list of bias positions
        if self.keywords.get('BIAS_POSITIONS') is None:
            self._exit("BIAS_POSITIONS needs to be specified")
        self.posbias = self.keywords.get('BIAS_POSITIONS').split(',')
        if len(self.posbias) != self.nreplicas:
            self._exit("Number of FORCE_CONSTANTS not equal to number of BIAS_POSITIONS")
        #simulation temperature
        if self.keywords.get('TEMPERATURE') is None:
            self._exit("TEMPERATURE is a required parameter")
        temperature = float(self.keywords.get('TEMPERATURE'))
        self.beta = 1./(0.0019872041*temperature)
#        #build parameters for the lambda states
#        self._buildBEDAMStates(lambdas)

#    def _buildBEDAMStates(self,lambdas):
#        self.stateparams = []
#        for lambd in lambdas:
#            st = {}
#            st['lambda'] = lambd
#            self.stateparams.append(st)
#        return len(self.stateparams)

    def _buildInpFile(self, replica):
        """
Builds input file for a BEDAM replica based on template input file
BASENAME.inp for the specified replica at lambda=lambda[stateid] for the
specified cycle.
"""
        basename = self.basename
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']

        restraint_template = "%s.RST" % basename
        restraint_file = "r%d/%s_%d.RST" % (replica, basename, cycle)
        rk = self.kbias[stateid]
        r0 = self.posbias[stateid]
        # read template buffer
        tfile = open(restraint_template, "r")
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace("@rk@",rk)
        tbuffer = tbuffer.replace("@r0@",r0)
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

        

        #parse more input files here ...


        # update the history status file
#        ofile = open("r%d/state.history" % replica, "a")
#        ofile.write("%d %d %s\n" % (cycle, stateid, lambd))
#        ofile.close()
        

    def _doExchange_pair(self,repl_a,repl_b):
        """
Performs exchange of lambdas for BEDAM replica exchange.        
"""
        cycle_a = self.status[repl_a]['cycle_current']
        sid_a = self.status[repl_a]['stateid_current']
        rk_a = self.kbias[sid_a]
        r0_a = self.posbias[sid_a]
        angle_a = self._extractLast_Angles(repl_a,cycle_a)

        cycle_b = self.status[repl_b]['cycle_current'] 
        sid_b = self.status[repl_b]['stateid_current']
        rk_b = self.kbias[sid_b]
        r0_b = self.posbias[sid_b]
        angle_b = self._extractLast_Angles(repl_b,cycle_b)

#        u_aa = self.bias_energy(angle_a, rk_a, r0_a)
#        u_ab = self.bias_energy(angle_b, rk_a, r0_a)
#        u_ba = self.bias_energy(angle_a, rk_b, r0_b)
#        u_bb = self.bias_energy(angle_b, rk_b, r0_b)
#        delta = (u_ab + u_ba) - (u_aa + u_bb)
        delta = 100.0
        
        if self.keywords.get('VERBOSE') == "yes":
            print "Pair Info:"
            print "replica = %d angle = %s  r0 = %s" % (repl_a, angle_a, r0_a)
            print "replica = %d angle = %s  r0 = %s" % (repl_b, angle_b, r0_b)

        csi = random.random()
        if math.exp(-self.beta*delta) > csi:
            if self.keywords.get('VERBOSE') == "yes":
                print "Accepted %f %f" % (math.exp(-self.beta*delta),csi)
                print (self.status[repl_a]['stateid_current'], self.status[repl_b]['stateid_current'])
            self.status[repl_a]['stateid_current'] = sid_b
            self.status[repl_b]['stateid_current'] = sid_a
            if self.keywords.get('VERBOSE') == "yes":
                print (self.status[repl_a]['stateid_current'], self.status[repl_b]['stateid_current'])
        else:
            if self.keywords.get('VERBOSE') == "yes":
                print "Rejected %f %f" % (math.exp(-self.beta*delta),csi)

    def _extractLast_Angles(self,repl,cycle):
        """
Extracts binding energy from Impact output
"""
        trace_file = "r%s/%s_%d.TRACE" % (repl,self.basename,cycle)
        datai = self._getAmberUSData(trace_file)
        print datai
        nf = len(datai[0])
        nr = len(datai)
        print datai[-1][1]
        return datai[-1][1]



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
