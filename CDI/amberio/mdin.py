"""
FILE: mdin.py - plugin for I/O of AMBER mdin files

DESCRIPTION: This module creates objects from AMBER mdin files. The objects are
essentially collections of Fortran namelists that also have knowledge of the
default values expected by AMBER (depending on the MD engine being used).

AUTHOR: Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>

REFERENCES: AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
            MMPBSA.py source code: ambermd.org/AmberTools-get.html
"""
__all__ = ['ReadAmberMdinFile','AmberMdin']

from namelist import *

def ReadAmberMdinFile(mdin_name, engine='sander'):
    """Read an AMBER mdin file and return an AmberMdin object.
    """
    # Read mdin to find namelists
    namelists = ReadNamelists(mdin_name)
    # Read again to find title and NMR variables
    otherLines = ReadEverythingExceptNamelists(mdin_name)
    title = ''
    nmrVarsPresent = False
    for line in otherLines:
        tokens = line.split('=')
        if len(tokens) > 1:
            # If '=' is present, assign the rhs value to the lhs variable.
            variable = tokens[0].strip()
            value = tokens[1].strip()
            if not nmrVarsPresent:
                # Create a new namelist with name = None.
                nmrVarsPresent = True
                namelists.append(AmberNamelist(None,{variable:value}))
            else:
                # Add a new variable to the nmr (None) namelist.
                namelists.GetFirstMatch(None)[variable] = value
        else:
            # Otherwise, append the line to the title.
            title += line + '\n'
    return AmberMdin(title=title,engine=engine,*namelists)

class AmberMdin(object):
    """
    Provides easy access to the data in an AMBER mdin file. The main
    attribute is a NamelistCollection (a dict derived object called "namelists")
    of the currently defined namelists. For potentially non-obvious reasons, it
    is best to "get" and "set" the variables in each namelist with the 
    provided functions. Convenient helper functions are also available to
    perform complicated tasks such as adding restraints, etc.
    """
    def __init__(self,*args,**kwargs):
        # First, check for title and engine information
        self.title = ''
        if kwargs.has_key('title'): self.title = kwargs.pop('title')
        engine = 'sander'
        if kwargs.has_key('engine'): engine = kwargs.pop('engine')
        # Now, collect the namelists:
        self.namelists = NamelistCollection()
        # Append any Namelist or Namelist-derived objects to the list.
        for arg in args:
            if isinstance(arg,Namelist):
                self.namelists.append(AmberNamelist(arg.name,arg))
        # Check for Namelists specified as named dicts.
        for name,namelist in kwargs.iteritems():
            if isinstance(namelist,dict):
                self.namelists.append(AmberNamelist(name,namelist))
        # Set default values of the namelists for the given engine
        # NB: This must be done last since the namelist values need to be given
        # before the remaining defaults can be set.
        self.SetDefaults(engine)

    def SetDefaults(self, engine):
        """Set the default variables for a particular MD engine (e.g. sander).
        """
        for nl in self.namelists: nl.SetDefaults(engine)

    def GetVariableValue(self, variable, namelist, wt_type=None):
        """
        Get the value of a variable from a given namelist. 
        
        Notes:
        - NMR variables belong to a namelist with name None.
        - wt namelists need an additional identifier for the 'type' attribute.
        - Unset variables with no default will return as None.
        - String values require double literals (e.g. 'str' --> "'str'").
        """
        for nl in self.namelists.GetAllMatches(namelist):
            if nl.has_key(variable):
                # Most namelists are uniquely defined by their name.
                if wt_type is None:
                    return nl[variable]
                # There may be multiple wt namelists with different type's
                elif nl['type'] == wt_type:
                    return nl[variable]
        # If there's no match, return the default
        # NB: wt namelists have no default
        if wt_type is None:
            return self.namelists.GetFirstMatch(namelist).defaults[variable]
        else:
            return None

    def SetVariableValue(self, variable, value, namelist, wt_type=None):
        """
        Set the value of a variable from a given namelist. Return 1 if
        successful, 0 otherwise.
        
        Notes:
        - NMR variables belong to a namelist with name None.
        - wt namelists need an additional identifier for the 'type' attribute.
        - String values require double literals (e.g. 'str' --> "'str'").
        """
        for nl in self.namelists.GetAllMatches(namelist):
            if wt_type is None or nl['type'] == wt_type:
                nl[variable] = value
                return 1
        # If the namelist of this variable does not exist, it must be made.
        self.namelists.append(AmberNamelist(namelist,{variable:value}))
        if wt_type is not None: 
            self.namelists[-1]['type'] = wt_type
            return 1
        return 0

    def SetSinglePoint(self):
        """Change settings to run a single point energy calculation.
        """
        # Use debug mode for more verbose output
        # self.SetVariableValue('do_debugf',1,'debug')
        # self.SetVariableValue('dumpfrc',1,'debug')
        # Use dynamics, not minimization (can report thermodynamic state info).
        self.SetVariableValue('imin',0,'cntrl')
        self.SetVariableValue('nstlim',0,'cntrl')
        # Just print an energy report for step 0, no coordinate output.
        self.SetVariableValue('irest',0,'cntrl')
        self.SetVariableValue('ntpr',1,'cntrl')
        self.SetVariableValue('ntwx',0,'cntrl')

    def SetRestart(self, isRestart=True):
        """Change whether or not this mdin file is a restart.
        """
        if isRestart:
            self.SetVariableValue('irest',1,'cntrl')
            self.SetVariableValue('ntx',5,'cntrl')
        else:
            self.SetVariableValue('irest',0,'cntrl')

    def AddRestraints(self, disang, dumpave='fort.35', print_step=0):
        """Add NMR restraints to this mdin file.
        """
        self.SetVariableValue('nmropt',1,'cntrl')
        self.SetVariableValue('DISANG',disang,None)
        self.SetVariableValue('DUMPAVE',dumpave,None)
        self.SetVariableValue('type',"'DUMPFREQ'",'wt',"'type'")
        self.SetVariableValue('istep1',print_step,'wt',"'DUMPFREQ'")
        self.SetVariableValue('type',"'END'",'wt',"'type'")

    def WriteAmberMdinFile(self, outfile, max_chars_per_line=40):
        """
        Write a new mdin file with the current namelist information. For 
        clarity, any values that are set to the default will be omitted.
        """
        import sys
        if hasattr(outfile,'write'):
            pass
        elif isinstance(outfile,str):
            outfile = open(outfile,'w')
        else:
            raise TypeError("'outfile' must be either a string or file object.")

        outfile.write(self.title)
        for name in ['debug','cntrl','ewald','qmmm','pb']:
            try:
                # GetFirstMatch raises an exception if no matches occur
                # TODO: fix this?
                nl = self.namelists.GetFirstMatch(name)
                outfile.write(nl.SprintFortran(max_chars_per_line))
            except:
                pass
        end = ''
        for nl in self.namelists.GetAllMatches('wt'):
            if nl['type'] not in ["'end'","'END'","'End'"]:
                outfile.write(nl.SprintFortran(max_chars_per_line))
            else:
                end = nl.SprintFortran(max_chars_per_line)
        outfile.write(end)
        # TODO: same as above?
        try:
            nl = self.namelists.GetFirstMatch(None)
            for key in nl.keys(): outfile.write(' %s=%s\n'%(key,nl[key]))
        except:
            pass

        if outfile is not sys.__stdout__: outfile.close()
                  
class AmberNamelist(Namelist):
    """
    An AmberNamelist is simply a Fortran namelist that also has defaults 
    specific to the MD engine (sander by default). Rather than always deal with
    very large namelists, only the variable keys that differ from the defaults 
    can also be returned.
    """
    def __init__(self,name,*args,**kwargs):
        engine = 'sander'
        if kwargs.has_key('engine'): engine = kwargs.pop('engine')
        super(AmberNamelist,self).__init__(name,*args,**kwargs)
        self.SetDefaults(engine)

    def SprintFortran(self, max_chars_per_line=40):
        """Return a formatted string of the namelist readable in Fortran.
        """
        sprint = ' &%s\n'%self.name.strip()
        line = '  '
        for key in self.GetNonDefaults():
            toAdd = '%s=%s, '%(key,self[key])
            # add to the line until it exceeds the max character count
            if len(line + toAdd) <= max_chars_per_line:
                line += toAdd
            else:
                sprint += line + '\n'
                line = '  ' + toAdd
        # print remainder, hack to to strip off extra space and comma
        sprint += line[:-2] + '\n /\n'
        return sprint
    
    def GetNonDefaults(self):
        """Return keys only for variables that differ from the default value.
        """
        non_defaults = []
        for key,value in self.iteritems():
            if self.defaults.has_key(key):
                if value != self.defaults[key]:
                    non_defaults.append(key)
            else:
                non_defaults.append(key)
        return non_defaults

    def SetDefaults(self, engine):
        """Set the default variables for a particular MD engine (e.g. sander).
        """
        sander_defaults = {
            'debug' : {
                'do_debugf' : 0, 'dumpfrc' : 0, 
                'zerochg' : 0, 'zerovdw' : 0, 'zerodip' : 0,
                'do_dir' : 1, 'do_rec' : 1, 'do_adj' : 1, 'do_self' : 1, 
                'do_bond' : 1, 'do_cbond' : 1, 'do_angle' : 1, 'do_ephi' : 1,
                'do_xconst' : 1, 'do_cap' : 1
                },
            'cntrl' : {
                'irest': 0, 'ibelly' : 0, 'ntx' : 1, 'ntxo' : 1,
                'ntcx' : 0, 'ig' : 71277, 'tempi' : 0.0, 'ntb' : 1, 'ntt' : 0,
                'nchain' : 1, 'temp0' : 300.0, 'tautp' : 1.0, 'ntp' : 0,
                'pres0' : 1.0, 'comp' : 44.6, 'taup' : 1.0, 'nscm' : 1000,
                'nstlim' : 1, 'dt' : 0.001, 'ntc' : 1, 'ntcc' : 0, 'nconp' : 0,
                'tol' : 0.00001, 'ntf' : 1, 'nsnb' : 25, 'cut' : 8.0, 
                'dielc' : 0, 'ntpr' : 50, 'ntwx' : 0, 'ntwv' : 0, 'ntwe' : 0, 
                'ntave' : 0, 'ntpp' : 0, 'ioutfm' : 0, 'ntr' : 0, 'nrc' : 0, 
                'ntrx' : 1, 'taur' : 0, 'nmropt' : 0, 'ivcap' : 0, 
                'cutcap' : 0.0, 'xcap' : 0.0, 'ycap' : 0.0, 'zcap' : 0.0, 
                'fcap' : 1.5, 'xlorth' : -1.0, 'ylorth' : -1.0, 'zlorth' : -1.0,
                'xorth' : 47114711.0, 'yorth' : 47114711.0, 
                'zorth' : 47114711.0, 'forth' : 1.5, 'imin' : 0, 
                'drms' : 1.0e-4, 'dele' : 0, 'dx0' : 0.01, 'pencut' : 0.1, 
                'ipnlty' : 1, 'iscale' : 0, 'scalm' : 100.0, 'noeskp' : 1, 
                'maxcyc' : 1, 'ncyc' : 10, 'ntmin' : 1, 'vlimit' : 20.0, 
                'mxsub' : 1, 'ipol' : 0, 'jfastw' : 0, 'watnam' : '    ', 
                'owtnm' : '    ', 'hwtnm1' : '    ', 'hwtnm2' : '    ', 
                'iesp' : 0, 'skmin' : 50, 'skmax' : 100, 'vv' : 0, 'vfac' : 0, 
                'tmode' : 1, 'ips' : 0, 'mipsx' : -1, 'mipsy' : -1, 
                'mipsz' : -1, 'mipso' : 4, 'gridips' : 2, 'raips' : -1.0, 
                'dvbips' : 1.0e-8, 'isgld' : 0, 'isgsta' : 1, 'isgend' : 0, 
                'tsgavg' : 0.2, 'sgft' : 0.0, 'tempsg' : 1.0, 'jar' : 0, 
                'iamoeba' : 0, 'numexchg' : 0, 'repcrd' : 1, 'numwatkeep' : -1,
                'hybridgb' : 0, 'ntwprt' : 0, 'tausw' : 0.1, 'ntwr' : 500, 
                'iyammp' : 0, 'imcdo' : -1, 'igb' : 0, 'alpb' : 0, 
                'arad' : 15.0, 'rgbmax' : 25.0, 'saltcon' : 0.0, 
                'offset' : -999999.0, 'gbsa' : 0, 'vrand' : 1000, 
                'surften' : 0.005, 'iwrap' : 0, 'nrespa' : 1, 'nrespai' : 1, 
                'gamma_ln' : 0.0, 'extdiel' : 78.5, 'intdiel' : 1.0, 
                'cut_inner' : 8.0, 'icfe' : 0, 'clambda' : 0.0, 'klambda' : 1, 
                'rbornstat' : 0, 'lastrst' : 1, 'lastist' : 1, 'itgtmd' : 0, 
                'tgtrmsd' : 0, 'tgtmdfrc' : 0, 'tgtfitmask' : '', 
                'tgtrmsmask' : '', 'idecomp' : 0, 'temp0les' : -1.0,
                'restraintmask' : '', 'restraint_wt' : 0, 'bellymask' : '',
                'noshakemask' : '', 'crgmask' : '', 'iwrap_mask' : '',
                'mmtsb_switch' : 0, 'mmtsb_iterations' : 100, 'rdt' : 0.0,
                'icnstph' : 0, 'solvph' : 7.0, 'ntcnstph' : 10, 'ifqnt' : 0,
                'ievb' : 0, 'ipimd' : 0, 'itimass' : 0, 'ineb' : 0,
                'profile_mpi' : 0, 'ilscivr' : 0, 'icorf_lsc' : 0, 'ipb' : 0,
                'inp' : 2, 'gbneckscale' : -999999.0, 'gbalphah' : 0.788440,
                'gbbetah' : 0.798699, 'gbgammah' : 0.437334,
                'gbalphac' : 0.733756, 'gbbetac' : 0.506378,
                'gbgammac' : 0.205844, 'gbalphan' : 0.503364,
                'gbbetan' : 0.316828, 'gbgamman' : 0.192915,
                'gbalphaos' : 0.867814, 'gbbetaos' : 0.876635,
                'gbgammaos' : 0.387882, 'gbalphap' : 1.0, 'gbbetap' : 0.8,
                'gbgammap' : 4.851, 'sh' : 1.425952, 'sc' : 1.058554,
                'sn' : 0.733599, 'so' : 1.061039, 'ss' : -0.703469, 'sp' : 0.5,
                't' : 0.0, 'ntn' : 0 , 'scalpha' : 0.5, 'scbeta' : 12.0,
                'ifsc' : 0, 'scmask' : '', 'logdvdl' : 0, 'dvdl_norest' : 0,
                'dynlmb' : 0, 'ifmbar' : 0, 'bar_intervall' : 100,
                'bar_l_min' : 0.1, 'bar_l_max' : 0.9, 'bar_l_incr' : 0.1,
                'idssp' : 0, 'irism' : 0, 'restart_cmd' : '.false.',
                'eq_cmd' : '.false.', 'adiab_param' : 1.0, 'dec_verbose' : 3 
                },
            'ewald' : {
                'dsum_tol' : 1.0e-5,'ew_coeff' : 0.0,'skinnb' : 0.0,
                'diptol' : 1.0e-4,'dipmass' : 0.33,'diptau' : 11.0,'nfft1' : 0,
                'nfft2' : 0,'nfft3' : 0,'order' : 4,'opt_infl' : 1,
                'verbose' : 0,'nbflag' : 0,'nbtell' : 0,'netfrc' : 12344321,
                'ew_type' : 0,'vdwmeth' : 1,'eedmeth' : 0,'ee_type' : 1,
                'eedtbdns' : 5000.0,'rsum_tol' : 5.0e-5,'maxexp' : 0.0,
                'mlimit' : "0,0,0",'use_pme' : 0,'maxiter' : 20,'indmeth' : 3,
                'irstdip' : 0,'nquench' : 0,'frameon' : 1,'chngmask' : 1,
                'scaldip' : 1,'gridpointers' : 1,'column_fft' : 1 
                },
            'qmmm' : {
                'qmcut' : -1,'iqmatoms' : '','qmmask' : '','qmgb' : 2,
                'qm_theory' : '', 'qmcharge' : 0,'qmqmdx' : 1,
                'verbosity' : 0,'tight_p_conv' : 0,'scfconv' : 1.0e-8,
                'errconv' : 1e-1,'ndiis_matrices' :6,'ndiis_attempts' : 0,
                'printcharges' : 0, 'printdipole' : 0,'peptide_corr' : 0,
                'itrmax' : 1000,'qmshake' : 1,'qmqm_erep_incore' : 1,
                'qmmmrij_incore' : 1,'lnk_dis' : 1.09,'lnk_atomic_no' : 1,
                'lnk_method' : 1,'spin' : 1, 'pseudo_diag' : 1,
                'pseudo_diag_criteria' : 0.05,'qm_ewald' : 1,'qm_pme' : 1,
                'kmaxqx' : 8,'kmaxqy' : 8,'kmaxqz' : 8,'ksqmaxq' : 100,
                'writepdb' : 0,'qmmm_int' : 1,'adjust_q' : 2,
                'diag_routine' : 1, 'density_predict' : 0,
                'fock_predict' : 0, 'fockp_d1' : 2.4,'fockp_d2' : -1.2,
                'fockp_d3' : -0.8,'fockp_d4' : 0.6,'idc' : 0,'divpb' : 0,
                'dftb_maxiter' : 70,'dftb_disper' : 0,
                'dftb_3rd_order' : 'NONE','dftb_chg' : 0,
                'dftb_telec' : 0,'dftb_telec_step' : 0,
                'qmmm_omp_max_threads' : 1, 'chg_lambda' : 1,
                'nearest_qm_solvent' : 0,'nearest_qm_solvent_fq' : 1,
                'nearest_qm_solvent_resname' : 'WAT ',
                'qmmm_switch' : 0, 'r_switch_hi' : -1, 'r_switch_lo' : -1
                },
            'pb' : {
                'epsin' : 1.0, 'epsout' : 80.0, 'smoothopt' : 1,
                'istrng' : 0.0, 'pbtemp' : 300.0, 'radiopt' : -1,
                'dprob' : 1.4, 'iprob' : 2.0, 'npbopt' : 0, 'solvopt' : 1,
                'accept' : 0.001, 'maxitn' : 100, 'fillratio' : 2.0,
                'space' : 0.5, 'nbuffer' : 0, 'nfocus' : 2, 'fscale' : 8,
                'npbgrid' : 1, 'arcres' : 0.25,'dbfopt' : 1,'bcopt' : 5,
                'scalec' : 0,'eneopt' : 2,'frcopt' : 0,'cutfd' : 5.0,
                'cutnb' : 0.0, 'nsnbr' : 1, 'nsnba' : 1,'phiout' : 0,
                'phiform' : 0, 'npbverb' : 0, 'npopt' : 2, 'decompopt' : 2,
                'use_rmin' : 1, 'sprob' : 0.557, 'vprob' : 1.3,
                'rhow_effect' : 1.129, 'use_sav' : 1,
                'cavity_surften' : 0.0378, 'cavity_offset' : -0.5692,
                'maxsph' : 400, 'maxarc' : 512, 'cutsa' : 9.0, 'ndofd' : 1,
                'ndosas' : 1, 'fmiccg' : 0.9, 'ivalence' : 1.0,
                'laccept' : 0.1, 'wsor' : 1.9, 'lwsor' : 1.95,
                'pbkappa' : 0, 'radinc' : 0.8, 'expthresh' : 0.2,
                'offx' : 0.0, 'offy' : 0.0, 'offz' : 0.0,
                'sepbuf' : 4.0, 'mpopt' : 0, 'lmax' : 80, 'maxarcdot' : 1500
                },
            'wt' : {
                'type' : None, 'istep1' : 0, 'istep2' : 0, 
                'value1' : None, 'value2' : None, 'iinc' : None, 'imult' : None,
                'bond' : None, 'angle' : None, 'torsion' : None, 
                'improp' : None, 'vdw' : None, 'hb' : None, 'elec' : None,
                'nb' : None, 'attract' : None, 'repulse' : None, 'rstar' : None
                },
            None : {'DISANG' : None, 'DUMPAVE' : 'fort.35', 'LISTIN' : None, 
                    'LISTOUT' : None}
            }
            
        pmemd_defaults = {
            'cntrl' : {
                'imin' : 0, 'nmropt' : 0, 'ntx' : 1, 'irest' : 0,
                'ntrx' : 1, 'ntxo' : 1, 'ntpr' : 50, 'ntave' : 0, 'ntwr' : 500,
                'iwrap' : 0, 'ntwx' : 0, 'ntwv' : 0, 'ntwe' : 0, 'ioutfm' : 0,
                'ntwprt' : 0, 'ntf' : 1, 'ntb' : 1, 'dielc' : 1.0, 'cut' : 0.,
                'nsnb' : 25, 'ipol' : 0, 'igb' : 0, 'intdiel' : 1.0,
                'extdiel' : 78.5, 'saltcon' : 0.0, 'rgbmax' : 25.0,
                'rbornstat' : 0, 'offset' : 0.09, 'gbsa' : 0, 'surften' : 0.005,
                'nrespai' : 1, 'cut_inner' : 8.0, 'ibelly' : 0, 'ntr' : 0,
                'bellymask' : 0, 'maxcyc' : 1, 'ncyc' : 10, 'ntmin' : 1,
                'dx0' : 0.01, 'drms' : 0.0001, 'nstlim' : 1, 'nscm' : 1000,
                't' : 0.0, 'dt' : 0.001, 'nrespa' : 1, 'ntt' : 0, 
                'temp0' : 300., 'tempi' : 0.0, 'ig' : 71277, 'tautp' : 1.0, 
                'gamma_ln' : 0.0, 'vrand' : 1000, 'vlimit' : 20.0, 'ntp' : 0, 
                'pres0' : 1.0, 'comp' : 44.6, 'taup' : 1.0, 'ntc' : 1, 
                'tol' : 0.00001, 'jfastw' : 0, 'hwtnm1' : '    ', 
                'hwtnm2' : '    ', 'owtnm' : '    ', 'watnam' : '    ', 
                'ivcap' : 0, 'fcap' : 1.5, 'pencut' : 0.1, 'ndfmin' : 0, 
                'jar' : 0, 'no_intermolecular_bonds' : 1, 
                'ene_avg_sampling' : -1, 'mdinfo_flush_interval' : 60, 
                'mdout_flush_interval' : 300, 'dbg_atom_redistribution' : 0, 
                'loadbal_verbose' : 0, 'es_cutoff' : 0.0, 'vdw_cutoff' : 0.0,
                'dtemp' : 0, 'dxm' : 0, 'heat' : 0, 'alpb' : 0, 'arad' : 15.0 
                },
            'ewald' : {
                'nfft1' : 0, 'nfft2' : 0, 'nfft3' : 0, 'order' : 4,
                'verbose' : 0, 'ew_type' : 0, 'dsum_tol' : 1.0e-5, 
                'rsum_tol' : 5.0e-5, 'mlimit' : "0,0,0", 'ew_coeff' : 0.0, 
                'nbflag' : 1, 'skinnb' : 2, 'nbtell' : 0, 'netfrc' : -1, 
                'use_pme' : 1, 'vdwmeth' : 1, 'eedmeth' : 1, 'eedtbdns' : 5000.,
                'ee_type' : 0, 'frameon' : 1, 'chngmask' : 1, 'alpha' : 0.0, 
                'beta' : 0.0, 'gamma' : 0.0, 'a' : 0.0, 'b' : 0.0, 'c' : 0.0, 
                'use_axis_opt' : -1, 'fft_grids_per_ang' : 1.0, 
                'block_fft' : -1, 'fft_blk_y_divisor' : -1, 'excl_recip' : -1, 
                'excl_master' : -1, 'atm_redist_freq' : -1 
                },
            'qmmm' : {},
            'pb' : {},
            'wt' : {
                'type' : None, 'istep1' : 0, 'istep2' : 0, 
                'value1' : None, 'value2' : None, 'iinc' : None, 'imult' : None,
                'bond' : None, 'angle' : None, 'torsion' : None, 
                'improp' : None, 'vdw' : None, 'hb' : None, 'elec' : None,
                'nb' : None, 'attract' : None, 'repulse' : None, 'rstar' : None
                },
            None : {'DISANG' : None, 'DUMPAVE' : 'fort.35', 'LISTIN' : None, 
                    'LISTOUT' : None}
            }

        valid_namelists = sander_defaults.keys()
        if self.name in valid_namelists:
            if engine in ['sander','sander.MPI']: 
                self.defaults = sander_defaults[self.name]
            elif engine in ['pmemd','pmemd.MPI']:  
                self.defaults = pmemd_defaults[self.name]
            else: 
                msg = 'Can\'t set defaults for unknown engine %s'%engine
                raise ValueError(msg)
        else:
            raise ValueError('Unrecognized AMBER mdin namelist: %s'%self.name)

if __name__ == '__main__':
    import sys,os
    print '=== AmberMdin Test Suite ==='
    clean = False
    if len(sys.argv) < 2:
        print 'No mdin test file specified, writing/reading a stock mdin file.'
        # Tried to fit as many "tricky" mdin settings as I could here
        mdin_test = (' --title1--\n'
                     ' &cntrl\n'
                     '  irest = 1, ntx = 5, ! restart flag\n'
                     '  ifqnt = 1, nmropt = 1\n'
                     ' /\n'
                     ' --title2--\n'
                     ' &ewald\n'
                     '  dsum_tol = 1.0e-6, order = 6, \n'
                     ' /\n'
                     ' &qmmm\n'
                     "  qmmask = ':1-2,@10-11', qmcharge = -1\n"
                     '  ! a comment on its own line!\n'
                     ' /\n'
                     " &wt type='DUMPFREQ', istep = 1 /\n"
                     " &wt type='END' /\n"
                     ' DISANG=RST\n'
                     )
        test_name = 'test.mdin'
        clean = True
        out = open(test_name,'w')
        out.write(mdin_test)
        out.close()
    else:
        test_name = sys.argv[1]
        print 'Reading mdin from %s'%test_name
        mdin_test = ''.join(open(test_name).readlines())
        
    print 'mdin test file:'
    print '======'
    print mdin_test
    print '======'
   
    try:
        print '>>> mdin_obj = ReadAmberMdinFile(%s)'%test_name
        mdin_obj = ReadAmberMdinFile(test_name)
        print '>>> mdin_obj.WriteAmberMdinFile(sys.stdout)'
        print '======'
        mdin_obj.WriteAmberMdinFile(sys.stdout)
        print '======'
    except KeyError:
        pass
    if clean:
        print 'Cleaning up %s'%test_name
        os.remove(test_name)
