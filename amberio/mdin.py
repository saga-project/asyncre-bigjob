"""
A module for I/O of AMBER mdin files

This module provides a class for AMBER mdin files. AmberMdin objects are
essentially collections of Fortran namelists that also have knowledge of the 
default values expected by AMBER (depending on the MD engine being used).

References: 
    AMBER 12 Manual: ambermd.org/doc12/Amber12.pdf
    MMPBSA.py source code: ambermd.org/AmberTools-get.html

Exported Classes:
    AmberNamelist A Namelist (see the namelist module) with default 
                  values dependent on the specified MD engine.
    AmberMdin     A list of AmberNamelists
"""
import sys

from namelist import Namelist, NamelistCollection
from amberio.ambertools import AmberMdout

__author__ = 'Brian K. Radak. (BKR) - <radakb@biomaps.rutgers.edu>'

__all__ = ['AmberMdin', 'AmberNamelist']


class AmberMdin(object):
    """
    List of AmberNamelist objects providing easy access to the data in an
    AMBER mdin file.
    """
    def __init__(self, title='', cntrl={}, ewald={}, qmmm={}, pb={}, wts=[], 
                 nmr_vars={}, engine='sander'):
        self.title = str(title)
        self.cntrl = AmberNamelist('cntrl',engine,**cntrl)
        self.ewald = AmberNamelist('ewald',engine,**ewald)
        self.qmmm = AmberNamelist('qmmm',engine,**qmmm)
        self.pb = AmberNamelist('pb',engine,**pb)
        self.wts = AmberWtList(*[AmberNamelist('wt',engine,**w) for w in wts])
        self.nmr_vars = AmberNamelist(None,engine,**nmr_vars)
        self.engine = str(engine)

    @classmethod
    def from_mdin(cls, mdin_name, engine='sander'):
        """Read an AMBER mdin file and return an AmberMdin object."""
        #     Separate the namelists as a NamelistCollection and convert to a 
        # list of AmberNamelists (with assigned defaults). Parse the remaining 
        # lines for title information and NMR variables (a NoneType 
        # AmberNamelist).
        nls,other_lines = NamelistCollection.separate_nls(mdin_name)
        if nls.first_match('cntrl') is not None:
            cntrl = nls.first_match('cntrl')
        else:
            cntrl = {}
        if nls.first_match('ewald') is not None:
            ewald = nls.first_match('ewald')
        else:
            ewald = {}
        if nls.first_match('qmmm') is not None:
            qmmm = nls.first_match('qmmm')
        else:
            qmmm = {}
        if nls.first_match('pb') is not None:
            pb = nls.first_match('pb')
        else:
            pb = {}
        wts = [m for m in nls.matches('wt')]

        title = ''
        nmr_vars = {}
        for line in other_lines:
            try:
                line = line[:line.index('!')]
            except ValueError:
                pass
            tokens = line.split('=')
            if len(tokens) > 1:
                nmr_vars[tokens[0].strip()] = tokens[1].strip()
            else:
                title += line.strip() + '\n'
        return cls(title,cntrl,ewald,qmmm,pb,wts,nmr_vars,engine)

    @classmethod
    def from_mdout(mdout_name, engine='sander'):
        """
        Return an AmberMdin object from specifications in an AMBER mdout file.
        
        WARNING! Not all specifications are fully/properly stored in mdout 
        files. Some information may be missing (e.g. weight changes).
        """
        mdout = AmberMdout(mdout_name)
        mdin = cls(engine=engine)
        for k,v in mdout.properties.iteritems():
            #    The AmberMdout 'properties' attribute is just a simple dict of
            # namelist key words and values; it does not distinguish which 
            # namelist the pairs come from. For flexibility, AmberMdin objects 
            # allow assignment of ANY namelist value and will raise a KeyError 
            # only if no default value exists. In order to determine which 
            # namelist a key word belongs to, check for a default value, if 
            # that does not raise a KeyError, then assign the value, otherwise 
            # move on to another namelist.
            #
            try:
                mdin.cntrl[k]
                mdin.cntrl[k] = v
            except KeyError:
                pass
            try:
                mdin.ewald[k]
                mdin.ewald[k] = v
            except KeyError:
                pass
            try:
                mdin.qmmm[k]
                mdin.qmmm[k] = v
            except KeyError:
                pass
            try:
                mdin.pb[k]
                mdin.pb[k] = v
            except KeyError:
                pass
        return mdin

    def __str__(self):
        # Format mdin file output as follows:
        # 1 - title lines
        # 2 - unique namelists (i.e. not &wt namelists)
        # 3 - &wt namelists, followed by &wt namelists with type='END'
        # 4 - variables that belong to no namelist (e.g. DISANG)
        #
        txt = '%s\n'%self.title.rstrip()
        txt += str(self.cntrl)
        txt += str(self.ewald)
        txt += str(self.qmmm)
        txt += str(self.pb)
        txt += str(self.wts)
        txt += str(self.nmr_vars)
        return txt

    def __setattr__(self, name, value):
        object.__setattr__(self,name,value)
        if name == 'engine':
            self.cntrl.engine = self.engine
            self.ewald.engine = self.engine
            self.qmmm.engine = self.engine
            self.pb.engine = self.engine
            for wt in self.wts:
                wt.engine = self.engine
            self.nmr_vars.engine = self.engine

    def modify_or_add_wt(self, type, position = 0, **kwargs):
        """
        If a wt with the specified type exists, modify it with the keyword
        values in kwargs. By default, the first match (position = 0) is 
        modified. If the specified wt does not exist, it will be made; this is
        the same as calling add_wt(-1,**kwargs).
        """
        wt_exists = False
        for i,match in enumerate(self.wts.type_matches(type)):
            if position == i:
                wt_exists = True
                for k,v in kwargs.iteritems():
                    match[k] = v
                break
        if not wt_exists:
            kwargs['type'] = type
            self.add_wt(**kwargs)

    def add_wt(self, type, position = -1, **kwargs):
        """
        Add a new wt namelist with keyword values in kwargs at the specified
        position (default is to append). 

        NB: 'position' ignores the existence of any 'END' sections and these
        will automatically be reset after the new wt has been added.
        """
        kwargs['type'] = type
        self.wts.insert(position,AmberNamelist('wt',self.engine,**kwargs))
        self.wts.set_end()

    def write_mdin(self, outfile, mode='w'):
        """
        Write a new mdin file with the current namelist information. For 
        clarity, any values that are set to the default will be omitted.
        """
        if hasattr(outfile,'write'):
            pass
        elif isinstance(outfile,str):
            outfile = open(outfile,mode)
        else:
            raise TypeError("'outfile' must be a string or file object.")
        outfile.write(str(self))
        if outfile != sys.stdout:
            outfile.close()

     
class AmberWtList(NamelistCollection):
    """
    A list of AMBER &wt namelists. Rather than distinguishing namelists by 
    their name (as in a NamelistCollection), the added methods use the 'type'
    attribute to search.
    """
    def append(self, item):
        if not isinstance(item,AmberNamelist) or item.name != 'wt':
            raise TypeError('AmberWtLists must contain &wt AmberNamelists!')
        NamelistCollection.append(self,item)

    def set_end(self):
        """Remove extra 'END' specifications and add a new one."""
        for i,nl in enumerate(self):
            if nl['type'] == "'END'":
                self.pop(i)
        self.append(AmberNamelist('wt','sander',**{'type': "'END'"}))

    def type_matches(self, *types):
        """
        Return a generator of all wt namelists whose type matches any of the
        arguments.
        """
        try:
            types[0]
        except IndexError:
            types = [types]
        for nl in self:
            for type in types:
                if nl['type'] == type:
                    yield nl

    def first_type_match(self, type):
        """
        Return the first result from type_matches(). If there are no matches,
        create a new namelist of that type and reset any 'END' specifications.
        """
        try:
            return next(self.type_matches(*type))
        except StopIteration:
            # This happens if type_matches() returns an empty generator.
            self.append(AmberNamelist('wt','sander',**{'type': "%s"%type}))
            self.set_end()
            # Not sure why calling self.first_type_match(type) here causes an
            # infinite recursion problem...
            #
            for nl in self:
                if nl['type'] == type:
                    return nl
           

class AmberNamelist(Namelist):
    """
    An AmberNamelist is simply a Fortran namelist that also has defaults
    specific to the MD engine.

    The current implementation only stores values that are explicitly set. If
    an unset value is requested, a default value is set and returned (this may 
    break if the default depends on another variable).

    In order to be extendible, there is no checking for compatibility with 
    different engines (e.g. you can set sander options for pmemd, even though
    the input will fail). However, engine dependent defaults are supported 
    (see above caveat).
    """
    def __init__(self, name, engine, *args, **kwargs):
        Namelist.__init__(self,name,line_prefix=' ',name_value_separator=' = ',
                          value_separator=', ',max_namevalues_per_line=72,
                          max_chars_per_line=72,*args,**kwargs)
        self.engine = engine
    
    def __str__(self):
        if len(self.keys()) > 0: # Don't print empty namelists.
            if self.name is not None:
                return Namelist.__str__(self)
            else:
                # Modify the Namelist print behavior for a NoneType namelist.
                # Namely, don't print the "&name" section and ending slash.
                return '\n'.join([' %s=%s'%(k,v) for k,v in self.iteritems()])
        else:
            return ''

    def __getitem__(self, key):
        try:
            return Namelist.__getitem__(self,key)
        except KeyError:
            self[key] = value = self.default(key)
            return value

    def default(self, key):
        """Set the default variables for the MD engine (e.g. sander)."""
        sander_defaults = {
            'cntrl': {
                'irest': 0, 'ibelly': 0, 'ntx': 1, 'ntxo': 1,
                'ntcx': 0, 'ig': 71277, 'tempi': 0.0, 'ntb': 1, 'ntt': 0,
                'nchain': 1, 'temp0': 300.0, 'tautp': 1.0, 'ntp': 0,
                'pres0': 1.0, 'comp': 44.6, 'taup': 1.0, 'nscm': 1000,
                'nstlim': 1, 'dt': 0.001, 'ntc': 1, 'ntcc': 0, 'nconp': 0,
                'tol': 0.00001, 'ntf': 1, 'nsnb': 25, 'cut': 8.0, 
                'dielc': 0, 'ntpr': 50, 'ntwx': 0, 'ntwv': 0, 'ntwe': 0, 
                'ntave': 0, 'ntpp': 0, 'ioutfm': 0, 'ntr': 0, 'nrc': 0, 
                'ntrx': 1, 'taur': 0, 'nmropt': 0, 'ivcap': 0, 
                'cutcap': 0.0, 'xcap': 0.0, 'ycap': 0.0, 'zcap': 0.0, 
                'fcap': 1.5, 'xlorth': -1.0, 'ylorth': -1.0, 'zlorth': -1.0,
                'xorth': 47114711.0, 'yorth': 47114711.0, 
                'zorth': 47114711.0, 'forth': 1.5, 'imin': 0, 
                'drms': 1.0e-4, 'dele': 0, 'dx0': 0.01, 'pencut': 0.1, 
                'ipnlty': 1, 'iscale': 0, 'scalm': 100.0, 'noeskp': 1, 
                'maxcyc': 1, 'ncyc': 10, 'ntmin': 1, 'vlimit': 20.0, 
                'mxsub': 1, 'ipol': 0, 'jfastw': 0, 'watnam': '    ', 
                'owtnm': '    ', 'hwtnm1': '    ', 'hwtnm2': '    ', 
                'iesp': 0, 'skmin': 50, 'skmax': 100, 'vv': 0, 'vfac': 0, 
                'tmode': 1, 'ips': 0, 'mipsx': -1, 'mipsy': -1, 
                'mipsz': -1, 'mipso': 4, 'gridips': 2, 'raips': -1.0, 
                'dvbips': 1.0e-8, 'isgld': 0, 'isgsta': 1, 'isgend': 0, 
                'tsgavg': 0.2, 'sgft': 0.0, 'tempsg': 1.0, 'jar': 0, 
                'iamoeba': 0, 'numexchg': 0, 'repcrd': 1, 'numwatkeep': -1,
                'hybridgb': 0, 'ntwprt': 0, 'tausw': 0.1, 'ntwr': 500, 
                'iyammp': 0, 'imcdo': -1, 'igb': 0, 'alpb': 0, 
                'arad': 15.0, 'rgbmax': 25.0, 'saltcon': 0.0, 
                'offset': -999999.0, 'gbsa': 0, 'vrand': 1000, 
                'surften': 0.005, 'iwrap': 0, 'nrespa': 1, 'nrespai': 1, 
                'gamma_ln': 0.0, 'extdiel': 78.5, 'intdiel': 1.0, 
                'cut_inner': 8.0, 'icfe': 0, 'clambda': 0.0, 'klambda': 1, 
                'rbornstat': 0, 'lastrst': 1, 'lastist': 1, 'itgtmd': 0, 
                'tgtrmsd': 0, 'tgtmdfrc': 0, 'tgtfitmask': '', 
                'tgtrmsmask': '', 'idecomp': 0, 'temp0les': -1.0,
                'restraintmask': '', 'restraint_wt': 0, 'bellymask': '',
                'noshakemask': '', 'crgmask': '', 'iwrap_mask': '',
                'mmtsb_switch': 0, 'mmtsb_iterations': 100, 'rdt': 0.0,
                'icnstph': 0, 'solvph': 7.0, 'ntcnstph': 10, 'ifqnt': 0,
                'ievb': 0, 'ipimd': 0, 'itimass': 0, 'ineb': 0,
                'profile_mpi': 0, 'ilscivr': 0, 'icorf_lsc': 0, 'ipb': 0,
                'inp': 2, 'gbneckscale': -999999.0, 'gbalphah': 0.788440,
                'gbbetah': 0.798699, 'gbgammah': 0.437334,
                'gbalphac': 0.733756, 'gbbetac': 0.506378,
                'gbgammac': 0.205844, 'gbalphan': 0.503364,
                'gbbetan': 0.316828, 'gbgamman': 0.192915,
                'gbalphaos': 0.867814, 'gbbetaos': 0.876635,
                'gbgammaos': 0.387882, 'gbalphap': 1.0, 'gbbetap': 0.8,
                'gbgammap': 4.851, 'sh': 1.425952, 'sc': 1.058554,
                'sn': 0.733599, 'so': 1.061039, 'ss': -0.703469, 'sp': 0.5,
                't': 0.0, 'ntn': 0 , 'scalpha': 0.5, 'scbeta': 12.0,
                'ifsc': 0, 'scmask': '', 'logdvdl': 0, 'dvdl_norest': 0,
                'dynlmb': 0, 'ifmbar': 0, 'bar_intervall': 100,
                'bar_l_min': 0.1, 'bar_l_max': 0.9, 'bar_l_incr': 0.1,
                'idssp': 0, 'irism': 0, 'restart_cmd': '.false.',
                'eq_cmd': '.false.', 'adiab_param': 1.0, 'dec_verbose': 3 
                },
            'ewald': {
                'dsum_tol': 1.0e-5,'ew_coeff': 0.0,'skinnb': 0.0,
                'diptol': 1.0e-4,'dipmass': 0.33,'diptau': 11.0,'nfft1': 0,
                'nfft2': 0,'nfft3': 0,'order': 4,'opt_infl': 1,
                'verbose': 0,'nbflag': 0,'nbtell': 0,'netfrc': 12344321,
                'ew_type': 0,'vdwmeth': 1,'eedmeth': 0,'ee_type': 1,
                'eedtbdns': 5000.0,'rsum_tol': 5.0e-5,'maxexp': 0.0,
                'mlimit': "0,0,0",'use_pme': 0,'maxiter': 20,'indmeth': 3,
                'irstdip': 0,'nquench': 0,'frameon': 1,'chngmask': 1,
                'scaldip': 1,'gridpointers': 1,'column_fft': 1 
                },
            'qmmm': {
                'qmcut': -1,'iqmatoms': '','qmmask': '','qmgb': 2,
                'qm_theory': '', 'qmcharge': 0,'qmqmdx': 1,
                'verbosity': 0,'tight_p_conv': 0,'scfconv': 1.0e-8,
                'errconv': 1e-1,'ndiis_matrices':6,'ndiis_attempts': 0,
                'printcharges': 0, 'printdipole': 0,'peptide_corr': 0,
                'itrmax': 1000,'qmshake': 1,'qmqm_erep_incore': 1,
                'qmmmrij_incore': 1,'lnk_dis': 1.09,'lnk_atomic_no': 1,
                'lnk_method': 1,'spin': 1, 'pseudo_diag': 1,
                'pseudo_diag_criteria': 0.05,'qm_ewald': 1,'qm_pme': 1,
                'kmaxqx': 8,'kmaxqy': 8,'kmaxqz': 8,'ksqmaxq': 100,
                'writepdb': 0,'qmmm_int': 1,'adjust_q': 2,
                'diag_routine': 1, 'density_predict': 0,
                'fock_predict': 0, 'fockp_d1': 2.4,'fockp_d2': -1.2,
                'fockp_d3': -0.8,'fockp_d4': 0.6,'idc': 0,'divpb': 0,
                'dftb_maxiter': 70,'dftb_disper': 0,
                'dftb_3rd_order': 'NONE','dftb_chg': 0,
                'dftb_telec': 0,'dftb_telec_step': 0,
                'qmmm_omp_max_threads': 1, 'chg_lambda': 1,
                'nearest_qm_solvent': 0,'nearest_qm_solvent_fq': 1,
                'nearest_qm_solvent_resname': 'WAT ',
                'qmmm_switch': 0, 'r_switch_hi': -1, 'r_switch_lo': -1
                },
            'pb': {
                'epsin': 1.0, 'epsout': 80.0, 'smoothopt': 1,
                'istrng': 0.0, 'pbtemp': 300.0, 'radiopt': -1,
                'dprob': 1.4, 'iprob': 2.0, 'npbopt': 0, 'solvopt': 1,
                'accept': 0.001, 'maxitn': 100, 'fillratio': 2.0,
                'space': 0.5, 'nbuffer': 0, 'nfocus': 2, 'fscale': 8,
                'npbgrid': 1, 'arcres': 0.25,'dbfopt': 1,'bcopt': 5,
                'scalec': 0,'eneopt': 2,'frcopt': 0,'cutfd': 5.0,
                'cutnb': 0.0, 'nsnbr': 1, 'nsnba': 1,'phiout': 0,
                'phiform': 0, 'npbverb': 0, 'npopt': 2, 'decompopt': 2,
                'use_rmin': 1, 'sprob': 0.557, 'vprob': 1.3,
                'rhow_effect': 1.129, 'use_sav': 1,
                'cavity_surften': 0.0378, 'cavity_offset': -0.5692,
                'maxsph': 400, 'maxarc': 512, 'cutsa': 9.0, 'ndofd': 1,
                'ndosas': 1, 'fmiccg': 0.9, 'ivalence': 1.0,
                'laccept': 0.1, 'wsor': 1.9, 'lwsor': 1.95,
                'pbkappa': 0, 'radinc': 0.8, 'expthresh': 0.2,
                'offx': 0.0, 'offy': 0.0, 'offz': 0.0,
                'sepbuf': 4.0, 'mpopt': 0, 'lmax': 80, 'maxarcdot': 1500
                },
            'wt': {'istep1': 0, 'istep2': 0},
            None: {}
            }
            
        pmemd_defaults = {
            'cntrl': {
                'imin': 0, 'nmropt': 0, 'ntx': 1, 'irest': 0,
                'ntrx': 1, 'ntxo': 1, 'ntpr': 50, 'ntave': 0, 'ntwr': 500,
                'iwrap': 0, 'ntwx': 0, 'ntwv': 0, 'ntwe': 0, 'ioutfm': 0,
                'ntwprt': 0, 'ntf': 1, 'ntb': 1, 'dielc': 1.0, 'cut': 0.,
                'nsnb': 25, 'ipol': 0, 'igb': 0, 'intdiel': 1.0,
                'extdiel': 78.5, 'saltcon': 0.0, 'rgbmax': 25.0,
                'rbornstat': 0, 'offset': 0.09, 'gbsa': 0, 'surften': 0.005,
                'nrespai': 1, 'cut_inner': 8.0, 'ibelly': 0, 'ntr': 0,
                'bellymask': 0, 'maxcyc': 1, 'ncyc': 10, 'ntmin': 1,
                'dx0': 0.01, 'drms': 0.0001, 'nstlim': 1, 'nscm': 1000,
                't': 0.0, 'dt': 0.001, 'nrespa': 1, 'ntt': 0, 
                'temp0': 300., 'tempi': 0.0, 'ig': 71277, 'tautp': 1.0, 
                'gamma_ln': 0.0, 'vrand': 1000, 'vlimit': 20.0, 'ntp': 0, 
                'pres0': 1.0, 'comp': 44.6, 'taup': 1.0, 'ntc': 1, 
                'tol': 0.00001, 'jfastw': 0, 'hwtnm1': '    ', 
                'hwtnm2': '    ', 'owtnm': '    ', 'watnam': '    ', 
                'ivcap': 0, 'fcap': 1.5, 'pencut': 0.1, 'ndfmin': 0, 
                'jar': 0, 'no_intermolecular_bonds': 1, 
                'ene_avg_sampling': -1, 'mdinfo_flush_interval': 60, 
                'mdout_flush_interval': 300, 'dbg_atom_redistribution': 0, 
                'loadbal_verbose': 0, 'es_cutoff': 0.0, 'vdw_cutoff': 0.0,
                'dtemp': 0, 'dxm': 0, 'heat': 0, 'alpb': 0, 'arad': 15.0 
                },
            'ewald': {
                'nfft1': 0, 'nfft2': 0, 'nfft3': 0, 'order': 4,
                'verbose': 0, 'ew_type': 0, 'dsum_tol': 1.0e-5, 
                'rsum_tol': 5.0e-5, 'mlimit': "0,0,0", 'ew_coeff': 0.0, 
                'nbflag': 1, 'skinnb': 2, 'nbtell': 0, 'netfrc': -1, 
                'use_pme': 1, 'vdwmeth': 1, 'eedmeth': 1, 'eedtbdns': 5000.,
                'ee_type': 0, 'frameon': 1, 'chngmask': 1, 'alpha': 0.0, 
                'beta': 0.0, 'gamma': 0.0, 'a': 0.0, 'b': 0.0, 'c': 0.0, 
                'use_axis_opt': -1, 'fft_grids_per_ang': 1.0, 
                'block_fft': -1, 'fft_blk_y_divisor': -1, 'excl_recip': -1, 
                'excl_master': -1, 'atm_redist_freq': -1 
                },
            'qmmm': {},
            'pb': {},
            'wt': {'istep1': 0, 'istep2': 0},
            None: {}
            }
        if self.engine in ['sander','sander.MPI']:
            defaults = sander_defaults
        else:
            defaults = pmemd_defaults
        try:
            defaults[self.name]
            try:
                return defaults[self.name][key]
            except KeyError:
                raise KeyError('No default for variable %s in namelist %s'
                               %(key,self.name))
        except KeyError:
            raise KeyError('No defaults for unknown namelist %s'%self.name)

if __name__ == '__main__':
    import sys
    import os

    print '=== AmberMdin Test Suite ==='
    clean = False
    if len(sys.argv) < 2:
        print 'No mdin test file specified, writing/reading a stock mdin file.'
        # Tried to fit as many "tricky" mdin settings as I could here
        mdin_test = (' --title1--\n'
                     ' &cntrl\n'
                     '  irest = 1, ntx = 5, ! restart flag\n'
                     '  ifqnt = 1, nmropt = 1,\n'
                     '  ntt = 3, temp0 = 300.0, gamma_ln = 1.0,\n'
                     ' /\n'
                     ' --title2--\n'
                     ' &ewald\n'
                     '  dsum_tol = 1e-6\n'
                     ' /\n'
                     ' &qmmm\n'
                     "  qmmask = ':1-2,@10-11', qmcharge = -1\n"
                     '  ! a comment on its own line!\n'
                     ' /\n'
                     " &wt type='DUMPFREQ', istep1 = 1 /\n"
                     " &wt type='END' /\n"
                     ' DISANG=foo.RST\n'
                     ' LISTIN=POUT')
        test_name = 'test.mdin'
        clean = True
        out = open(test_name,'w')
        out.write(mdin_test)
        out.close()
    else:
        test_name = sys.argv[1]
        print 'Reading mdin from %s'%test_name
        mdin_test = ''.join(open(test_name).readlines())
        
    print '=== mdin test file:'
    print mdin_test
    print '==================='

    try:
        print '>>> mdin_obj = AmberMdin.from_mdin(%s)'%test_name
        mdin_obj = AmberMdin.from_mdin(test_name)
        mdin_obj.add_wt("'TEMP0'",0,**{'istep1': 0})
        mdin_obj.modify_or_add_wt("'DUMPFREQ'",0,**{'istep1': 10})
        mdin_obj.modify_or_add_wt("'TEMP0'",1,**{'istep1': 30000})

        mdin_obj.nmr_vars['DISANG'] = 'bar.RST'
        print '>>> print mdin_obj'
        print mdin_obj
    finally:
        if clean:
            os.remove(test_name)
