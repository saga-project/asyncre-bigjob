from amberrun import AmberRun
class AmberHamiltonian(AmberRun):
    def __init__(self,mode=None,basename=None,**filenames):
        AmberRun.__init__(self,mode=None,basename=None,**filenames)
        
