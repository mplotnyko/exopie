import numpy as np
from scipy.interpolate import interpn
from exopie.property import exoplanet, load_Data
from exopie.tools import chemistry
import warnings

class rocky(exoplanet):
    def __init__(self, Mass=[1,0.001], Radius=[1,0.001],  N=50000, **kwargs):
        super().__init__(N, Mass, Radius,**kwargs)
        xSi = kwargs.get('xSi', [0,0.2])
        xFe = kwargs.get('xFe', [0,0.2])
        self.set_xSi(a=xSi[0], b=xSi[1])
        self.set_xFe(a=xFe[0], b=xFe[1])

    def run(self):
        get_R = lambda x: interpn(PointsRocky, Radius_DataRocky, x) # x=cmf,Mass,xSi,xFe 
        self._check_rocky(get_R)
        args = np.asarray([self.Radius,self.Mass,self.xSi,self.xFe]).T
        residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],*param[1:]]).T))**2/1e-4 
        self.CMF = self._run_MC(residual,args)
        self.FeMF,self.SiMF,self.MgMF = chemistry(self.CMF,xSi=self.xSi,xFe=self.xFe)
    
class ocean(exoplanet):
    def __init__(self, Mass=[1,0.001], Radius=[1,0.001], N=50000, **kwargs):
        super().__init__(N, Mass, Radius, **kwargs)
        CMF = kwargs.get('CMF', [0.325,0.325])
        self.set_CMF(a=CMF[0], b=CMF[1])

    def run(self):
        get_R = lambda x: interpn(PointsWater, Radius_DataWater, x) # x=wmf,Mass,cmf   
        self._check_ocean(get_R)
        args = np.asarray([self.Radius,self.Mass,self.CMF]).T
        residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],param[1],param[2]*(1-x[0])]).T))**2/1e-4 
        self.WMF = self._run_MC(residual,args)

class envelope(exoplanet):
    def __init__(self, Mass=[1,0.001], Radius=[1,0.001], xSi=[0,0.], xFe=[0.,0.], atm_height=[20,30], N=50000, **kwargs):
        super().__init__(N, Mass, Radius, xSi, xFe, atm_height, **kwargs)
        self.set_xSi(a=xSi[0], b=xSi[1])
        self.set_xFe(a=xFe[0], b=xFe[1])
        self.set_atm_height(a=atm_height[0], b=atm_height[1])

    def run(self):
        get_R = lambda x: interpn(PointsRocky, Radius_DataRocky, x[:4].T)+x[4]/6.371e3 # x=cmf,Mass,xSi,xFe,atm_h
        pos = (self.Mass>10**-0.5) & (self.Mass<10**1.3)
        for item in  ['Mass','Radius','xSi','xFe','atm_height']:
            setattr(self, item, getattr(self, item)[pos])
        args = np.asarray([self.Radius,self.Mass,self.xSi,self.xFe,self.atm_height])
        residual = lambda x,param: np.sum(param[0]-_get_R(np.asarray([x[0],*param[1:]])))**2/1e-4 
        self.CMF = self._run_MC(residual,args)
        self.FeMF,self.SiMF,self.MgMF = chemistry(self.CMF,xSi=self.xSi,xFe=self.xFe)
    
PointsRocky,Radius_DataRocky,PointsWater,Radius_DataWater = load_Data() # load interpolation fits

def get_radius(M,cmf=0.325,wmf=0,xSi=0,xFe=0.1):
    '''
    Find the Radius of a planet, given mass and interior parameters.
    Mass (M) has to be in Earth masses and a float.
    '''
    n,M = [len(M),np.asarray(M)] if isinstance(M, (list, np.ndarray)) else [1,np.full(1,M)]
    wmf = np.full(n,wmf)
    cmf = np.full(n,cmf)
    xSi = np.full(n,xSi)
    xFe = np.full(n,xFe)
    

    if (wmf>0).any():
        xi = np.asarray([wmf,M,cmf*(1-wmf)]).T
        return interpn(PointsWater, Radius_DataWater, xi)[:n]
    xi = np.asarray([cmf,M,xSi,xFe]).T
    return interpn(PointsRocky, Radius_DataRocky, xi)[:n]