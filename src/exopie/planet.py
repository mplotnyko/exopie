import numpy as np
from scipy.interpolate import interpn
from exopie.property import exoplanet, load_Data
from exopie.tools import chemistry
import warnings

PointsRocky,Radius_DataRocky,PointsWater,Radius_DataWater = load_Data() # load interpolation fits

class rocky(exoplanet):
    def __init__(self, Mass=[1,0.001], Radius=[1,0.001],  N=50000, **kwargs):
        '''
        Rocky planet method. CMF, WMF and atm_height will be ignored.
        '''
        super().__init__(N, Mass, Radius,**kwargs)
        xSi = kwargs.get('xSi', [0,0.2])
        xFe = kwargs.get('xFe', [0,0.2])
        self.set_xSi(a=xSi[0], b=xSi[1])
        self.set_xFe(a=xFe[0], b=xFe[1])
        self._save_parameters = ['Mass','Radius','CMF','xSi','xFe','FeMF','SiMF','MgMF']

    def run(self):
        '''
        Run the rocky planet model.

        Attributes:
        --------
        self.CMF: array
            Core mass fraction 
        self.FeMF: array
            Iron mass fraction
        self.SiMF: array
            Silicon mass fraction
        self.MgMF: array
            Magnesium mass fraction
        '''
        get_R = lambda x: interpn(PointsRocky, Radius_DataRocky, x) # x=cmf,Mass,xSi,xFe 
        self._check_rocky(get_R)
        args = np.asarray([self.Radius,self.Mass,self.xSi,self.xFe]).T
        residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],*param[1:]]).T))**2/1e-4 
        self.CMF = self._run_MC(residual,args)
        self.FeMF,self.SiMF,self.MgMF = chemistry(self.CMF,xSi=self.xSi,xFe=self.xFe)
    
class water(exoplanet):
    def __init__(self, Mass=[1,0.001], Radius=[1,0.001], N=50000, **kwargs):
        '''
        Water planet method. xSi, xFe and atm_height will be ignored.
        '''
        super().__init__(N, Mass, Radius, **kwargs)
        CMF = kwargs.get('CMF', [0.325,0.325])
        self.set_CMF(a=CMF[0], b=CMF[1])
        self._save_parameters = ['Mass','Radius','WMF','CMF']

    def run(self):
        '''
        Run the water planet model.

        Attributes:
        --------
        self.WMF: array
            Water mass fraction
        self.CMF: array
            Rocky core mass fraction (cmf = rcmf*(1-wmf))
        '''
        get_R = lambda x: interpn(PointsWater, Radius_DataWater, x) # x=wmf,Mass,cmf   
        self._check_water(get_R)
        args = np.asarray([self.Radius,self.Mass,self.CMF]).T
        residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],param[1],param[2]*(1-x[0])]).T))**2/1e-4 
        self.WMF = self._run_MC(residual,args)

class envelope(exoplanet):
    def __init__(self, Mass=[1,0.001], Radius=[1,0.001], atm_height=[20,30], N=50000, **kwargs):
        '''
        Envelope planet method (beta). CMF and WMF will be ignored.
        '''
        super().__init__(N, Mass, Radius, xSi, xFe, atm_height, **kwargs)
        xSi = kwargs.get('xSi', [0,0.])
        xFe = kwargs.get('xFe', [0,0.])
        self.set_xSi(a=xSi[0], b=xSi[1])
        self.set_xFe(a=xFe[0], b=xFe[1])
        self.set_atm_height(a=atm_height[0], b=atm_height[1])
        self._save_parameters = ['Mass','Radius','CMF','atm_height']

    def run(self):
        '''
        Run the envelope planet model.

        Attributes:
        --------
        self.CMF: array
            Core mass fraction
        self.atm_height: array
            Height of the atmosphere (km)
        '''
        get_R = lambda x: interpn(PointsRocky, Radius_DataRocky, x[:4].T)+x[4]/6.371e3 # x=cmf,Mass,xSi,xFe,atm_h
        pos = (self.Mass>10**-0.5) & (self.Mass<10**1.3)
        for item in  ['Mass','Radius','xSi','xFe','atm_height']:
            setattr(self, item, getattr(self, item)[pos])
        args = np.asarray([self.Radius,self.Mass,self.xSi,self.xFe,self.atm_height])
        residual = lambda x,param: np.sum(param[0]-_get_R(np.asarray([x[0],*param[1:]])))**2/1e-4 
        self.CMF = self._run_MC(residual,args)
    
def get_radius(M,cmf=0.325,wmf=None,xSi=0,xFe=0.1):
    '''
    Find the Radius of a planet, given mass and interior parameters.
    
    Parameters:
    -----------
    M: float or array
        Mass of the planet in Earth masses, 
        if array the same interior parameters will be used for all masses.
    cmf: float
        Core mass fraction. 
    wmf: float
        Water mass fraction.
        xSi and xFe will be ignored and cmf corresponds to rocky portion only (rcmf).
        Thus rcmf is will keep the mantle to core fraction constant, rather than the total core mass.
    xSi: float
        Molar fraction of silicon in the core (between 0-0.2).
    xFe: float
        Molar fraction of iron in the mantle (between 0-0.2).
    
    Returns:
    --------
    Radius: float or array
        Radius of the planet in Earth radii.
    '''
    if isinstance(M, (list, np.ndarray)):
        n = len(M)
        wmf = np.full(n,wmf)
        cmf = np.full(n,cmf)
        xSi = np.full(n,xSi)
        xFe = np.full(n,xFe)
    if wmf is None:
        xi = np.asarray([cmf, M, xSi, xFe]).T
        result = interpn(PointsRocky, Radius_DataRocky, xi)
    else:
        xi = np.asarray([wmf, M, cmf * (1 - wmf)]).T
        result = interpn(PointsWater, Radius_DataWater, xi)
    return result if isinstance(M, (list, np.ndarray)) else result[0]

def get_cmf(M,R,xSi=0,xFe=0.1):
    '''
    Find the Core Mass Fraction of a planet, given mass and radius.
    
    Parameters:
    -----------
    M: float or array
        Mass of the planet in Earth masses, 
    R: float or array
        Radius of the planet in Earth radii.
    xSi: float
        Molar fraction of silicon in the core (between 0-0.2).
    xFe: float
        Molar fraction of iron in the mantle (between 0-0.2).
    
    Returns:
    --------
    cmf: float or array
        Core mass fraction of the planet.
    '''
    residual = lambda x,param: (param[0]-get_radius(param[1],cmf=x,xSi=param[2],xFe=param[3]))**2/1e-4
    if isinstance(M, (list, np.ndarray)):
        res = []
        for i in range(M):
            args = [R[i],M[i],xSi,xFe]
            res.append(minimize(residual,0.325,args=args,bounds=[[0,1]]).x[0])
    else:
        args = [R,M,xSi,xFe]
        res = minimize(residual,0.325,args=args,bounds=[[0,1]]).x[0]
    return res

def get_wmf(M,R,cmf=0.325):
    '''
    Find the Water Mass Fraction of a planet, given mass and radius.

    Parameters:
    -----------
    M: float or array
        Mass of the planet in Earth masses,
    R: float or array
        Radius of the planet in Earth radii.
    cmf: float
        Core mass fraction (rocky portion only).
    
    Returns:
    --------
    wmf: float or array
        Water mass fraction of the planet.
    '''
    residual = lambda x,param: (param[0]-get_radius(param[1],cmf=param[2],wmf=x[0]))**2/1e-4
    if isinstance(M, (list, np.ndarray)):
        res = []
        for i in range(M):
            args = [R[i],M[i],cmf]
            res.append(minimize(residual,0,args=args,bounds=[[0,1]]).x[0])
    else:
        args = [R,M,cmf]
        res = minimize(residual,0,args=args,bounds=[[0,1]]).x[0]
    return res

def get_mass(R,cmf=0.325,wmf=None,xSi=0,xFe=0.1):
    '''
    Find the Mass of a planet, given radius and interior parameters.

    Parameters:
    -----------
    R: float or array
        Radius of the planet in Earth radii.
        if array the same interior parameters will be used for all masses.
    cmf: float
        Core mass fraction. 
    wmf: float
        Water mass fraction.
        xSi and xFe will be ignored and cmf corresponds to rocky portion only (rcmf).
        Thus rcmf is will keep the mantle to core fraction constant, rather than the total core mass.
    xSi: float
        Molar fraction of silicon in the core (between 0-0.2).
    xFe: float
        Molar fraction of iron in the mantle (between 0-0.2).
    
    Returns:
    --------
    Mass: float or array
        Mass of the planet in Earth masses.
    '''
    residual = lambda x,param: (param[0]-get_radius(x[0],cmf=param[1],wmf=param[2],xSi=param[3],xFe=param[4]))**2/1e-4
    if isinstance(R, (list, np.ndarray)):
        res = []
        for i in range(R):
            args = [R[i],cmf,wmf,xSi,xFe]
            res.append(minimize(residual,1,args=args,bounds=[[10**-0.5,10**1.3]]).x[0])
    else:
        args = [R,cmf,wmf,xSi,xFe]
        print(args)
        res = minimize(residual,1,args=args,bounds=[[10**-0.5,10**1.3]]).x[0]
    return res
