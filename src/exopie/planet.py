import numpy as np
from scipy.interpolate import interpn
from scipy.optimize import minimize
from exopie.property import exoplanet, host_star, load_Data
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

    def run(self,star=None,ratio=None,star_norm=None):
        '''
        Run the rocky planet model.

        Parameters:
        -----------
        star: list [Fe/H, Mg/H, Si/H]
            Stellar abundances of Fe/H, Mg/H, Si/H.
        ratio: str
            ratio to constrain the planet chemistry to the star.
            e.g. ratio='Fe/Si,Fe/Mg,Mg/Si' (default=None)
        star_norm: list [Fe/H, Mg/H, Si/H]
            Normalization reference for the stellar abundances. 
        
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
        if star is None:
            residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],*param[1:]]).T))**2/1e-4 
            self.CMF = self._run_MC(residual,args)
            self.FeMF,self.SiMF,self.MgMF,_,_,_ = chemistry(self.CMF,xSi=self.xSi,xFe=self.xFe)
        elif ratio is None:
            warnings.warn('No target ratio provided. Running without stellar constraint.')
            residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],*param[1:]]).T))**2/1e-4
            self.CMF = self._run_MC(residual,args)
            self.FeMF,self.SiMF,self.MgMF,_,_,_ = chemistry(self.CMF,xSi=self.xSi,xFe=self.xFe)
        else:
            if star_norm is None:
                star_norm = [7.46,7.55,7.51] # Fe, Mg, Si Asplund 2021
            mu = [55.85e-3,28.09e-3,24.31e-3] # Fe, Mg, Si atomic masses
            star_w = [10**(star[i]+star_norm[i]-12)*mu[i] for i in range(3)]
            
            ratio_split = ratio.lower().split(',') # Fe/Si, Fe/Mg, Mg/Si ....
            dr_star = {'fe': star_w[0],'mg': star_w[1], 'si': star_w[2]}
            print('Using stellar constraints:',end=' ')
            [print(item+f' ({eval(item, dr_star.copy()):.2f})',end=', ') for item in ratio_split]
            print()
            
            def residual(x, param):
                radius_residual = np.sum(param[0] - get_R(np.asarray([x[0],param[1],x[1],x[2]]).T))**2 / 1e-4
                data = chemistry(x[0], xSi=x[1], xFe=x[2],xWu=x[3])
                dr_planet = {'fe': data[0], 'si': data[1], 'mg': data[2]}
                chem_residual = 0
                for item in ratio_split:
                    chem_residual += np.sum(eval(item, dr_star.copy())-eval(item, dr_planet.copy()))**2/1e-4
                return radius_residual + chem_residual
            
            args = np.asarray([self.Radius,self.Mass]).T
            self.CMF,self.xSi,self.xFe,self.xWu = self._run_MC(residual,args,
                                xi=[0.325,0.1,0.1,0.2],bounds=[[0,1],[0,0.2],[0,0.2],[0,0.5]]).T
            self.FeMF,self.SiMF,self.MgMF,_,_,_ = chemistry(self.CMF,xSi=self.xSi,xFe=self.xFe,xWu=self.xWu)
    
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
            Rocky core mass fraction (rcmf = (1-wmf)/cmf)
        self.FeMF: array
            Iron mass fraction
        self.SiMF: array
            Silicon mass fraction
        self.MgMF: array
            Magnesium mass fraction
        '''
        get_R = lambda x: interpn(PointsWater, Radius_DataWater, x) # x=wmf,Mass,cmf   
        self._check_water(get_R)
        args = np.asarray([self.Radius,self.Mass,self.CMF]).T
        residual = lambda x,param: np.sum(param[0]-get_R(np.asarray([x[0],param[1],param[2]*(1-x[0])]).T))**2/1e-4 
        self.WMF = self._run_MC(residual,args)
        self.FeMF,self.SiMF,self.MgMF,_,_,_ = chemistry(self.CMF,xSi=0.,xFe=0.,xWu=0.2)

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
    rocky = True if wmf is None else False
    
    if isinstance(M, (list, np.ndarray)):
        n = len(M)
        wmf = np.full(n,wmf)
        cmf = np.full(n,cmf)
        xSi = np.full(n,xSi)
        xFe = np.full(n,xFe)
    
    if rocky:
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

def get_rhoe(M,R, **kwargs):
    '''
    Find the planet density normalized to Earth-like planet for the same mass.

    Parameters:
    -----------
    M: float or array
        Mass of the planet in Earth masses.
    R: float or array
        Radius of the planet in Earth radii.
    **kwargs: dict {'cmf': float, 'xSi': float, 'xFe': float}
        Optional parameters for radius calculation, default is Earth-like values.
    Returns:
    --------
    rhoe: float or array
        rho_bulk/rho_earth(M).
    '''
    if not kwargs:
        kwargs = {'cmf': 0.325, 'xSi': 0.2, 'xFe': 0.05}
    if isinstance(M, (list, np.ndarray)):
        rhoe = np.zeros(len(M))
        for i in range(len(M)):
            r_earth = get_radius(M[i],**kwargs)
            rhoe[i] = (r_earth/R[i])**3
    else:
        r_earth = get_radius(M,**kwargs)
        rhoe = (r_earth/R)**3
    return rhoe

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
        res = minimize(residual,1,args=args,bounds=[[10**-0.5,10**1.3]]).x[0]
    return res

def star_to_planet(Fe,Si,Mg,Ca=-2,Al=-2,Ni=-2,Sun=[7.46,7.51,7.55,6.20,6.30,6.43],
                   xSi=[0,0.2],xFe=[0,0.2],xCore_trace=0.02,tol=1e-8):
    '''
    Convert stellar abundances to planet abundances, using Monte Carlo sampling.
    Stellar abundances (X/H in dex) are assumed to be normalized to Asplund 2021 solar abundances.
    
    Parameters:
    -----------
    Fe, Mg, Si: array
        Stellar abundances of Fe, Si, Mg ratative to solar in dex.
    Ca, Al, Ni: array, optional
        Stellar abundances of Ca, Al, Ni ratative to solar in dex.
        Default is to set to -2 dex.
    Sun: list
        Solar abundances of Fe, Si, Mg, Ca, Al, Ni in log scale.
    xSi: list, array
        Range or array of silicon molar fraction in the core.
    xFe: list, array
        Range or array of iron molar fraction in the mantle.
    xCore_trace: float
        Molar fraction of trace metals in the core.
    tol: float
        Tolerance for the optimization.
    Returns:
    --------
    host_star: object
        Host star object with all the properties.
    '''
    if not isinstance(Fe, (np.ndarray)):
        Fe = np.array([Fe])
    N = len(Fe)
    # Assume Asplund 2021 solar abundances
    star = np.empty([6,N])
    mu = [55.85e-3,28.09e-3,24.31e-3,40.08e-3,26.98e-3,58.69e-3] # Fe, Si, Mg, Ca, Al, Ni atomic masses
    star_abundance = [Fe,Si,Mg,Ca,Al,Ni]
    for i in range(6):
        star[i] = 10**(star_abundance[i]+Sun[i])*mu[i]
    
    if len(xFe)==2:
        xfe = np.random.uniform(xFe[0],xFe[1],N)
    if len(xSi)==2:
        xsi = np.random.uniform(xSi[0],xSi[1],N)

    Fe2Si,Fe2Mg,Mg2Si = star[0]/star[1],star[0]/star[2],star[2]/star[1]
    # if Ni,Ca,Al are not provided, set to very low value so xCa, xAl, xNi are zero
    Ni2Fe,Ca2Mg,Al2Mg = star[5]/star[0],star[3]/star[2],star[4]/star[2]
    star_ratio = [Fe2Si,Fe2Mg,Mg2Si,Ni2Fe,Ca2Mg,Al2Mg]

    def residual(x,param):
        # residual function for Monte Carlo sampling
        cmf, Xmgsi, xNi, xAl, xCa = x
        Fe2Si,Fe2Mg,Mg2Si,Ni2Fe,Ca2Mg,Al2Mg,xFe,xSi = param
        
        femf,simf,mgmf,camf,almf,nimf = chemistry(cmf,xSi=xSi,xFe=xFe,trace_core=xCore_trace,
                                   xNi=xNi,xAl=xAl,xCa=xCa,xWu=0,xSiO2=0)
        xSiO2, xWu = (0, Xmgsi) if Mg2Si > mgmf / simf else (Xmgsi, 0)
        femf,simf,mgmf,camf,almf,nimf = chemistry(cmf,xSi=xSi,xFe=xFe,trace_core=xCore_trace,
                                   xNi=xNi,xAl=xAl,xCa=xCa,xWu=xWu,xSiO2=xSiO2)
        res = ( (femf/simf - Fe2Si)**2/1e-4 + (femf/mgmf - Fe2Mg)**2/1e-4 + (mgmf/simf - Mg2Si)**2/1e-4 +
                (nimf/femf - Ni2Fe)**2/1e-3 + (camf/mgmf - Ca2Mg)**2/1e-3 + (almf/mgmf - Al2Mg)**2/1e-3 )
        return res
    
    model_param = np.zeros([N,8])
    planet_data = np.zeros([N,6])
    for i in range(N):
        xFe,xSi = xfe[i],xsi[i]
        param = [Fe2Si[i],Fe2Mg[i],Mg2Si[i],Ni2Fe[i],Ca2Mg[i],Al2Mg[i],xFe,xSi]
        res = minimize(residual,[0.325,0.2,0,0,0],args=param,tol=tol,
                       bounds=[[1e-15,1-1e-15],[1e-15,0.5],[0,0.2],[0,0.2],[0,0.2]])
        
        if res.success:
            cmf, Xmgsi, xNi, xAl, xCa = res.x
            femf,simf,mgmf,nimf,camf,almf = chemistry(cmf,xSi=xSi,xFe=xFe,trace_core=xCore_trace,
                                   xNi=xNi,xAl=xAl,xCa=xCa,xWu=0,xSiO2=0)
            xSiO2, xWu = (0, Xmgsi) if Mg2Si[i] > mgmf / simf else (Xmgsi, 0)
            data = chemistry(cmf,xSi=xSi,xFe=xFe,trace_core=xCore_trace,
                                    xNi=xNi,xAl=xAl,xCa=xCa,xWu=xWu,xSiO2=xSiO2)
    
            model_param[i] = cmf,xSi,xFe,xNi,xAl,xCa,xWu,xSiO2
            planet_data[i] = data               
            
        else:
            model_param[i] = np.repeat(np.nan,8)
            planet_data[i] = np.repeat(np.nan,6)
    return host_star(star_abundance,star_ratio,planet_data,model_param)
    
