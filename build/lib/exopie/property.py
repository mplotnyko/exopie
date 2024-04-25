import os
import pickle
import numpy as np
from scipy.optimize import minimize


class planet_property:
    def __init__(self, N=50000, Mass=None, Radius=None, CMF=None, FeMF=None, WMF=None, xSi=None, xFe=None, atm_height=None):
        self._N = N
        self._Mass = np.array(Mass)
        self._Radius = np.array(Radius)
        self._xSi = np.array(xSi)
        self._xFe = np.array(xFe)
        self._CMF = np.array(CMF)
        self._FeMF = np.array(FeMF)
        self._WMF = np.array(WMF)
        self._atm_height = np.array(atm_height)
    

    @property
    def N(self):
        return self._N
    @property
    def Mass(self):
        return self._Mass
    @Mass.setter
    def Mass(self, value):
        self._Mass = value
    @property
    def Radius(self):
        return self._Radius
    @Radius.setter
    def Radius(self, value):
        self._Radius = value
    @property
    def xSi(self):
        return self._xSi
    @xSi.setter
    def xSi(self, value):
        self._xSi = value
    @property
    def xFe(self):
        return self._xFe
    @xFe.setter
    def xFe(self, value):
        self._xFe = value
    @property
    def CMF(self):
        return self._CMF
    @CMF.setter
    def CMF(self, value):
        self._CMF = value
    @property
    def FeMF(self):
        return self._FeMF
    @FeMF.setter
    def FeMF(self, value):
        self._FeMF = value
    @property
    def WMF(self):
        return self._WMF
    @WMF.setter
    def WMF(self, value):
        self._WMF = value
    @property
    def atm_height(self):
        return self._atm_height
    @atm_height.setter
    def atm_height(self, value):
        self._atm_height = value

class exoplanet(planet_property):
    '''
    Find interior structure errors (cmf, Fe-mf) for a given exoplanet. 
    This code uses data from an interior structure model SUPEREARTH 
    deveopled by Valencia et al. 2006 and updated in Plotnykov & Valencia 2020.
    Interpolating the data, we construct a fit to find radius given mass and 
    interior parameters (cmf, wmf, xSi, xFe) of a planet.
    The methodoly is described in Plotnykov & Valencia 2024.

    Parameters
    ----------
    Mass: real
        Mass of the planet in Earth masses. Valid only in range [0.1,10].
    Radius: real
        Radius of the planet in Earth radii.
    Mass_error: list
        Mass error of the planet in either [upper,lower] or [error] format.
    Radius_error: list
        Radius error of the planet in either [upper,lower] or [error] format.
    planet: string
        Planet type, either 'rocky' (find cmf) or 'ocean' (find wmf) planet.
    '''

    def __init__(self, N=50000, Mass=[1,0.001], Radius=[1,0.001], **kwargs):
        N = len(Mass) if len(Mass)>3 else int(N)
        N = len(Radius) if len(Radius)>3 else int(N)      
        super().__init__(N, Mass, Radius, **kwargs)
        if len(Mass)<=3:
            self.set_Mass(mu=Mass[0], sigma=Mass[1:])
        if len(Radius)<=3:
            self.set_Radius(mu=Radius[0], sigma=Radius[1:])

    def set_Mass(self, mu=1, sigma=0.001):
        self.Mass = self._set_parameter(mu,sigma)

    def set_Radius(self, mu=1, sigma=0.001):
        self.Radius = self._set_parameter(mu,sigma)
    
    def set_xSi(self, a=0, b=0.2):
        self.xSi = np.random.uniform(a,b,self.N)
    
    def set_xFe(self, a=0, b=0.2):
        self.xFe = np.random.uniform(a,b,self.N)
    
    def set_CMF(self, a=0, b=1):
        self.CMF = np.random.uniform(a,b,self.N)

    def set_atm_height(self, a=20, b=30):
        self.atm_height = np.random.uniform(a,b,self.N)

    def _set_parameter(self, mu, sigma):
        if type(sigma) == np.ndarray or type(sigma) == list:
            try:
                return np.random.choice(self._skewposterior(mu,sigma[0],sigma[1],self.N),self.N)
            except:
                return np.random.normal(mu, sigma[0], self.N)    
        else:
            return np.random.normal(mu, sigma, self.N)
    
    def _skewposterior(self, mu, sigma_up, sigma_lw, N):
        UP = np.random.normal(0,abs(sigma_up),size=N)
        LW = np.random.normal(0,abs(sigma_lw),size=N)
        return mu + np.concatenate([UP[UP>0],LW[LW<0]])
    
    def _check_rocky(self,get_R):
        # check if parameters are in bounds
        if (self.Mass==None).any():
            raise Exception('Mass must be set.')
        elif (self.Radius==None).any():
            raise Exception('Radius must be set.')
        if len(self.Mass) != len(self.Radius):
            raise Exception('Mass and Radius must be of same length.')

        pos = (self.Mass>10**-0.5) & (self.Mass<10**1.3) & (0<=self.xSi) & (self.xSi<=0.2) & (0<=self.xFe) & (self.xFe<=0.2)
        for item in  ['Mass','Radius','xSi','xFe']:
            setattr(self, item, getattr(self, item)[pos])

        pos = ( (self.Radius > get_R(np.asarray([np.repeat(1,sum(pos)),self.Mass,self.xSi,self.xFe]).T)) &
                (self.Radius < get_R(np.asarray([np.repeat(0,sum(pos)),self.Mass,self.xSi,self.xFe]).T)))
        self._N = sum(pos)
        if sum(pos)==0:
            raise Exception('Wrong planet type, no M-R pair in bounds')
        for item in  ['Mass','Radius','xSi','xFe']:
            setattr(self, item, getattr(self, item)[pos])        
        
    def _check_ocean(self,get_R):
        pos = (self.Mass>10**-0.5) & (self.Mass<10**1.3) & (0<=self.CMF) & (self.CMF<=1)
        for item in  ['Mass','Radius','CMF']:
            setattr(self, item, getattr(self, item)[pos])
            
        pos = ( (self.Radius < get_R(np.asarray([np.repeat(1,self.N),self.Mass,np.repeat(0,self.N)]).T)) &
                (self.Radius > get_R(np.asarray([np.repeat(0,self.N),self.Mass,self.CMF]).T))) 
        self._N = sum(pos)
        if self.N==0:
            raise Exception('Wrong planet type, no M-R pair in bounds')
        for item in  ['Mass','Radius','CMF']:
            setattr(self, item, getattr(self, item)[pos])

    def _run_MC(self,residual,args,bounds=[[0,1]],xi=0.3):
        res = []
        for i in range(self.N):
            res.append(minimize(residual,xi,args=args[i],bounds=bounds).x)
        if not isinstance(xi, (list, np.ndarray)):
            return np.asarray(res).flatten()
        return np.asarray(res).reshape(-1,self.N).T
    
    def save_data(self,filename=None):
        if filename is None:
            filename = 'data_run_0.pkl'
            i=0
            while os.path.exists(filename):
                i+=1
                filename = f'data_run_{i}.pkl'
        with open(filename,'wb') as f:
            dic = {}
            for item in  ['Mass','Radius','CMF','FeMF','WMF','xSi','xFe','atm_height']:
                dic[item] = getattr(self, item)
            pickle.dump(dic,f)

def load_Data():
    package_dir = os.path.dirname(__file__)
    # load rocky data
    with open(package_dir+'/Data/MRdata_rocky.pkl','rb') as f:
        Data = pickle.load(f)
        PointsRocky = [Data['CMF'],Data['Mass'],Data['xSi'],Data['xFe']]
        Radius_DataRocky = Data['Radius_total'] # tuple of radius data in Re
    # load water data
    with open(package_dir+'/Data/MRdata_water.pkl','rb') as f:
        Data = pickle.load(f)
        PointsWater = [Data['WMF'],Data['Mass'],Data['CMF']]
        Radius_DataWater = Data['Radius_total'] # tuple of radius data in Re
    return PointsRocky,Radius_DataRocky,PointsWater,Radius_DataWater