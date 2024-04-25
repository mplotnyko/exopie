import os
import numpy as np
from scipy.interpolate import interpn
import warnings
import pickle

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

def radius(M,cmf=0.325,wmf=0,xSi=0,xFe=0.1):
    '''
    Find the Radius of a planet, given mass and interior parameters.
    Mass (M) has to be in Earth masses and a float.
    '''
    if wmf>0:
        if cmf+wmf>1:
            warnings.warn('CMF+WMF>1, decrease cmf or wmf. Using cmf=1-wmf')
            cmf=1-wmf
        return interpn(_PointsWater, _Radius_DataWater, [wmf,M,cmf])[0]
    return interpn(_PointsRocky, _Radius_DataRocky, [cmf,M,xSi,xFe])[0]

def sigmacmf(M,R,dM,dR):
    '''
    Analytical function to find ..
    Eq 2 in the paper Pl & Valencia 2023
    '''
    M = M*5.97219e24
    R = R*6.371e6
    rho_bulk = M/(4/3*np.pi*R**3)/1e3
    # find_rho = lambda rho:  (cmf - (rho[1]/rho[0] - (rho[1]/rho_bulk)) / (rho[1]/rho[0] - 1))**2
    # res = optimize.minimize(find_rho,[3,10],bounds=[[3,8],[7,30]]) # rho = rho_m,rho_c
    # rho_m,rho_c = res.x
    rho_m,rho_c = 5,10
    return (rho_c/rho_bulk)*np.sqrt(9*dR**2+dM**2)/(rho_c / rho_m - 1)

def deltacmf(M,R,dM,dR):
    '''
    Analytical function to find ..
    Eq 2 in the paper Pl & Valencia 2023
    '''
    rho_m,rho_c = 5,10
    rho_bulk = M*5.97219e24/(4/3*np.pi*(R*6.371e6)**3)/1e3 
    return (rho_c/rho_bulk)*(dM/M-3*dR/R)/(rho_c / rho_m - 1)

def chemistry(cmf,xSi=0,xFe=0.1,trace_core=0.02,
                xNi=None,xAl=0.04,xCa=0.05,xWu=0.2,wmf=0):
    '''
    
    '''
    mmf = 1-cmf-wmf # mantle mass fraction
    xPv = 1-xWu-xCa #molar fraction of porovskite in the mantle
    Fe,Ni,Si,Mg,Ca,Al,O,XCore = [55.85e-3,58.69e-3,28.09e-3,24.31e-3,
                                    40.078e-3,26.98e-3,16e-3,50e-3] # atmoic masses, Xcore stands for other metals in the core
    if xNi is None:
        xNi = (1-xSi-trace_core)/(16*Ni/Fe+1) #based on McDonough&Sun 1995 Fe/Ni ratio is 16 [w]
    xfe_core = 1-xSi-xNi-trace_core #molar fraction of Fe in core
    core_mol = xfe_core*Fe+Si*xSi+Ni*xNi+trace_core*XCore
    man_mol = ( ((1-xFe)*Mg+xFe*Fe+O)*xWu + xCa*(Ca+Si+O*3) 
               +((1-xFe-xAl)*Mg+xFe*Fe+Si*(1-xAl)+xAl*Al*2+O*3)*xPv ) # assume lower mantle scpecies only

    fe_core = cmf*xfe_core*Fe/core_mol
    fe_man = mmf*(xPv+xWu)*xFe*Fe/man_mol 
    
    si_core = cmf*xSi*Si/core_mol
    si_man = mmf*(xPv*(1-xAl)+xCa)*Si/man_mol
    
    fe_mass = fe_core+fe_man
    si_mass = si_core+si_man
    mg_mass = mmf*(xPv*(1-xFe-xAl)+xWu*(1-xFe))*Mg/man_mol
    return fe_mass,si_mass,mg_mass
    
def magnisium_number(xFe,xWu,xCa,xAl):
    xPv = 1-xWu-xCa # calculate the mole fraction of Pv
    Fe_m = xPv*xFe+xWu*xFe #moles of Fe
    Mg_m = xPv*(1-xFe-xAl)+xWu*(1-xFe) #moles of Mg
    return Mg_m/(Fe_m+Mg_m)
    
_PointsRocky, _Radius_DataRocky, _PointsWater, _Radius_DataWater = load_Data() # load grid data for fitting interior

