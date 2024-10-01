def sigma_cmf(M,R,dM,dR):
    '''
    Analytical function to find the error in core mass fraction.
    '''
    M = M*5.97219e24
    R = R*6.371e6
    rho_bulk = M/(4/3*np.pi*R**3)/1e3
    rho_m,rho_c = 5,10
    return (rho_c/rho_bulk)*np.sqrt(9*dR**2+dM**2)/(rho_c / rho_m - 1)

def delta_cmf(M,R,dM,dR):
    '''
    Analytical function to find difference in core mass fraction.
    '''
    rho_m,rho_c = 5,10
    rho_bulk = M*5.97219e24/(4/3*np.pi*(R*6.371e6)**3)/1e3 
    return (rho_c/rho_bulk)*(dM/M-3*dR/R)/(rho_c / rho_m - 1)

def chemistry(cmf,xSi=0,xFe=0.1,trace_core=0.02,
                xNi=None,xAl=0.04,xCa=0.05,xWu=0.2):
    '''
    Calculate the interior chemistry of a planet.

    Parameters:
    ----------
    cmf: float or array
        Core mass fraction.
    xSi: float or array
        Molar fraction of silicon in the core.
    xFe: float or array
        Molar fraction of iron in the mantle.
    trace_core: float or array
        Molar fraction of trace metals in the core.
    xNi: float or array
        Molar fraction of nickel in the core. If None, assume Fe/Ni ratio is 16.
    xAl: float or array
        Molar fraction of aluminum in the mantle.
    xCa: float or array
        Molar fraction of calcium in the mantle.
    xWu: float or array
        Molar fraction of wustite in the mantle.
    
    Returns:
    --------
    FeMF, SiMF, MgMF: list
        Mass fractions of Fe, Si, and Mg.
    '''
    mmf = 1-cmf # mantle mass fraction
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
    '''
    Calculate the magnesium number in the mantle.

    Parameters:
    ----------
    xFe: float
        Molar fraction of iron in the mantle.
    xWu: float
        Molar fraction of wustite in the mantle.
    xCa: float
        Molar fraction of calcium in the mantle.
    xAl: float
        Molar fraction of aluminum in the mantle.
    '''
    xPv = 1-xWu-xCa # calculate the mole fraction of Pv
    Fe_m = xPv*xFe+xWu*xFe #moles of Fe
    Mg_m = xPv*(1-xFe-xAl)+xWu*(1-xFe) #moles of Mg
    return Mg_m/(Fe_m+Mg_m)