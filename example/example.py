import exopie
import numpy as np

# find the rocky threshold radius (RTR, cmf=0) at 5 Earth masses
R = exopie.get_radius(5,xSi=0,xFe=0,cmf=0.)
print(f'Radius Threshold (at 5Me) = {R[0]:.2f} Re')

# Find the interior of M=5 +/-0.5, R=1.4 +/-0.05 planet
planet = exopie.rocky(N=50000,Mass=[5,0.5],Radius=[1.4,0.05]) #model input
planet.run() # run the model

print(f'CMF = {np.mean(planet.CMF):.2f}', 
      f'FeMF = {np.mean(planet.FeMF):.2f}')

# Case where no Si in the core (xSi=0) and no Fe in the mantle (xFe=0)
planet = exopie.rocky(N=50000,Mass=[5,0.5],Radius=[1.4,0.05],
                        xSi=[0,0],xFe=[0,0]) # model input
planet.run() # run the model
print(f'CMF = {np.mean(planet.CMF):.2f}', 
      f'FeMF = {np.mean(planet.FeMF):.2f}')

planet.save_data('rocky.pkl') # save the model


# To make a nice corner plot use (need corner.py)
# fig,axs = planet.corner()

