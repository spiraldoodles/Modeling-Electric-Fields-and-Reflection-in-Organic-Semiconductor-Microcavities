import sys
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
from scipy.signal import argrelmin
plt.rcParams.update({'font.size':12})

h=6.62607015 * 10**-34
c=3*10**8*1000000 #c in micrometers/second

L = .100 #Width of TIPS-Tc
lam0 = 0.532 #Wavelength

wls = np.arange(lam0,lam0+1,1) #Wavelengths
aoi = np.arange(0,90,1)#Angles of Incidence
spec_wls = np.arange(0.4,0.65,0.05)
spec_wls= np.append(spec_wls, lam0)

WLS = np.hstack([wls for _ in aoi])
AOI = np.hstack([len(wls)*[_] for _ in aoi])

#Build the system:
L1_Material=eml.catalog.dielectrics.TiO2
L1_Depth= lam0/(4*L1_Material.n(lam0)[0]) #lambda/4n - n depends on lambda and material
L2_Material=eml.catalog.dielectrics.SiO2
L2_Depth= lam0/(4*L2_Material.n(lam0)[0])#lambda/4n - n depends on lambda and material
L3_Material=eml.catalog.dielectrics.TIPSTc
L3_Depth= L

system = eml.layers.System(
        [
            eml.layers.Layer(L1_Material, L1_Depth),
            eml.layers.Layer(L2_Material, L2_Depth)
        ],
        transmission_material = eml.catalog.dielectrics.BK7 # Can specify material for transmission, injection and gap layers. If not specified, set to air, air & vaccum, respectively.
).compile()
print(system)
system.solve(WLS,AOI*np.pi/180) #Solve for varying WLS and AOI


###     Reflectance vs AOI Graphs        ###
###     Constant lamda, varying aoi             ###

plt.figure()
#P-Polarized
R,T = system.get_RT("p") #Get reflection and transmission coefficent, p-polarized light
R = R.reshape(len(aoi),len(wls))
T = T.reshape(len(aoi),len(wls))

plt.plot(aoi, R[:,0], label="p-polarized")

#S-Polarized
R,T = system.get_RT("s") #Get reflection and transmission coefficent, s-polarized light
R = R.reshape(len(aoi),len(wls))
T = T.reshape(len(aoi),len(wls))


plt.plot(aoi, R[:,0], label="s-polarized")

plt.ylim(0,1)
plt.xlabel(r"AOI")
plt.ylabel("Reflectance")
plt.legend()

###     Transmission vs AOI Graphs        ###
###     Constant lamda, varying aoi             ###
# plt.figure()
# #P-Polarized
# R,T = system.get_RT("p") #Get reflection and transmission coefficent, p-polarized light
# R = R.reshape(len(aoi),len(wls))
# T = T.reshape(len(aoi),len(wls))

# normc=np.max(T[:,0])
# print(normc)
# plt.plot(aoi, T[:,0], label="p-polarized")

# #S-Polarized
# R,T = system.get_RT("s") #Get reflection and transmission coefficent, s-polarized light
# R = R.reshape(len(aoi),len(wls))
# T = T.reshape(len(aoi),len(wls))

# normc=np.max(T[:,0])
# print(normc)
# plt.plot(aoi, T[:,0], label="s-polarized")

# plt.ylim(0,1)
# plt.xlabel(r"AOI")
# plt.ylabel("Transmittance")
# plt.legend()



plt.show()