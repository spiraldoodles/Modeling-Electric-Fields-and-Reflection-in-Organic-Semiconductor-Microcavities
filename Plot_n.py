import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
from scipy.signal import argrelmax
import sys
plt.rcParams.update({'font.size':20})
#Plot n vs lambda for matierlas used

wls = np.arange(0.4,0.9,0.001)

#SiO
# plt.figure()
n1=np.array(eml.catalog.dielectrics.SiO2.n(wls)[0])

# plt.plot(wls, n1, label="SiO2")
#TiO2
n2=np.array(eml.catalog.dielectrics.TiO2.n(wls)[0])

#Ag
n3=np.array(eml.catalog.metals.Ag.n(wls)[0])

#Air
n4=np.array(eml.catalog.dielectrics.Air.n(wls)[0])

#Glass
n5=np.array(eml.catalog.dielectrics.BK7.n(wls)[0])

#TIPS-Tc
ren= np.zeros(len(wls))
imn= np.zeros(len(wls))

for wl_index in range(len(wls)):
    ren[wl_index], imn[wl_index]=eml.catalog.dielectrics.TIPSTc.n(wls[wl_index])

# Plots
# plt.figure()

# plt.plot(wls, n1, label="SiO2")
# plt.plot(wls, n2, label="TiO2")
# plt.plot(wls, n3, label="Ag")
# plt.plot(wls, n4, label="Air")
# plt.plot(wls, n5, label="Glass")
# plt.plot(wls, ren, label="TIPS-Tc, Real")

# plt.axvline(x=.532, color= "grey", linestyle='dashed', label='532 nm')
# plt.xlabel("Wavelength $(\mu m)$")
# plt.ylabel("Refractive Index n")
# plt.ylim(0,2.75)
# plt.legend( loc="upper right")

#TIPS plots

plt.figure()
plt.plot(wls, ren)

plt.xlabel("Wavelength $(\mu m)$")
plt.ylabel("Refractive Index n")

plt.figure()
plt.plot(wls, imn)
print(type(imn))

maxs=argrelmax(imn)
maxs=np.transpose(maxs)
for m in maxs:
    plt.axvline(x=wls[m], color= "grey", linestyle='dashed', label= str(round(*wls[m], 3))+" nm")

plt.xlabel("Wavelength $(\mu m)$")
plt.ylabel("Attenuation $(\kappa)$")
plt.legend( loc="upper right")
plt.show()