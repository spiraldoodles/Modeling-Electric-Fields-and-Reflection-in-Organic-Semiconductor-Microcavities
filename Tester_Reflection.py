import sys
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
from scipy.signal import argrelmin
plt.rcParams.update({'font.size':20})

h=6.62607015 * 10**-34
c=3*10**8*1000000 #c in micrometers/second

L = .100 #Width of TIPS-Tc
lam0 = .532 #Wavelength

wls = np.arange(0.4,0.9,0.001) #Wavelengths
aoi = np.arange(0,81,20)#Angles of Incidence
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
            eml.layers.Layer(L2_Material, L2_Depth),

            eml.layers.Layer(L1_Material, L1_Depth),
            eml.layers.Layer(L2_Material, L2_Depth),

            eml.layers.Layer(L1_Material, L1_Depth),
            eml.layers.Layer(L2_Material, L2_Depth),

            eml.layers.Layer(L1_Material, L1_Depth),
            eml.layers.Layer(L2_Material, L2_Depth),

            eml.layers.Layer(L1_Material, L1_Depth),
            eml.layers.Layer(L2_Material, L2_Depth),

            # eml.layers.Layer(L1_Material, L1_Depth),
            # eml.layers.Layer(L2_Material, L2_Depth),

            # eml.layers.Layer(L1_Material, L1_Depth),
            # eml.layers.Layer(L2_Material, L2_Depth),

            # eml.layers.Layer(L1_Material, L1_Depth),
            # eml.layers.Layer(L2_Material, L2_Depth),

            # eml.layers.Layer(L1_Material, L1_Depth),
            # eml.layers.Layer(L2_Material, L2_Depth),

            # eml.layers.Layer(L1_Material, L1_Depth),
            # eml.layers.Layer(L2_Material, L2_Depth),            

            eml.layers.Layer(L1_Material, L1_Depth),
        ],
        transmission_material = eml.catalog.dielectrics.BK7 # Can specify material for transmission, injection and gap layers. If not specified, set to air, air & vaccum, respectively.
).compile()
print(system)
z = np.arange(-0.05, system.left_edges[-1]+0.05, 0.01) #Depth position in z direction along stack, with 0 at the left edge of the silver mirror
system.solve(WLS,AOI*np.pi/180) #Solve for varying WLS and AOI


###     Reflectance vs Wavelength Graphs        ###
###     Constant lamda, varying aoi             ###


#P-Polarized
# R,T = system.get_RT("p") #Get reflection and transmission coefficent, p-polarized light
# R = R.reshape(len(aoi),len(wls))
# T = T.reshape(len(aoi),len(wls))

# plt.figure()
# for aindex in range(len(aoi)):
#     normc=np.max(R[aindex])
#     plt.plot(wls,R[aindex,:]/normc,label=r"$%.0f\degree$" % aoi[aindex])
#     # plt.plot(wls,R[aindex,:],label=r"$%.0f\degree$" % aoi[aindex])

# plt.legend(title="AOI", loc='upper right')
# plt.ylim(0,1)
# plt.xlabel(r"Wavelength $(\mu m)$")
# plt.ylabel("Reflectance (Normalized)")

# #S-Polarized
# R,T = system.get_RT("s") #Get reflection and transmission coefficent, s-polarized light
# R = R.reshape(len(aoi),len(wls))
# T = T.reshape(len(aoi),len(wls))

# plt.figure()
# for aindex in range(len(aoi)):
#     normc=np.max(R[aindex])
#     plt.plot(wls,R[aindex,:]/normc,label=r"$%.0f\degree$" % aoi[aindex])
#     # plt.plot(wls,R[aindex,:],label=r"$%.0f\degree$" % aoi[aindex])
  
# plt.legend(title="AOI", loc='upper right')
# plt.ylim(0,1)
# plt.xlabel(r"Wavelength $(\mu m)$")
# plt.ylabel("Reflectance (Normalized)")


###     Reflectance vs Wavelength Graphs SHIFTED    ###
###     Constant lamda, varying aoi                 ###

#P-Polarized
R,T = system.get_RT("p") #Get reflection and transmission coefficent, p-polarized light
R = R.reshape(len(aoi),len(wls))
T = T.reshape(len(aoi),len(wls))

plt.figure()
for aindex in range(len(aoi)):
    normc=np.max(R[aindex])
    plt.plot(wls,R[aindex,:]/normc+aindex,label=r"$%.0f\degree$" % aoi[aindex])
plt.axvline(x=lam0, color= "grey", linestyle='dashed')

plt.legend(title="AOI", loc='upper right')
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance (Shifted)")
plt.tick_params(left = False, labelleft = False)

#S-Polarized
R,T = system.get_RT("s") #Get reflection and transmission coefficent, s-polarized light
R = R.reshape(len(aoi),len(wls))
T = T.reshape(len(aoi),len(wls))

plt.figure()
for aindex in range(len(aoi)):
    normc=np.max(R[aindex])
    plt.plot(wls,R[aindex,:]/normc+aindex,label=r"$%.0f\degree$" % aoi[aindex])
plt.axvline(x=lam0, color= "grey", linestyle='dashed')
  
plt.legend(title="AOI")
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance (Shifted)")
plt.tick_params(left = False, labelleft = False)

# ###     Reflectance Dips as a function of Angle        ###
# ###     Constant lamda, varying aoi                    ###

# #P-Polarized
# R,T = system.get_RT("p") #Get reflection and transmission coefficent, p-polarized light
# R = R.reshape(len(aoi),len(wls))
# T = T.reshape(len(aoi),len(wls))

# plt.figure()
# for aindex in range(len(aoi)):
#     mins=argrelmin(R[aindex,:])
#     mins= np.transpose(mins)
#     E_mins=np.zeros(len(mins))
#     for min_index in range(len(mins)):
#         E_mins[min_index]=h*c/wls[mins[min_index]]*6.2415*10**18 
#     anglearray=np.full(len(mins),aoi[aindex])
#     plt.scatter(anglearray, E_mins)

# plt.xlabel("AOI")
# plt.ylabel(r"Energy (eV)")

# #S-Polarized
# R,T = system.get_RT("s") #Get reflection and transmission coefficent, p-polarized light
# R = R.reshape(len(aoi),len(wls))
# T = T.reshape(len(aoi),len(wls))

# plt.figure()
# for aindex in range(len(aoi)):
#     mins=argrelmin(R[aindex,:])
#     mins= np.transpose(mins)
#     E_mins=np.zeros(len(mins))
#     for min_index in range(len(mins)):
#         E_mins[min_index]=h*c/wls[mins[min_index]]*6.2415*10**18 
#     anglearray=np.full(len(mins),aoi[aindex])
#     plt.scatter(anglearray, E_mins)

# plt.xlabel("AOI")
# plt.ylabel(r"Energy (eV)")

plt.show()