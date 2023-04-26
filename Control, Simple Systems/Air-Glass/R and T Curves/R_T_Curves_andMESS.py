import sys
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
from scipy.signal import argrelmin
plt.rcParams.update({'font.size':12})

h=6.62607015 * 10**-34
c=3*10**8*1000000 #c in micrometers/second

L = .100 #Width of TIPS-Tc
lam0 = 0.5876 #Wavelength

wls = np.arange(0.4,0.9,0.1) #Wavelengths
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

        ],
        injection_material = eml.catalog.dielectrics.BK7 # Can specify material for transmission, injection and gap layers. If not specified, set to air, air & vaccum, respectively.
).compile()
print(system)
system.solve(lam0,AOI*np.pi/180) #Solve for varying WLS and AOI
z = np.arange(-0.05, system.left_edges[-1]+0.05, 0.005)

###     Reflectance vs AOI Graphs               ###
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


###     E vs Aoi vs Depth, Contours   ###  

#P-polarized
plt.figure()
E_vals= np.zeros((len(aoi),len(z)), dtype="complex")
for a_index in range(0,len(aoi)):
    system.solve(lam0, aoi[a_index]*np.pi/180, save_field=True)
    E= system.get_field(z, "p")
    E_net=np.multiply(np.conjugate(E),E)
    E_net=np.sum(E_net, axis=2)
    E_net=E_net[:,0]
    E_vals[a_index,:]=E_net

Z, A = np.meshgrid(z, aoi)
cs= plt.contourf(Z, A, E_vals, cmap="viridis", levels=100)
plt.xlabel(r"Depth $(\mu m)$")
plt.ylabel(r"AOI")
plt.colorbar(cs, label="Electric Field Density")

#S-polarized
plt.figure()
E_vals= np.zeros((len(aoi),len(z)), dtype="complex")
for a_index in range(0,len(aoi)):
    system.solve(lam0, aoi[a_index]*np.pi/180, save_field=True)
    E= system.get_field(z, "s")
    E_net=np.multiply(np.conjugate(E),E)
    E_net=np.sum(E_net, axis=2)
    E_net=E_net[:,0]
    E_vals[a_index,:]=E_net

Z, A = np.meshgrid(z, aoi)
cs= plt.contourf(Z, A, E_vals, cmap="viridis", levels=100)
plt.xlabel(r"Depth $(\mu m)$")
plt.ylabel(r"AOI")
plt.colorbar(cs, label="Electric Field Density")

#Components
# plt.figure()
# E_valsx= np.zeros((len(aoi),len(z)), dtype="complex")
# for a_index in range(0,len(aoi)):
#     system.solve(lam0, aoi[a_index]*np.pi/180, save_field=True)
#     E= system.get_field(z, "s")
#     E_net=np.multiply(np.conjugate(E),E)
#     print(np.shape(E_net))
#     E_netx=E_net[:,0,0]
#     E_valsx[a_index,:]=E_netx

# Z, A = np.meshgrid(z, aoi)
# cs= plt.contourf(Z, A, E_valsx, cmap="viridis", levels=100)
# plt.xlabel(r"Depth $(\mu m)$")
# plt.ylabel(r"AOI")
# plt.colorbar(cs, label="Electric Field Density")

# plt.figure()
# E_valsy= np.zeros((len(aoi),len(z)), dtype="complex")
# for a_index in range(0,len(aoi)):
#     system.solve(lam0, aoi[a_index]*np.pi/180, save_field=True)
#     E= system.get_field(z, "s")
#     E_net=np.multiply(np.conjugate(E),E)
#     E_nety=E_net[:,0,1]
#     E_valsy[a_index,:]=E_nety

# Z, A = np.meshgrid(z, aoi)
# cs= plt.contourf(Z, A, E_valsy, cmap="viridis", levels=100)
# plt.xlabel(r"Depth $(\mu m)$")
# plt.ylabel(r"AOI")
# plt.colorbar(cs, label="Electric Field Density")

# plt.figure()
# E_valsz= np.zeros((len(aoi),len(z)), dtype="complex")
# for a_index in range(0,len(aoi)):
#     system.solve(lam0, aoi[a_index]*np.pi/180, save_field=True)
#     E= system.get_field(z, "s")
#     E_net=np.multiply(np.conjugate(E),E)
#     E_netz=E_net[:,0,2]
#     E_valsz[a_index,:]=E_netz

# Z, A = np.meshgrid(z, aoi)
# cs= plt.contourf(Z, A, E_valsz, cmap="viridis", levels=100)
# plt.xlabel(r"Depth $(\mu m)$")
# plt.ylabel(r"AOI")
# plt.colorbar(cs, label="Electric Field Density")

plt.show()