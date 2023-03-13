import sys
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
from scipy.signal import argrelmin
plt.rcParams.update({'font.size':20})

h=6.62607015 * 10**-34
c=3*10**8*1000000 #c in micrometers/second

L = .100 #Width of TIPS-Tc
lam0 = 0.532 #Wavelength


spec_aoi = np.arange(0,81,20)
spec_wls = np.arange(0.5,0.6,0.01)
spec_wls= np.append(spec_wls, lam0)


#Build the system:
L1_Material=eml.catalog.dielectrics.TiO2
L1_Depth= lam0/(4*L1_Material.n(0.532)[0]) #lambda/4n - n depends on lambda and material
L2_Material=eml.catalog.dielectrics.SiO2
L2_Depth= lam0/(4*L2_Material.n(0.532)[0])#lambda/4n - n depends on lambda and material
L3_Material=eml.catalog.dielectrics.TIPSTc
L3_Depth= L

system = eml.layers.System(
        [
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

            # eml.layers.Layer(L1_Material, L1_Depth),
        ],
        transmission_material = eml.catalog.dielectrics.BK7 # Can specify material for transmission, injection and gap layers. If not specified, set to air, air & vaccum, respectively.
).compile()
print(system)
print("TiO2, n="+str(L1_Material.n(lam0)[0]))
print("SiO2, n="+str(L2_Material.n(lam0)[0]))
print("Air, n="+str(eml.catalog.dielectrics.Air.n(lam0)[0]))
print("Glass, n="+str(eml.catalog.dielectrics.BK7.n(lam0)[0]))

z = np.arange(-0.05, system.left_edges[-1]+0.05, 0.001)

# Color-coding for easy reading
def colorcode():
    for numlayers in range(len(system)):
        if system[numlayers].mat == eml.catalog.metals.Ag:
            plt.axvspan(system.left_edges[numlayers-1], system.left_edges[numlayers], color = "grey", alpha = 0.2)
        elif system[numlayers].mat == eml.catalog.dielectrics.TiO2:
            plt.axvspan(system.left_edges[numlayers -1], system.left_edges[numlayers], color = "blue", alpha = 0.2)
        elif system[numlayers].mat == eml.catalog.dielectrics.SiO2:
            plt.axvspan(system.left_edges[numlayers -1], system.left_edges[numlayers], color = "purple", alpha = 0.2)
        elif system[numlayers].mat == eml.catalog.dielectrics.TIPSTc:
            plt.axvspan(system.left_edges[numlayers -1], system.left_edges[numlayers], color = "green", alpha = 0.2)
        elif system[numlayers].mat == system.trn.mat:
            plt.axvspan(system.left_edges[-1], system.left_edges[-1]+0.1, color = "yellow", alpha = 0.2)
        elif system[numlayers].mat == system.inj.mat:
            plt.axvspan(system.left_edges[0]-0.1, system.left_edges[0], color = "red", alpha = 0.2)
        else:
            print("uh oh")


###     E Field Density vs AOI                  ###
###     Constant lamda, varying aoi             ###   

#P-pol
plt.figure()
for a in spec_aoi:
    system.solve(lam0, a*np.pi/180, save_field=True)
    E= system.get_field(z, "p")
    E_net=np.multiply(np.conjugate(E),E)
    E_net=np.sum(E_net, axis=2)
    plt.plot(z, E_net, label= r"$%.0f\degree$" % a)

#Plot setup
plt.legend(title="AOI", loc='upper right')
plt.xlabel(r"Depth $(\mu m)$")
plt.ylabel(r"Electric Field Density")
plt.xlim(-0.05, system.left_edges[-1]+0.05)
colorcode()

#S-pol
plt.figure()
for a in spec_aoi:
    system.solve(lam0, a*np.pi/180, save_field=True)
    E= system.get_field(z, "s")
    E_net=np.multiply(np.conjugate(E),E)
    E_net=np.sum(E_net, axis=2)
    plt.plot(z, E_net, label= r"$%.0f\degree$" % a)

#Plot setup
plt.legend(title="AOI", loc='upper right')
plt.xlabel(r"Depth $(\mu m)$")
plt.ylabel(r"Electric Field Density")
plt.xlim(-0.05, system.left_edges[-1]+0.05)
colorcode()


###     E Field Density vs Wavelength                  ###
###     Constant aoi, varying lambda             ###   

#P pol 
# plt.figure()
# for w in spec_wls:
#     system.solve(w, 60*np.pi/180, save_field=True)
#     E= system.get_field(z, "p")
#     E_net=np.multiply(np.conjugate(E),E)
#     E_net=np.sum(E_net, axis=2)
#     plt.plot(z, E_net, label= r"$\lambda$" + "=" + str(round(w,3)))

# #Plot setup
# plt.legend(title="Wavelength")
# plt.xlabel(r"Depth $(\mu m)$")
# plt.ylabel(r"Electric Field Density")
# plt.xlim(-0.05, system.left_edges[-1]+0.05)
# colorcode()

plt.show()