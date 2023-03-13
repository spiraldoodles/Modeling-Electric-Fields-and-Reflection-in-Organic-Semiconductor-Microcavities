from multiprocessing.resource_sharer import stop
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
import math
from scipy.interpolate import UnivariateSpline
plt.rcParams.update({'font.size':20})

###     This code is for comparison with experimental data of DBR Mirror        ###


lam0 = .582 #Wavelength

aoi = np.array([20, 40, 60]) #Angles of Incidence
wls = np.arange(0.4,0.901,0.001) #Wavelengths

WLS = np.hstack([wls for _ in aoi])
AOI = np.hstack([len(wls)*[_] for _ in aoi])

#Build the system:
L1_Material=eml.catalog.dielectrics.TiO2
L1_Depth= 0.0579 #57.9 nm
L2_Material=eml.catalog.dielectrics.SiO2
L2_Depth= 0.0916 #91.6 nm

system = eml.layers.System(
        [
                eml.layers.Layer(L1_Material, L1_Depth), #10.5 Pairs of TiO2 and SiO2 to make DBR Mirror
                eml.layers.Layer(L2_Material, L2_Depth),
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
                eml.layers.Layer(L1_Material, L1_Depth),
                eml.layers.Layer(L2_Material, L2_Depth),
                eml.layers.Layer(L1_Material, L1_Depth),
                eml.layers.Layer(L2_Material, L2_Depth),
                eml.layers.Layer(L1_Material, L1_Depth),
                eml.layers.Layer(L2_Material, L2_Depth),
                eml.layers.Layer(L1_Material, L1_Depth),
                eml.layers.Layer(L2_Material, L2_Depth),
                eml.layers.Layer(L1_Material, L1_Depth)
        ],
        transmission_material = eml.catalog.dielectrics.BK7 # Can specify material for transmission, injection and gap layers. If not specified, set to air, air & vaccum, respectively.
).compile()

###     Reflectance vs Wavelength Graphs        ###
###     Constant lamda, varying aoi             ###

system.solve(WLS,AOI*np.pi/180) #Solve for varying WLS and AOI

R,T = system.get_RT("s") #Get reflection and transmission coefficent, s-polarized light
R = R.reshape(len(aoi),len(wls))
T = T.reshape(len(aoi),len(wls))

#AOI= 20
plt.figure()
normt=np.max(R[0])
plt.plot(wls,R[0]/normt,label="Theoretical")

f2 = open("C:\\Users\\spira\\OneDrive\\Desktop\\Research\\Control, Experimental Comparison\\TwentyAOI.txt", 'r') 
lines = f2.readlines()
x = []
y = [] 
for line in lines:
    p = line.split()
    x.append(float(p[0]))
    y.append(float(p[1]))
f2.close() 

x= np.array(x)
y= np.array(y)
print(x)
print(y)

norme=np.max(y)
plt.plot(x/1000,y/norme, label="Experimental")

# plt.title("Reflectance vs Wavelength, S-Polarized Light" )
plt.legend(loc='upper right')
plt.ylim(0,1)
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance")

#AOI= 40
plt.figure()
normt=np.max(R[1])
plt.plot(wls,R[1]/normt,label="Theoretical")

f2 = open("C:\\Users\\spira\\OneDrive\\Desktop\\Research\\Control, Experimental Comparison\\FourtyAOI.txt", 'r') 
lines = f2.readlines()
x = []
y = [] 
for line in lines:
    p = line.split()
    x.append(float(p[0]))
    y.append(float(p[1]))
f2.close() 

x= np.array(x)
y= np.array(y)
print(x)
print(y)

norme=np.max(y)
plt.plot(x/1000,y/norme, label="Experimental")

# plt.title("Reflectance vs Wavelength, S-Polarized Light" )
plt.legend(loc='upper right')
plt.ylim(0,1)
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance")

#AOI= 60
plt.figure()
normt=np.max(R[2])
plt.plot(wls,R[2]/normt,label="Theoretical")

f2 = open("C:\\Users\\spira\\OneDrive\\Desktop\\Research\\Control, Experimental Comparison\\SixtyAOI.txt", 'r') 
lines = f2.readlines()
x = []
y = [] 
for line in lines:
    p = line.split()
    x.append(float(p[0]))
    y.append(float(p[1]))
f2.close() 

x= np.array(x)
y= np.array(y)
print(x)
print(y)

norme=np.max(y)
plt.plot(x/1000,y/norme, label="Experimental")

# plt.title("Reflectance vs Wavelength, S-Polarized Light" )
plt.legend(loc='upper right')
plt.ylim(0,1)
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance")

plt.show()

