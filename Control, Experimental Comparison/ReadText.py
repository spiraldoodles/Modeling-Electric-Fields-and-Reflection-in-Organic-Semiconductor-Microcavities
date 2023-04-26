import csv
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
f2 = open("C:\\Users\\spira\\Desktop\\Research\\Control, Experimental Comparison\\TwentyAOI.txt", 'r') 
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

normc=np.max(y)
plt.plot(x,y/normc)

plt.title("Reflectance vs Wavelength, S-Polarized Light")    
plt.ylim(0,1)
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance (Normalized)")
plt.show()