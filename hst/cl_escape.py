#calculates escape velocity from stellar masses obtained from prospector and radii obtained from galfit. does so for 475 and 814 using the best value prosepector gave us for stellar mass
import numpy as np
#need stellar mass in kg and radius in meters

G =6.674*(10**(-11)) #m^3 k^-1 s^-2

#getting stellar mass in kg
masslogsol = [10.48,10.58,10.12,10.61,10.54,11.09,10.29,10.23,10.98,11.35,10.67,10.11] #in standard galactic order
masssol = np.zeros(12)
masskg = np.zeros(12)
for w in range(0,12):
    masssol[w] = 10**(masslogsol[w])
for w in range(0,12):
    masskg[w] = masssol[w]*(1.98855*(10**30))

#getting radius in meters from pc
rem = np.zeros([12,2])
repc = [[84,134],[175,149],[86,83],[158,74],[62,179],[164,239],[170,81],[77,108],[215,1032],[560,886],[98,178],[91,189]]
for w in range(0,12):
    for i in range(0,2):
        rem[w][i] = (repc[w][i])*3.086*(10**(16))

#now plugging into escape velocity formula to get v_esc in km/s
esc = np.zeros([12,2])
for w in range(0,12):
    for i in range(0,2):
        esc[w][i] = (((2*G*masskg[w])/(rem[w][i]))**(1/2))/1000
