#Calculation of momentum input in dynes using dp/dt = 5*10^(33) * (SFR in solar masses per year)
import numpy as np
#SFR in order of J0826,J0901,J0905,J0944,J1107,J1219,J1341,J1506,J1558,J1613,J2116,J2140
SFR = [287,119,339,184,183,269,515,638,98,429,149,302]
momin = np.zeros(12)
for w in range(0,12):
    momin[w] = 5*10**(33)*SFR[w]
