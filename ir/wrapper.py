import Templates
import importlib
importlib.reload(Templates)
import numpy as np

# can update these from Table 2 of Petter et al. 20202
names =  [          'J0106','J0826','J0827','J0905','J0908','J0944','J1039','J1107','J1125','J1219','J1229','J1232','J1248','J1341','J1506','J1613','J2116','J2118','J2140','J2256']
redshifts = np.array([0.454,  0.603,  0.681,  0.711,  0.502,  0.514,  0.634,  0.467,  0.519,  0.451,  0.614,  0.401,  0.632,  0.661,  0.437,  0.449,  0.728,  0.459,  0.751,  0.727])


for i in range(0, len(names)):

    #test = Templates.IR_SFRs(0.454, 'J0106')
    print(names[i])
    test = Templates.IR_SFRs(redshifts[i], names[i])
    print(test)
