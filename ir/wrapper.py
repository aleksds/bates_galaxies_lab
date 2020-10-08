import Templates
import importlib
importlib.reload(Templates)
import numpy as np

# can update these from Table 2 of Petter et al. 20202
names =  [          'J0106','J0826','J0827','J0905','J0908']
redshifts = np.array([0.454,  0.603,  0.681,  0.711,  0.502])


for i in range(0, len(names)):

    #test = Templates.IR_SFRs(0.454, 'J0106')
    print(names[i])
    test = Templates.IR_SFRs(redshifts[i], names[i])
    print(test)
