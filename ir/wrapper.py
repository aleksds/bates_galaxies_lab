import Templates
import importlib
importlib.reload(Templates)
import numpy as np

# can update these from Table 2 of Petter et al. 20202
names = ['J0106', 'J0826']
redshifts = np.array([0.454, 0.603])


for i in range(0, len(names)):

    #test = Templates.IR_SFRs(0.454, 'J0106')
    test = Templates.IR_SFRs(redshifts[i], names[i])
    print(test)
