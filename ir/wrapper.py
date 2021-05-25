import Templates
from Templates import templates
from brown_seds_clean import files
import importlib
import glob
import os
importlib.reload(Templates)
import numpy as np



projpath = os.getcwd()+'/Brown2014/An_Atlas_of_Galaxy_SEDs/An_Atlas_of_Galaxy_SEDs/Converted/'
br_tems = glob.glob(projpath + '*.dat')
print("br_tems",br_tems)





#print("This is wrapper")
# can update these from Table 2 of Petter et al. 2020
#added J1558 for kicks
names = ['J0106', 'J0826', 'J0827', 'J0905', 'J0908', 'J0944', 'J1039', 'J1107', 'J1125', 'J1219', 'J1229', 'J1232',
         'J1248', 'J1341', 'J1506', 'J1613', 'J2116', 'J2118', 'J2140', 'J2256','J1558']
redshifts = np.array(
    [0.454, 0.603, 0.681, 0.711, 0.502, 0.514, 0.634, 0.467, 0.519, 0.451, 0.614, 0.401, 0.632, 0.661, 0.437, 0.449,
     0.728, 0.459, 0.751, 0.727,.403])

'''
for i in range(0, len(names)):

    #test = Templates.IR_SFRs(0.454, 'J0106')
    print(names[i])
    test = Templates.IR_SFRs(redshifts[i], names[i])
    print(test)

'''

# These are the 5 gals  for the JWST proposal
#jwst_redshifts = []  # We'll populate these redshifts using the for loop down later on
jwst_names = ['J1107', 'J1219', 'J1506', 'J1613', 'J2118']
jwst_redshifts = [0.467, 0.451, 0.608, 0.449, 0.459]
best_fits = [""]

poster_names = ['J1558']
poster_redshifts = [0.403]
# For special names
'''
for i in range(0,len(names)):
    if names[i] in special_names:
        print(names[i])
        #Call this one for best worst?
        sed = Templates.IR_SFRs(redshifts[i],names[i],'special_test',tems = special_temps[names[i]])

        #Call this one for all of them?
        #brown_sed = Templates.IR_SFRs(redshifts[i],names[i], tems = files)
        special_redshifts.append(redshifts[i])
        # print(sed)
        # print(brown_sed)

'''
print("files:", files)

def plot_seds(names_subset):
    redshift_list = []
    for i in range(0, len(names)):
        if names[i] in names_subset:
            print("name", names[i])
            print('redshifts[i]', redshifts[i])
            # IR_SFRs returns two variables, one avg SFR and the other Stdev SFR
            #What if tems doesnt equal files?
            #names_subset_sfr, names_subset_stdv_sfr, sed_list = Templates.IR_SFRs(redshifts[i], names[i], '_Analyze', tems=br_tems)  # Does this need to be an assignment statement?
            Templates.IR_SFRs(redshifts[i], names[i], '_poster',calc_SFR = True, tems=templates)
            #redshift_list.append(redshifts[i])
            #Templates.plot_all(names[i]+"All", sed_list)

            # print(sed)
            # print(brown_sed)

   # return redshift_list


#Write this one later to consolidate the messiness of each of the plots, choose which one to plot here
def plot_choice(choice):
    if choice == "JWST":
        jwst_names = ['J1107', 'J1219', 'J1506', 'J1613', 'J2118', 'J1558']
        jwst_redshifts = [0.467, 0.451, 0.608, 0.449, 0.459, 0.403]
        return jwst_names,jwst_redshifts

    elif choice == "Special":
        # These are for the original 8 galaxies of interest
        special_names = ['J0826', 'J0905', 'J0944', 'J1107', 'J1219', 'J1341', 'J1613', 'J2118']
        special_redshifts = []
        # This is for the best/worst SED pdf
        special_temps = {'J0826': [templates[1], templates[7]], 'J0905': [templates[1], templates[7]],
                         'J0944': [templates[1],
                                   templates[7]],
                         'J1107': [templates[1], templates[10]], 'J1219': [templates[1], templates[3]],
                         'J1341': [templates[1], templates[7]], 'J1613': [templates[0], templates[3]],
                         'J2118': [templates[1], templates[8]]}
        return special_names, special_redshifts,special_temps


    '''
    Choose which of the data set to plot, special 8, JWST, all gals etc.
    '''


#works for simplest case
"""
for i in range(0, len(names)):

    #test = Templates.IR_SFRs(0.454, 'J0106')
    print(names[i])
    test = Templates.IR_SFRs(redshifts[i], names[i],"simplest_case")
    print(test)

"""

#special_redshifts = plot_seds(special_names)

# For JWST names

#plot_seds('J1107')
plot_seds(poster_names)
"""
for i in range(0, len(names)):
    if names[i] in jwst_names:
        print("name", names[i])
        print('redshifts[i]', redshifts[i])
        print("jwst[i]", names[i])
        jwst_sed = Templates.IR_SFRs(redshifts[i], names[i], 'fool_test', tems=files)
        jwst_redshifts.append(redshifts[i])
        # print(sed)
        # print(brown_sed)
"""
#print("Special redshifts:", special_redshifts)




