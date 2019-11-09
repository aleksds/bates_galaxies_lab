import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.cosmology import WMAP9 as cosmo
from scipy.constants import parsec as prsc
from scipy.constants import c as lightspeed
import astropy.units as u



# Importing the Umeh table
table = ascii.read('C:/Users/Chris/Documents/GitHub/bates_galaxies_lab/hst/umeh_table.dat')
umeh_names = np.array(table.field('Galaxy'))
umeh_w4mag = np.array(table.field('w4_mag'))
umeh_w4unc = np.array(table.field('w4_unc'))
# Redshifts for the galaxies in Umeh table taken from Tremonti table
umeh_Z_T = np.array([0.6031527, 0.45889184, 0.711421, 0.51383984, 0.46653742, 0.45092207, 0.66115844, 0.6080914,
                     0.402196, 0.44924328, 0.72798246, 0.7509297])
umeh_shortnames = np.array(['J0826+4305', 'J0901+0314', 'J0905+5759', 'J0944+0930', 'J1107+0417', 'J1219+0336',
                            'J1341-0321', 'J1506+5402', 'J1558+3957', 'J1613+2834', 'J2116-0634', 'J2140+1209'])

# Creating arrays with the Tremonti table
ChristyTablePath = 'files/hizea_wind_ancillary.fit'
ChristyTable = fits.open(ChristyTablePath)
ShortNames = ChristyTable[1].data['SHORT_NAME']
ChristyMass = ChristyTable[1].data['MASS']
ChristySFR = ChristyTable[1].data['SFR']
ChristyV50 = ChristyTable[1].data['VAVG']
ChristyV98 = ChristyTable[1].data['VMAX']
ChristyAge = ChristyTable[1].data['LW_AGE']
ChristyZ = ChristyTable[1].data['Z']

class OurGalaxy:
    def __init__(self, name, sfr, sfr_u, logmass, mass_upunc, mass_lowunc, v50, v98):
        self.name = name
        self.sfr = sfr
        self.sfr_u = sfr_u
        self.logmass = logmass
        self.mass_upunc = mass_upunc
        self.mass_lowunc = mass_lowunc
        self.v50 = v50
        self.v98 = v98

 # galaxies  J1125, J1232, J1450, J1713, and J2118 have no mass from pospector
J0826 = OurGalaxy('J0826+4305',  10.6,   1.2,    10.71,  0.02,   0.01,   850.21,    1220.36)
J0901 = OurGalaxy('J0901+0314',  6.5,    0.5,    10.59,  0.01,   0.01,   0,         0)
J0905 = OurGalaxy('J0905+5759',  14.3,   1.5,    10.75,  0.02,   0.02,   2476.77,   2913.18)
J0944 = OurGalaxy('J0944+0930',  5.2,    0.7,    10.74,  0.04,   0.03,   1241,      1747.79)
J1107 = OurGalaxy('J1107+0417',  5.8,    0.4,    10.68,  0.02,   0.02,   1626.14,   1995.31)
J1125 = OurGalaxy('J1125-0145',  2.5,    0.8,    0,      0,      0,      0,         0)
J1219 = OurGalaxy('J1219+0336',  3.8,    0.5,    11.34,  0.02,   0.02,   1586.07,   1818.08)
J1232 = OurGalaxy('J1232+0723',  2.48,   0.32,   0,      0,      0,      44.93,     347.02)
J1341 = OurGalaxy('J1341-0321',  19.7,   1.6,    10.61,  0.02,   0.02,   711.96,    1700.5)
J1450 = OurGalaxy('J1450+4621',  5,      4,      0,      0,      0,      406.4,     1189.78)
J1506 = OurGalaxy('J1506+5402',  21.3,   1.1,    10.65,  0.02,   0.03,   1998.43,   1499.38)
J1558 = OurGalaxy('J1558+3957',  6.61,   0.34,   10.97,  0.02,   0.02,   571.05,    1010.25)
J1613 = OurGalaxy('J1613+2834',  5.9,    0.4,    11.34,  0.02,   0.01,   1876.04,   2518.87)
J1713 = OurGalaxy('J1713+2817',  0.9,    0.7,    0,      0,      0,      414.55,    302.87)
J2116 = OurGalaxy('J2116-0634',  9.4,    2.7,    11.6,   0.04,   0.05,   0,         0)
J2118 = OurGalaxy('J2118+0017',  4.02,   0.34,   0,      0,      0,      0,         0)
J2140 = OurGalaxy('J2140+1209',  4.2,    2,      11.12,  0.05,   0.05,   282.11,    584.13)

# Organizing our data
ourgals = [J0826, J0901, J0905, J0944, J1107, J1125, J1219, J1232, J1341,
           J1450, J1506, J1558, J1613, J1713, J2116, J2118, J2140]
ourgalnames = [x.name for x in ourgals]
ourmass = [x.logmass for x in ourgals]
ourmasserr = [[x.mass_lowunc for x in ourgals], [x.mass_upunc for x in ourgals]]
ourV50 = np.array([x.v50 for x in ourgals])
ourV98 = np.array([x.v98 for x in ourgals])
bradna_Cmask = [True if x in ourgalnames else False for x in ShortNames]

# Data from Tremonti Table
TremontiGalaxies = np.array([[ShortNames[i], ChristyMass[i], ChristySFR[i], ChristyV50[i], ChristyV98[i]]
                             for i in range(0, len(ShortNames)) if ShortNames[i] in ourgalnames])
mass_T = [float(x[1]) for x in TremontiGalaxies]
HiresV98 = [ChristyV98[i] for i in range(0, len(ChristyV98)) if ShortNames[i] in ourgalnames]
HiresAge = [ChristyAge[i] for i in range(0, len(ChristyAge)) if ShortNames[i] in ourgalnames]

# The galaxies for which mass is 0 above have no mass calculated by prospector
# The following statement fills in those blanks with mass from the Tremonti table
for i in range(0, len(ourmass)):
    if ourmass[i] == 0:
        ourmass[i] = mass_T[i]



# Getting the flux #######

# function to return a flux in maggies given an AB magnitude
def flux(mag):
    flux = 3631 * 10. ** (mag / (-2.5))
    return flux

# Function to get luminosity (taken from cb galaxy fits sfr analysis tools)
# But we need vSv. WHat's v?

def get_fv(LV,Z):

    lum_dis = (cosmo.luminosity_distance(umeh_Z_T))*(1/u.Mpc)*(prsc*10**6)*(10**6)
    luminosity = (LV)*4*np.pi*(lum_dis**2)

    return luminosity

umeh_w4flux = flux(umeh_w4mag)

# Get frequencies for each galaxy
freq = ((lightspeed*(10**6))/22.194)*(1+umeh_Z_T)

# Get vobs fobs
vf = get_fv(umeh_w4flux*freq, umeh_Z_T)

#UmehMask
ChristyUmehMask = [True if i in umeh_shortnames else False for i in ShortNames]
ChrstUmehSFR = ChristySFR[ChristyUmehMask]
##### Got Flux and ARrays to Plot ##########


# Plotting...
plt.close('all')

# fig = plt.figure(figsize=(12, 9))
fig = plt.figure()
plt.tight_layout()

# --------------- PLOT 1
# Plotting Mass Tremonti vs Mass Pospector

ax = fig.add_subplot(2, 2, 1)

ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_ylim(10.4, 11.4)
ax.set_xlim(10.4, 11.4)

ax.set_xlabel("Mass Prosp")
ax.set_ylabel("Mass Tremonti")

# Use these instead once we get all mass values
# plt.ylim(np.amin([x.logmass for x in ourgals]+[x.logmass_T for x in ourgals])-0.3,
#          np.amin([x.logmass for x in ourgals]+[x.logmass_T for x in ourgals])+0.3)
# plt.xlim(np.amin([x.logmass for x in ourgals]+[x.logmass_T for x in ourgals])-0.3,
#          np.amin([x.logmass for x in ourgals]+[x.logmass_T for x in ourgals])+0.3)

# plot y=x line
ax.plot(np.linspace(9, 12, 16), np.linspace(9, 12, 16), "--", lw=0.5, color="black", alpha=0.3)
ax.scatter(ourmass, mass_T)
ax.errorbar(ourmass, mass_T, xerr=ourmasserr, ls='none')

# --------------- PLOT 2
# Plotting Stellar Age vs VMAX from Tremonti Data

for i in range(0, len(ourgals)):
    ax.annotate(ourgals[i].name, (ourmass[i]-0.02, mass_T[i]+0.02))

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_ylabel("VMax Tremonti")
ax2.set_xlabel("Stellar Age Tremonti")

ax2.scatter(ChristyAge, ChristyV98, c='blue')
ax2.scatter(HiresAge, HiresV98, c='red')
for i in range(0, len(ourgals)):
    if ShortNames[i] in ourgalnames:
        ax2.annotate(ourgals[i].name, (HiresAge[i]+0.02, HiresV98[i]-0.02))

# ------------------ PLOT 3
# plotting LV vs SFR
ax3 = fig.add_subplot(2, 2, 3)
ax3.scatter(ChrstUmehSFR, vf)
ax3.set_xlabel("Christy SFR")
ax3.set_ylabel(r"$\nu$ $L_{\nu}$")

# -------------------- PLOT 4
# Plotting my velocities versus Christy's

#Old plot plotting VBradna vs VTremonti
ax4 = fig.add_subplot(2, 2, 4)
#ax4.plot(np.linspace(-3100, 0, 100), np.linspace(-3100, 0, 100), "--", lw=0.5, color="black", alpha=0.3)
v50_zero_mask = (ourV50 != 0)
v98_zero_mask = (ourV98 != 0)
# ax4.scatter([(-1)*x for x in ourV50 if x != 0],
#             [ChristyV50[bradna_Cmask][i] for i in range(0, len(ourV50)) if ourV50[i] != 0], c='skyblue')
# ax4.scatter([(-1)*x for x in ourV98 if x != 0],
#             [ChristyV98[bradna_Cmask][i] for i in range(0, len(ourV98)) if ourV98[i] != 0], c='blue')
ax4.scatter(((-1)*ourV50[v50_zero_mask])-ChristyV50[bradna_Cmask][v50_zero_mask],
            ChristyV50[bradna_Cmask][v50_zero_mask], c='skyblue')
ax4.scatter(((-1)*ourV98[v98_zero_mask])-ChristyV98[bradna_Cmask][v98_zero_mask],
            ChristyV98[bradna_Cmask][v98_zero_mask], c='blue')



ax4.set_ylabel("Christy's Velocity")
ax4.set_xlabel("Bradna Velocity")



#plt.savefig("comparison.png", bbox_inches="tight", dpi=1800)
plt.show(block=False)



