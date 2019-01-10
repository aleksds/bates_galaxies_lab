from astropy.io import ascii
import numpy as np

ad_values = ascii.read('ad_mag_size_table.dat')

re_avg = np.array(ad_values['re'] * 0.025)
re_475 = np.array(ad_values['re_small'] * 0.025)
re_814 = np.array(ad_values['re_large'] * 0.025)

j2140 = np.array([re_475[11], re_avg[11], re_814[11]])
re_avg[11] = j2140[0]
re_475[11] = j2140[1]
re_814[11] = j2140[2]

print(re_avg)

# How much smaller are the 475 sizes?
dif_475 = re_475 - re_avg
print(dif_475)

per_475 = (dif_475)/re_avg

print(per_475)
print(np.median(per_475))

# How much larger are the 814 sizes?
dif_814 = re_814 - re_avg
print(dif_814)

per_814 = (re_814-re_avg)/re_avg

print(per_814)
print(np.median(per_814))

# How much smaller are the 475 sizes relative to the 814 sizes?
per_814_475 = (re_475 - re_814) / re_814

print(per_814_475)
print(np.median(per_814_475))
print(np.mean(per_814_475))
