# distributing 131 galaxies over 10,000 square degrees

ngal = 131
area_sdss = 1e4

density = ngal/area_sdss
print('galaxies per square degree: ', density)

# COSMOS
area_cosmos = 1.8

p_cosmos = area_cosmos / area_sdss
print('probability of one galaxy being in the COSMOS footprint', p_cosmos)

# probability of 0 galaxies in cosmos
p_0 = (1-p_cosmos)**ngal
print('probability that none of the galaxies are in COSMOS: ', p_0)
print('probability of one or more galaxies in COSMOS', 1-p_0)

# probability of 1 galaxy in cosmos

p_1 = ngal * p_cosmos * (1-p_cosmos)**(ngal-1)
print('probability of 1 galaxies in COSMOS: ', p_1)

# probability of 2 galaxies is cosmos

p_2 = 131 * 130. / (2) * p_cosmos**2 *  (1-p_cosmos)**(ngal-2)
print('probability of 2 galaxies in cosmos: ', p_2)

# CANDELS
area_candels = 800. / 60**2
p_candels = area_candels / area_sdss
p_0_candels = (1-p_candels)**ngal
print('probability that none of the galaxies are in CANDELS: ', p_0_candels)
print('probability that one or more galaxies are in CANDELS: ', 1-p_0_candels)

# probability of 1 galaxy in candels

p_1_candels = ngal * p_candels * (1-p_candels)**(ngal-1)
print('probability that 1 galaxy is in CANDELS: ', p_1_candels)
