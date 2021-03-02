# Author: Yuewen FANG
# Date: 2021 March 2nd
# The used muc in EPW code of estiamtion of Tc is up to 0.2,
# Here, I used larger muc to estimate Tc
#logavg and l_a2F are obtained from epw.out
import math
logavg =    0.0022479 # modify
l_a2F=0.9736146 # modify


####constants####
ryd2ev = 13.6056981  
kelvin2eV = 8.6173427909E-05
###############

for i in range(1,25):
    mu = 0.1 + 0.02 * (i - 1)
    tc = logavg/1.2 * math.exp( -1.04*(1.0 + l_a2F)/(l_a2F - mu*(1+0.62*l_a2F)))
    tc = tc * ryd2ev / kelvin2eV
    print("mu = %.2f, Tc = %.8f" % (mu, tc))
