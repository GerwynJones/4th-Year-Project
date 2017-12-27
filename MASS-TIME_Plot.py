#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import sys

import matplotlib.pyplot as plt

from arepo_read_master import *


### THE MEAT OF THE MATTER ###

filecount = 1

scalecount = 1

thin_factor = 100

v_log = 1



print '\nAvailable files:\n'

path = inputfile = raw_input("Enter path: ")

files = []

inputdir = raw_input("Enter filename: ")

for i in os.listdir(path):

    if os.path.isfile(os.path.join(path, i)) and inputdir in i:
        files.append(i)

print files



decision = raw_input("Continue ? : ")

if decision == "yes" or decision == "y":

    print "Continue"

else:

    sys.exit()

TIME = np.zeros(len(files))

TOTMASS_n = np.zeros((2, len(files)))

TOTMASS_T = np.zeros((2, len(files)))


for j in range(len(files)):

    file_dir = path + files[j]

    print file_dir

    DATA = aread(file_dir)

    pos = DATA[0]

    mass = DATA[3]

    VEL = DATA[1]

    ID = DATA[2]

    U = DATA[4]

    RHO = DATA[5]

    POT = DATA[6]

    DIVV = DATA[7]

    ACCE = DATA[8]

    DUST = DATA[9]

    TSTP = DATA[10]

    BFLD = DATA[11]

    DIVB = DATA[12]

    SOFT = DATA[13]

    CHEM = DATA[14]

    time = DATA[15][0]

    nsink = DATA[16]

    ngas = DATA[17]

    SINKPOS = DATA[18]

    SINKVEL = DATA[19]

    SINKMASS = DATA[20]

    SINKID = DATA[21]

    nspecies = DATA[22]

    POS = pos[:ngas, :]

    MASS = mass[:ngas]

    ### TEMP ###

    ABHE = 0.1

    n = RHO * rho_cu / ((1.0 + 4.0 * ABHE) * m_p)

    Energy = (U * RHO * ergs_cu) / (dist_cu**3)

    ntot = (1.0 + ABHE - CHEM[0, :] + CHEM[1, :]) * n

    TEMP = 2.0 * Energy / (3.0 * ntot * k_B)


    #### MASS and Abundances plots and calculations ####










