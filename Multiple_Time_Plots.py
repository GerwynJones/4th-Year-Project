#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import sys

import matplotlib.pyplot as plt

from Arepo_Read_Master import *


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

if decision == "yes" or decision == "y" or decision == "":

    print "Continue"

else:

    sys.exit()


TIME = np.zeros(len(files))

TOTMASS_n = np.zeros((2, len(files)))

TOTMASS_T = np.zeros((2, len(files)))

CUM_MASS = np.zeros(len(files))

Peak_Abundances = np.zeros((1, len(files)))

Density = np.zeros(len(files))

Temperature = np.zeros(len(files))


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

    ### TEMP AND DENSITY plots and calculations ###

    ABHE = 0.1

    n = RHO * rho_cu / ((1.0 + 4.0 * ABHE) * m_p)

    Energy = (U * RHO * ergs_cu) / (dist_cu**3)

    ntot = (1.0 + ABHE - CHEM[0, :] + CHEM[1, :]) * n

    TEMP = 2.0 * Energy / (3.0 * ntot * k_B)


    Density[j] = np.max(ntot.value)

    Temperature[j] = np.min(TEMP.value)


    #### MASS and Abundances plots and calculations ####

    TIME[j] = time

    Peak_H2 = np.max(CHEM[0, :])

    Peak_Abundances[0, j] = Peak_H2

    n_lessthan1e2 = n.value[n.value <= 1e2]

    n_lessthan1e3 = n.value[n.value <= 1e3]

    T_lessthan30 = TEMP.value[TEMP.value <= 30]

    T_between30and350 = TEMP.value[(TEMP.value > 30) & (TEMP.value < 350)]

    TOTMASS_n[0, j] = np.sum(n_lessthan1e2)

    TOTMASS_n[1, j] = np.sum(n_lessthan1e3)

    TOTMASS_T[0, j] = np.sum(T_lessthan30)

    TOTMASS_T[0, j] = np.sum(T_between30and350)



plt.figure()

plt.semilogy(TIME, Density, marker='x', linestyle='None', label="Density")

plt.semilogy(TIME, Temperature, marker='x', linestyle='None', label="Temperature")

plt.xlabel(r'$Time \/ (CODE)$')

plt.ylabel(r'$ Density/Temperature $')

plt.legend(loc='best')


plt.figure()

plt.semilogy(TIME, Peak_Abundances[0, :], marker='x', linestyle='None', label="H2")

plt.xlabel(r'$Time \/ (CODE)$')

plt.ylabel(r'$Peak Abundances \/$')

plt.legend(loc='best')

plt.figure()

plt.plot(TIME, TOTMASS_n[0, :], marker='x', linestyle='None')

plt.plot(TIME, TOTMASS_n[1, :], marker='x', linestyle='None')

plt.xlabel(r'$Time \/ (CODE)$')

plt.ylabel(r'$Mass \/ $')

plt.figure()

plt.plot(TIME, TOTMASS_T[0, :], marker='x', linestyle='None')

plt.plot(TIME, TOTMASS_T[1, :], marker='x', linestyle='None')

plt.xlabel(r'$Time \/ (CODE)$')

plt.ylabel(r'$Mass \/ $')



plt.show()

