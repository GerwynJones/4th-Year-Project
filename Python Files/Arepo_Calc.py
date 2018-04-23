#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import matplotlib.pyplot as plt

from Arepo_Read_Master import *

### THE MEAT OF THE MATTER ###

filecount = 1

scalecount = 1

thin_factor = 100


while filecount > 0:

    inputdir = raw_input("Enter directory: ")

    print '\nAvailable files:\n'

    print os.listdir(str(inputdir))

    os_dir = os.listdir(str(inputdir))

    inputfile = raw_input("Enter filename: ")

    file_dir = inputdir+inputfile

    DATA = aread(file_dir)

    POS = DATA[0]

    VEL = DATA[1]

    ID = DATA[2]

    MASS = DATA[3]

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

    POS = POS[:ngas, :]

    MASS = MASS[:ngas]

    ### TEMP ###

    ABHE = 0.1

    n = RHO * rho_cu / ((1.0 + 4.0 * ABHE) * m_p)

    Energy = (U * RHO * ergs_cu) / (dist_cu**3)

    ntot = (1.0 + ABHE - CHEM[0, :] + CHEM[1, :]) * n

    TEMP = 2.0 * Energy / (3.0 * ntot * k_B)

    # print np.min(RHO)
    # print np.max(RHO)
    #
    # print np.min(n)
    # print np.max(n)
    #
    # print np.min(ntot)
    # print np.max(ntot)
    #
    # print np.amin(U)
    # print np.amax(U)
    #
    # print np.min(energy)
    # print np.max(energy)
    #
    # print np.min(TEMP)
    # print np.max(TEMP)

    parameter = raw_input("Enter Calculation that you want to find (Time, CHEM, RHO, TEMP) : ")

    if parameter == 'time' or parameter == 'TIME':

        print time

    if parameter == 'chem' or parameter == 'CHEM':

        X = np.int(raw_input("Enter density (cm^{-3}) : "))

        Y = np.int(raw_input("Enter density range (cm^{-3}) : "))

        print "Options for chem species", np.arange(nspecies)

        Z = np.int(raw_input("Enter chem species (See options above) : "))

        CHEM_Range = CHEM[Z, :][(ntot >= (X - Y)) & (ntot <= (X + Y))]

        print np.mean(CHEM_Range)

        print "min chem of species", np.min(CHEM[Z, :])
        print "max chem of species", np.max(CHEM[Z, :])

    if parameter == 'rho' or parameter == 'RHO':

        RHO_TRUE = (RHO * rho_cu.value)

    if parameter == 'temp' or parameter == 'TEMP':

        X = np.int(raw_input("Enter density (cm^{-3}) : "))

        Y = np.int(raw_input("Enter density range (cm^{-3}) : "))

        TEMP_Range = TEMP[(ntot >= (X - Y)) & (ntot <= (X + Y))]

        print np.mean(TEMP_Range)

    if parameter == 'no-of-sinks' or parameter == 'NO-of-Sinks' or parameter == 'sinks':

        print nsink


    file_decision = raw_input("Make a new decision?: ")

    if file_decision == 'y':

        continue

    else:

        break
