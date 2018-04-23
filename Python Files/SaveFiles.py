#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""
from __future__ import division

import sys
import re
import os

import cPickle as pickle

from Arepo_Read_Master import *


### BEGIN ###

v_log = 1


def keyFunc(afilename):
    nondigits = re.compile("\D")
    return int(nondigits.sub("", afilename))


print('\nAvailable files:\n')

inputdir = raw_input("Enter directory: ")

print(" ")

files = []
files_sorted = []

inputfile = raw_input("Enter relevant letters of files: ")

print(" ")

for i in os.listdir(inputdir):

    if os.path.isfile(os.path.join(inputdir, i)) and inputfile in i:
        files.append(i)

for x, name in enumerate(sorted(files, key=keyFunc)):

    files_sorted.append(name)

Files_sorted = files_sorted[1 : ]

print(Files_sorted)

print(" ")

decision = raw_input("Continue ? : ")

if decision == "yes" or decision == "y" or decision == "":

    print("Continuing")

else:

    sys.exit()


init_file_dir = inputdir + Files_sorted[0]

_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, num_species = aread(init_file_dir)


TIME = np.zeros(len(Files_sorted))

TOTMASS_n = np.zeros((4, len(Files_sorted)))

TOTMASS_T = np.zeros((4, len(Files_sorted)))

CUM_MASS = np.zeros(len(Files_sorted))

Peak_Abundances = np.zeros((num_species, len(Files_sorted)))

TotMassCO = np.zeros(len(Files_sorted))

TotMassCplus = np.zeros(len(Files_sorted))

TotMassH2 = np.zeros(len(Files_sorted))

MassFraction = np.zeros((num_species, len(Files_sorted)))

Density_max = np.zeros(len(Files_sorted))

Temperature_min = np.zeros(len(Files_sorted))


for j in range(len(Files_sorted)):

    file_dir = inputdir + Files_sorted[j]

    DATA = aread(file_dir)

    TPOS = DATA[0]*(dist_cu/(ap.pc.to('cm')))  # Converting to parsecs

    TMASS = DATA[3]

    VEL = DATA[1]

    ID = DATA[2]

    U = DATA[4]

    RHO = DATA[5]

    POT = DATA[6]

    DIVV = DATA[7]

    ACCE = DATA[8]

    DUSTTEMP = DATA[9]

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


    POS = TPOS[:ngas, :]

    MASS = TMASS[:ngas]


    #### TEMP AND DENSITY plots and calculations ####

    ABHE = 0.1

    n = RHO * rho_cu / ((1.0 + 4.0 * ABHE) * m_p)

    Energy = (U * RHO * ergs_cu) / (dist_cu**3)

    ntot = (1.0 + ABHE - CHEM[0, :] + CHEM[1, :]) * n

    TEMP = 2.0 * Energy / (3.0 * ntot * k_B)


    a = 21

    n_sort = np.sort(ntot.value)

#    print n_sort[-a:]

    T_sort = np.sort(TEMP.value)

#    print T_sort[:a]

    Density_max[j] = np.mean(n_sort[-a:-11])

    Temperature_min[j] = np.mean(T_sort[11:a])


    #### MASS and Abundances plots and calculations ####

    TIME[j] = time * (time_cu.value/au.megayear.to('s'))


    #### Chemistry for NL97 and NL99 ####

    Max_Carbon = 1.4e-4

    Max_Hydrogen = 1

    Totmass = 2e4

    if num_species == 3:

        Peak_H2 = np.max(CHEM[0, :])

        Peak_Hplus = np.max(CHEM[1, :])

        Peak_CO = np.max(CHEM[2, :])

        Peak_Cplus = Max_Carbon - np.max(CHEM[2, :])


        Peak_Abundances[0, j] = Peak_H2

        Peak_Abundances[1, j] = Peak_Hplus

        Peak_Abundances[2, j] = Peak_CO

        Peak_Abundances[3, j] = Peak_Cplus


        TotMassCO[j] = np.sum(CHEM[2, :] * MASS)

        TotMassCplus[j] = np.sum((Max_Carbon - CHEM[2, :]) * MASS)

        TotMassH2[j] = np.sum(CHEM[0, :] * MASS)

        MassFraction[0, j] = TotMassCO[j]/(Max_Carbon*Totmass)

        MassFraction[1, j] = TotMassCplus[j]/(Max_Carbon*Totmass)

        MassFraction[2, j] = TotMassH2[j] / (Max_Hydrogen * Totmass)

    if num_species == 9:

        Peak_H2 = np.max(CHEM[0, :])

        Peak_Hplus = np.max(CHEM[1, :])

        Peak_Cplus = np.max(CHEM[2, :])

        Peak_CHx = np.max(CHEM[3, :])

        Peak_OHx = np.max(CHEM[4, :])

        Peak_CO = np.max(CHEM[5, :])

        Peak_HCOplus = np.max(CHEM[6, :])

        Peak_HEplus = np.max(CHEM[7, :])

        Peak_Mplus = np.max(CHEM[8, :])


        Peak_Abundances[0, j] = Peak_H2

        Peak_Abundances[1, j] = Peak_Hplus

        Peak_Abundances[2, j] = Peak_Cplus

        Peak_Abundances[3, j] = Peak_CHx

        Peak_Abundances[4, j] = Peak_OHx

        Peak_Abundances[5, j] = Peak_CO

        Peak_Abundances[6, j] = Peak_HCOplus

        Peak_Abundances[7, j] = Peak_HEplus

        Peak_Abundances[8, j] = Peak_Mplus

        TotMassCO[j] = np.sum(CHEM[5, :] * MASS)

        TotMassCplus[j] = np.sum((CHEM[2, :]) * MASS)

        TotMassH2[j] = np.sum(CHEM[0, :] * MASS)

        MassFraction[0, j] = TotMassCO[j]/(Max_Carbon*Totmass)

        MassFraction[1, j] = TotMassCplus[j]/(Max_Carbon*Totmass)

        MassFraction[2, j] = TotMassH2[j] / (Max_Hydrogen * Totmass)


    n_greaterthan10 = MASS[n.value >= 10]

    n_greaterthan1e2 = MASS[n.value >= 1e2]

    n_greaterthan1e3 = MASS[n.value >= 1e3]

    n_greaterthan1e4 = MASS[n.value >= 1e4]

    n_greaterthan1e5 = MASS[n.value >= 1e5]


    A = 20
    B = 100
    C = 500
    D = 2000


    T_lessthanA = MASS[TEMP.value <= A]

    T_betweenAandB = MASS[(TEMP.value > A) & (TEMP.value < B)]

    T_betweenBandC = MASS[(TEMP.value > C) & (TEMP.value < D)]

    T_overC = MASS[TEMP.value >= D]


    TOTMASS_n[0, j] = np.sum(n_greaterthan1e2)

    TOTMASS_n[1, j] = np.sum(n_greaterthan1e3)

    TOTMASS_n[2, j] = np.sum(n_greaterthan1e4)

    TOTMASS_n[3, j] = np.sum(n_greaterthan1e5)


    TOTMASS_T[0, j] = np.sum(T_lessthanA)

    TOTMASS_T[1, j] = np.sum(T_betweenAandB)

    TOTMASS_T[2, j] = np.sum(T_betweenBandC)

    TOTMASS_T[3, j] = np.sum(T_overC)


# InputDir = raw_input("Enter either ../ or ./ : ")

Input = raw_input("Enter name of pkl file location : ")

File = Input.replace('/', '_')

if File[0] == '_':

    File = File[1:]

if File[-1] == '_':

    File = File[:-1]

Filename = "../" + File + ".pkl"


data = [num_species,TIME, Density_max, Temperature_min, Peak_Abundances, MassFraction, TOTMASS_n, TOTMASS_T, A, B, C, D]

with open(Filename, 'wb') as outfile:
    pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)



