#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import matplotlib
matplotlib.use('TkAgg')

import sys
import re
import os
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 15})

from Arepo_Read_Master import *

### BEGIN ###

filecount = 1

v_log = 1

def keyFunc(afilename):
    nondigits = re.compile("\D")
    return int(nondigits.sub("", afilename))

files = []
files_sorted = []

print(" ")

inputdir = raw_input("Enter directory: ")

file = raw_input("Enter relevant letters of files: ")

print(" ")

for i in os.listdir(inputdir):

    if os.path.isfile(os.path.join(inputdir, i)) and file in i:
        files.append(i)

for x, name in enumerate(sorted(files, key=keyFunc)):

    files_sorted.append(name)

Files_sorted = files_sorted[1 : ]

print('\nAvailable files:\n')

print(Files_sorted)

print(" ")

inputfile = raw_input("Enter filename: ")

file_dir = inputdir+inputfile

DATA = aread(file_dir)

TPOS = DATA[0]

VEL = DATA[1]

ID = DATA[2]

TMASS = DATA[3]

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

TIME = time * (time_cu.value / au.megayear.to('s'))


### TEMP ###

ABHE = 0.1

n = RHO * rho_cu / ((1.0 + 4.0 * ABHE) * m_p)

Energy = (U * RHO * ergs_cu) / (dist_cu**3)

ntot = (1.0 + ABHE - CHEM[0, :] + CHEM[1, :]) * n

TEMP = 2.0 * Energy / (3.0 * ntot * k_B)


### ###

x_min = np.min(POS[:, 0])
x_max = np.max(POS[:, 0])

y_min = np.min(POS[:, 1])
y_max = np.max(POS[:, 1])

z_min = np.min(POS[:, 2])
z_max = np.max(POS[:, 2])

while filecount > 0:

    ### ###

    parameter = raw_input("Enter data tag to plot (Default = Density; options = Dusttemp, Chemplot, Chem, Mass, Temp, PDF or 0 to change file or exit) : ")

    print(" ")

    if parameter == '0':

        break

    elif parameter == '' or parameter == 'density' or parameter == 'Density':

        Density = n.value

        v_min, v_max = np.min(Density), np.max(Density)

        PD = np.hstack((POS, Density.reshape(Density.size, 1)))

        PDX = PD[np.argsort(PD[:, 3])]

        ax1_decision = raw_input("Choose Axis (X=0, Y=1, Z=2) : ")
        ax2_decision = raw_input("Choose Axis (X=0, Y=1, Z=2) : ")

        ax1 = np.int(ax1_decision)
        ax2 = np.int(ax2_decision)

        fig, ax = plt.subplots(figsize=(10, 8))

        if v_log == 1:

            plt.scatter(PDX[:, ax1], PDX[:, ax2], c=np.log10(PDX[:, 3]),
                        vmin=np.log10(v_min), vmax=np.log10(v_max), marker='.', s=1, cmap='jet')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$log_{10}(density)$")

            if nsink > 0:

                plt.scatter(SINKPOS[:, ax1], SINKPOS[:, ax2], c='k', marker='*', s=3)

        else:

            plt.scatter(PDX[:, ax1], PDX[:, ax2], c=PDX[:, 3],
                        vmin=v_min, vmax=v_max, marker='.', s=1, cmap='jet')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$density$")

            if nsink > 0:

                plt.scatter(SINKPOS[:, ax1], SINKPOS[:, ax2], c='k', marker='*', s=3)

        plt.xlim((x_min, x_max))

        plt.ylim((y_min, y_max))

        if ax1 == 0:

            lab1 = "x"

        elif ax1 == 1:

            lab1 = "y"

        elif ax1 == 2:

            lab1 = "z"

        if ax2 == 0:

            lab2 = "x"

        elif ax2 == 1:

            lab2 = "y"

        elif ax2 == 2:

            lab2 = "z"

        plt.xlabel(r"$%s\ [pc]$" % lab1)

        plt.ylabel(r"$%s\ [pc]$" % lab2)

        plt.show()

        xy_decision = raw_input("Change scale? : ")

        if xy_decision == 'y':

            x_min = input("xmin: ")

            x_max = input("xmax: ")

            y_min = input("ymin: ")

            y_max = input("ymax: ")

            v_log = input("Colour scale (0 = linear, 1 = log10) : ")

            continue

        else:

            break

    elif parameter == 'chemplot' or parameter == 'CHEMPLOT':

        print("Options for chem species", np.arange(nspecies))

        ### 0 for H2    1 for H+    2 for CO ###

        X = raw_input("Choose which chem species to plot (See above for options) : ")

        v_min, v_max = np.min(CHEM[np.int(X), :]), np.max(CHEM[np.int(X), :])

        POSCHEM = np.hstack((POS, CHEM[np.int(X), :].reshape(CHEM[np.int(X), :].size, 1)))

        POSCHEMX = POSCHEM[np.argsort(POSCHEM[:, 3])]

        ax1_decision = raw_input("Choose Axis (X=0, Y=1, Z=2) : ")
        ax2_decision = raw_input("Choose Axis (X=0, Y=1, Z=2) : ")

        ax1 = np.int(ax1_decision)
        ax2 = np.int(ax2_decision)

        plt.figure(figsize=(10, 8))

        if v_log == 1:

            plt.scatter(POSCHEMX[:, ax1], POSCHEMX[:, ax2], c=np.log10(POSCHEMX[:, 3] + 1e-7), marker='.', s=1, cmap='jet')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$log_{10}(CO Abundance)$")

            if nsink > 0:

                plt.scatter(SINKPOS[:, ax1], SINKPOS[:, ax2], c='k', marker='*', s=3)

        else:

            plt.scatter(POSCHEMX[:, ax1], POSCHEMX[:, ax2], c=(POSCHEMX[:, 3] + 1e-7), marker='.', s=1, cmap='jet')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$CO Abundance$")

            if nsink > 0:

                plt.scatter(SINKPOS[:, ax1], SINKPOS[:, ax2], c='k', marker='*', s=3)

        plt.xlim((x_min, x_max))

        plt.ylim((y_min, y_max))

        if ax1 == 0:

            lab1 = "x"

        elif ax1 == 1:

            lab1 = "y"

        elif ax1 == 2:

            lab1 = "z"

        if ax2 == 0:

            lab2 = "x"

        elif ax2 == 1:

            lab2 = "y"

        elif ax2 == 2:

            lab2 = "z"

        plt.xlabel(r"$%s\ [pc]$" % lab1)

        plt.ylabel(r"$%s\ [pc]$" % lab2)

        plt.show()

        xy_decision = raw_input("Change scale? : ")

        if xy_decision == 'y':

            x_min = input("xmin: ")

            x_max = input("xmax: ")

            y_min = input("ymin: ")

            y_max = input("ymax: ")

            v_log = input("Colour scale (0 = linear, 1 = log10) : ")

            continue

        else:

            break

    elif parameter == 'dusttemp' or parameter == 'DUSTTEMP':

        v_min, v_max = np.min(DUSTTEMP), np.max(DUSTTEMP)

        POSDUSTTEMP = np.hstack((POS, DUSTTEMP.reshape(DUSTTEMP.size, 1)))

        POSDUSTTEMPX = POSDUSTTEMP[np.argsort(POSDUSTTEMP[:, 3])]

        ax1_decision = raw_input("Choose Axis (X=0, Y=1, Z=2) : ")
        ax2_decision = raw_input("Choose Axis (X=0, Y=1, Z=2) : ")

        ax1 = np.int(ax1_decision)
        ax2 = np.int(ax2_decision)

        plt.figure(figsize=(10, 8))

        if v_log == 1:

            plt.scatter(POSDUSTTEMPX[:, ax1], POSDUSTTEMPX[:, ax2], c=np.log10(POSDUSTTEMPX[:, 3]),
                        vmin=np.log10(v_min), vmax=np.log10(v_max), marker='.', s=1, cmap='jet')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$log_{10}(Temp\ K)$")

            if nsink > 0:

                plt.scatter(SINKPOS[:, ax1], SINKPOS[:, ax2], c='k', marker='*', s=3)

        else:

            plt.scatter(POSDUSTTEMPX[:, ax1], POSDUSTTEMPX[:, ax2], c=POSDUSTTEMPX[:, 3],
                        vmin=v_min, vmax=v_max, marker='.', s=1, cmap='jet')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$Temp\ K$")

            if nsink > 0:

                plt.scatter(SINKPOS[:, ax1], SINKPOS[:, ax2], c='k', marker='*', s=3)

        plt.xlim((x_min, x_max))

        plt.ylim((y_min, y_max))

        if ax1 == 0:

            lab1 = "x"

        elif ax1 == 1:

            lab1 = "y"

        elif ax1 == 2:

            lab1 = "z"

        if ax2 == 0:

            lab2 = "x"

        elif ax2 == 1:

            lab2 = "y"

        elif ax2 == 2:

            lab2 = "z"

        plt.xlabel(r"$%s\ [pc]$" % lab1)

        plt.ylabel(r"$%s\ [pc]$" % lab2)

        plt.show()

        xy_decision = raw_input("Change scale? : ")

        if xy_decision == 'y':

            x_min = input("xmin: ")

            x_max = input("xmax: ")

            y_min = input("ymin: ")

            y_max = input("ymax: ")

            v_log = input("Colour scale (0 = linear, 1 = log10) : ")

            continue

        else:

            break


    elif parameter == 'chem' or parameter == 'CHEM':

        print("Options for chem species", np.arange(nspecies))

        ### 0 for H2    1 for H+    2 for CO ###

        X = raw_input("Choose which chem species to plot (See above for options) : ")

        plt.figure(figsize=(10, 8))

        if v_log == 1:

            plt.scatter(np.log10(n.value), np.log10(TEMP.value), c=np.log10(CHEM[np.int(X), :] + 1e-7), marker='.', s=0.5,
                        cmap='jet_r')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$log_{10}(CO\ Abundance)$")

        else:

            plt.scatter(n.value, TEMP.value, c=(CHEM[np.int(X), :] + 1e-7), marker='.', s=0.5,
                        cmap='jet_r')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$CO\ Abundance$")

        plt.xlabel(r'$n\ [cm^{-3}]$')

        plt.ylabel(r'$Temperature\ [K]$')

        plt.show()

        xy_decision = raw_input("Change plot? : ")

        if xy_decision == 'y':

            v_log = 0

            continue

        else:

            break

    elif parameter == 'mass' or parameter == 'MASS':

        plt.figure(figsize=(10, 8))

        if v_log == 1:

            plt.scatter(np.log10(n.value), np.log10(TEMP.value), c=MASS, marker='.', s=0.5, cmap='jet_r')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$log_{10}(Mass)$")

        else:

            plt.scatter(n.value, TEMP.value, c=MASS, marker='.', s=0.5, cmap='jet_r')

            cbar = plt.colorbar()

            cbar.ax.set_ylabel(r"$Mass$")

        plt.xlabel(r'$n\ [cm^{-3}]$')

        plt.ylabel(r'$Temperature\ [K]$')

        plt.show()

        xy_decision = raw_input("Change plot? : ")

        if xy_decision == 'y':

            v_log = 0

            continue

        else:

            break


    elif parameter == 'pdf' or parameter == 'PDF':

        plt.figure(figsize=(10, 8))

        n_pdf = n.value[n.value >= 1]

        nmin = np.min(n_pdf)

        nmax = np.max(n_pdf)

        nbins = 150

        dlogn = np.log10(nmax - nmin) / nbins

        nbinmin = nmin

        Massbin = np.zeros(nbins)

#        nbinmax = np.zeros(nbins)

        nbincent = np.zeros(nbins)

        for i in xrange(1, nbins + 1):

            nbinmax = 10 ** (i * dlogn)

            nbincent[i-1] = nbinmin + (nbinmax-nbinmin)/2

            ParticleID = np.argwhere((n_pdf >= nbinmin) & (n_pdf <= nbinmax))

            MassID = MASS[ParticleID]

            Massbin[i-1] = np.sum(MassID)

            nbinmin = nbinmax

#        print(nbincent[-1])

#        print(nbinmax)

        MassPDF = Massbin/np.sum(Massbin)

        plt.loglog(nbincent, MassPDF)

        plt.xlabel(r'$n\ [cm^{-3}]$')

        plt.ylabel(r'$Mass-weighted\ PDF$')

        plt.show()

        xy_decision = raw_input("Change plot? : ")

        if xy_decision == 'y':

            v_log = 0

            continue

        else:

            break

        # elif parameter == 'POS':

    else:

        print('\n       NOT YET DEFINED')

        continue

