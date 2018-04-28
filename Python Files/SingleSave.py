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

def Ensure_Dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

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

### Calculate ###

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

    nbincent[i - 1] = nbinmin + (nbinmax - nbinmin) / 2

    ParticleID = np.argwhere((n_pdf >= nbinmin) & (n_pdf <= nbinmax))

    MassID = MASS[ParticleID]

    Massbin[i - 1] = np.sum(MassID)

    nbinmin = nbinmax

#        print(nbincent[-1])

#        print(nbinmax)

MassPDF = Massbin / np.sum(Massbin)

CDF = np.cumsum(MassPDF)

print("Time = %.3f Myr \n" % TIME)


Ensure_Dir("../Plot_Files/")

File = "SingleSave_%.3f" % TIME

Filename = "../Plot_Files/" + File + ".pkl"


data = [MassPDF, CDF, TIME, nbincent]

with open(Filename, 'wb') as outfile:
    pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)



