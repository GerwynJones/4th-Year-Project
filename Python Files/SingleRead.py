#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import sys
import os
import numpy as np

import cPickle as pickle

import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 15})


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


### BEGIN ###

InputDir = "../Analysis/"  # raw_input("Enter either ../ or ./ : ")

T = 0

Size = 1

result = []

while T < 2:

    Input = raw_input("Enter name of pkl file location: ")

    File = raw_input("Enter name of pkl file: ")

    SaveLocation = InputDir + Input + "Plot_Files/"

    Filelocation = SaveLocation + File + ".pkl"

    with open(Filelocation, 'rb') as infile:
        result.append(pickle.load(infile))

    decision = raw_input("New File ? (y or n): ")

    if decision == 'y' or decision == 'yes':

        Size += 1

        T = 1

    else:

        if T == 1:

            File_Name = raw_input("Input Filename : ")

            if File_Name[0] == '/':
                File_Name = File_Name[1:]

            if File_Name[-1] == '/':
                File_Name = File_Name[:-1]

            SaveLocation = InputDir + "Multiple_File_Plot/" + File_Name + "/"

            ensure_dir(SaveLocation)

            T = 4

        else:

            T = 4



Max_Time = 45

print ("")

print("Here is what is Plotted \n")

print("There is %s Files Used" % Size)

print("MassPDF, CDF")


fig1 = plt.figure(figsize=(10, 8))

fig2 = plt.figure(figsize=(10, 8))

for i in xrange(Size):

    if len(result[i]) == 4:

        MassPDF, CDF, TIME, nbincent = result[i]

    else:

        print("ERROR : Incorrect File Used")

        sys.exit()

    ### Maximum Time = 45 code unit ###

    plt.loglog(nbincent, MassPDF)

    plt.xlabel(r'$n\ [cm^{-3}]$')

    plt.ylabel(r'$Mass-weighted\ PDF$')

    PDF_Save = SaveLocation + "Mass-weighted PDF.png"


    plt.loglog(nbincent, CDF)

    plt.xlabel(r'$n\ [cm^{-3}]$')

    plt.ylabel(r'$Cumulative\ PDF$')

    CDF_Save = SaveLocation + "Cumulative.png"


fig1.savefig(PDF_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

fig2.savefig(CDF_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
