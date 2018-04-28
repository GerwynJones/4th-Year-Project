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

InputDir = "../Analysis/"    # raw_input("Enter either ../ or ./ : ")

T = 0

Size = 1

result = []

while T < 2:

    Input = raw_input("Enter name of pkl file location: ")

    File = Input.replace('/', '_')

    if Input[0] == '/':

        Input = Input[1:]

    if Input[-1] == 's':

        Input = Input + '/'

    if File[0] == '_':

        File = File[1:]

    if File[-1] == '_':

        File = File[:-1]

    SaveLocation = InputDir+Input+"Plot_Files/"

    Filelocation = SaveLocation+File+".pkl"

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


"""
Result = np.zeros(Size)

for i in xrange(len(Result)):

    print len(result[i])

    N = len(result[i])

    Result[i] = np.pad(Result[i], (0, N), 'constant')

    print Result[i]

    Result[i] = result[i]
"""

Max_Time = 45


print ("")

print("Here is what is Plotted \n")

print("There is %s Files Used" % Size)

print("Species, Time, Max_Density, Min_Temperature, Peak_Abundances, Mass_Fraction, Total_Mass_density, Total_Mass_temperature, No. Of. Sinks")



fig1, ax1 = plt.subplots()

fig2, ax2 = plt.subplots()

fig3, ax3 = plt.subplots()

ax4 = ax3.twinx()

fig4, ax5 = plt.subplots()

fig5, ax6 = plt.subplots()

fig6, ax7 = plt.subplots()


linestyle = ['solid', 'dashdot']

linecolor = ['red', 'green', 'blue', 'orange']

Blank = None

labeln1 = [r"$n\ >\ 1x10^2\ cm^{-3}$"]
labeln2 = [r"$n\ >\ 1x10^3\ cm^{-3}$"]
labeln3 = [r"$n\ >\ 1x10^4\ cm^{-3}$"]
labeln4 = [r"$n\ >\ 1x10^5\ cm^{-3}$"]

labelT1 = ['']
labelT2 = ['']
labelT3 = ['']
labelT4 = ['']

labelCO = ['CO']
labelCplus = ['C+']
labelH2 = ['H2']

Max_Massn = np.zeros(Size)
Max_MassT = np.zeros(Size)

MaxTemp = np.zeros(Size)
MinTemp = np.zeros(Size)


for i in xrange(Size):

    if len(result[i]) == 14:

        num_species, Time_Code, Time, Density_max, Temperature_min, Peak_Abundances, MassFraction, TOTMASS_n, TOTMASS_T, A, B, C, D, NSink = result[i]

    else:

        print("ERROR : Incorrect File Used")

        sys.exit()

    ### Maximum Time = 45 code unit ###

    ID = np.argwhere(Time_Code <= Max_Time)

    Time = Time[ID]; Density_max = Density_max[ID]; Temperature_min = Temperature_min[ID]; Peak_Abundances = Peak_Abundances[:, ID]

    MassFraction = MassFraction[:, ID]; TOTMASS_n = TOTMASS_n[:, ID]; TOTMASS_T = TOTMASS_T[:, ID]; NSink = NSink[ID]


    Max_Massn[i] = np.max(TOTMASS_n)
    Max_MassT[i] = np.max(TOTMASS_T)
    MaxTemp[i] = np.max(Temperature_min)
    MinTemp[i] = np.min(Temperature_min)

    labelT1[0] = r"$T\ <\ %s K$" % A
    labelT2[0] = r"$%s K\ < T\ <\ %s K$" % (A, B)
    labelT3[0] = r"$%s K\ <\ T\ <\ %s K$" % (C, D)
    labelT4[0] = r"$T\ >\ %s K$" % D


    if i > 0:

        labelT1.append(Blank)
        labelT2.append(Blank)
        labelT3.append(Blank)
        labelT4.append(Blank)

        labeln1.append(Blank)
        labeln2.append(Blank)
        labeln3.append(Blank)
        labeln4.append(Blank)

        labelT1.append(Blank)
        labelT2.append(Blank)
        labelT3.append(Blank)
        labelT4.append(Blank)

        labelCO.append(Blank)
        labelCplus.append(Blank)
        labelH2.append(Blank)

    else :
        print("")

    #### Sorting out the ordering ####

    if num_species == 3:

        ax1.semilogy(Time, Peak_Abundances[0, :], marker='None', linestyle=linestyle[i], color=linecolor[0], label=labelH2[i])

        ax1.semilogy(Time, Peak_Abundances[2, :], marker='None', linestyle=linestyle[i], color=linecolor[1], label=labelCO[i])

        """
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height])
    
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

        lgd = ax1.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

        ax1.set_xlabel(r'$Time\ [Myr]$')

        ax1.set_ylabel(r'$Peak\ Abundances$')

        Peak_Save = SaveLocation + "Peak Abundances.png"


        ax2.semilogy(Time, MassFraction[0, :], marker='None', linestyle=linestyle[i], color=linecolor[0], label=labelCO[i])

        ax2.semilogy(Time, MassFraction[1, :], marker='None', linestyle=linestyle[i], color=linecolor[1], label=labelCplus[i])

        ax2.semilogy(Time, 2*MassFraction[2, :], marker='None', linestyle=linestyle[i], color=linecolor[2], label=labelH2[i])


        """
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height])
    
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

        lgd = ax2.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

        ax2.set_xlabel(r'$Time\ [Myr]$')

        ax2.set_ylabel(r'$Mass\ Fraction$')

        Fraction_Save = SaveLocation + "Mass Fraction.png"


    if num_species == 9:

        ax1.semilogy(Time, Peak_Abundances[0, :], marker='None', linestyle=linestyle[i], color=linecolor[0], label=labelH2[i])

        ax1.semilogy(Time, Peak_Abundances[5, :], marker='None', linestyle=linestyle[i], color=linecolor[1], label=labelCO[i])

        """
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height])
    
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

        lgd = ax1.legend(loc='center right', bbox_to_anchor=(1.4, 0.5))

        ax1.set_xlabel(r'$Time\ [Myr]$')

        ax1.set_ylabel(r'$Peak\ Abundances$')

        Peak_Save = SaveLocation + "Peak Abundances.png"


        ax2.semilogy(Time, MassFraction[0, :], marker='None', linestyle=linestyle[i], color=linecolor[0], label=labelCO[i])

        ax2.semilogy(Time, MassFraction[1, :], marker='None', linestyle=linestyle[i], color=linecolor[1], label=labelCplus[i])

        ax2.semilogy(Time, 2*MassFraction[2, :], marker='None', linestyle=linestyle[i], color=linecolor[2], label=labelH2[i])

        """
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height])
    
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

        lgd = ax2.legend(loc='center right', bbox_to_anchor=(1.4, 0.5))

        ax2.set_xlabel(r'$Time\ [Myr]$')

        ax2.set_ylabel(r'$Mass\ Fraction$')

        Fraction_Save = SaveLocation + "Mass Fraction.png"


    ax3.semilogy(Time, Density_max, color='b', marker='None', linestyle=linestyle[i])

    ax3.set_xlabel(r'$Time\ [Myr]$')

    ax3.set_ylabel(r'$ Density\ \left[ cm^{-3} \right]$')

    #ax3.tick_params('y', colors='b')


    ax4.semilogy(Time, Temperature_min, color='r', marker='None', linestyle=linestyle[i])

    #ax4.tick_params('y', colors='r')

    fig3.tight_layout()

    DandT_Save = SaveLocation + "Density and Temperature Plot.png"


    ax5.semilogy(Time, TOTMASS_n[0, :], marker='None', linestyle=linestyle[i], color=linecolor[0], label=labeln1[i])

    ax5.semilogy(Time, TOTMASS_n[1, :], marker='None', linestyle=linestyle[i], color=linecolor[1], label=labeln2[i])

    ax5.semilogy(Time, TOTMASS_n[2, :], marker='None', linestyle=linestyle[i], color=linecolor[2], label=labeln3[i])

    ax5.semilogy(Time, TOTMASS_n[3, :], marker='None', linestyle=linestyle[i], color=linecolor[3], label=labeln4[i])

    """
    # Shrink current axis by 20%def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    lgd = ax5.legend(loc='center right', bbox_to_anchor=(1.58, 0.5))

    ax5.set_xlabel(r'$Time\ [Myr]$')

    ax5.set_ylabel(r'$Mass\ [M_{\odot}]$')

    Density_Save = SaveLocation + "Mass evolution of Gas with Density.png"

    ax6.semilogy(Time, TOTMASS_T[0, :], marker='None', linestyle=linestyle[i], color=linecolor[0], label=labelT1[i])

    ax6.semilogy(Time, TOTMASS_T[1, :], marker='None', linestyle=linestyle[i], color=linecolor[1], label=labelT2[i])

    ax6.semilogy(Time, TOTMASS_T[2, :], marker='None', linestyle=linestyle[i], color=linecolor[2], label=labelT3[i])

    ax6.semilogy(Time, TOTMASS_T[3, :], marker='None', linestyle=linestyle[i], color=linecolor[3], label=labelT4[i])

    """
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    lgd = ax6.legend(loc='center right', bbox_to_anchor=(1.62, 0.5))

    ax6.set_xlabel(r'$Time\ [Myr]$')

    ax6.set_ylabel(r'$Mass\ [M_{\odot}]$')

    Temperature_Save = SaveLocation + "Mass evolution of Gas with Temperature.png"


    ax7.semilogy(Time, NSink, marker='None', linestyle=linestyle[i], color=linecolor[0])

    """
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    ax7.set_xlabel(r'$Time\ [Myr]$')

    ax7.set_ylabel(r'$Number\ of\ Sinks$')

    NOS_Save = SaveLocation + "Number of Sinks.png"


MaxT = np.max(Max_MassT)
Maxn = np.max(Max_Massn)

a = 1.5

if Maxn > MaxT:

    MaxT = Maxn*a
    Maxn = Maxn*a

else:

    Maxn = MaxT*a
    MaxT = MaxT*a

Minn = 1.0e-2
MinT = 1.0e-2


T_Q = (Time[-1] - Time[0]) / 4

ax3.text(T_Q, np.min(Density_max), r'$max\ \rho$')

ax4.text(T_Q, np.max(Temperature_min), r'$min\ T$')


ax4.set_ylabel('Temperature [K]')

ax5.set_ylim(Minn, Maxn)

ax6.set_ylim(MinT, MaxT)


fig1.savefig(Peak_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

fig2.savefig(Fraction_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

fig3.savefig(DandT_Save, dpi=300, format='png', bbox_inches='tight')

fig4.savefig(Density_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

fig5.savefig(Temperature_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

fig6.savefig(NOS_Save, dpi=300, format='png', bbox_inches='tight')

