#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import os


import cPickle as pickle

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 15})

### BEGIN ###

InputDir = "../Analysis/"    # raw_input("Enter either ../ or ./ : ")

Input = raw_input("Enter name of pkl file location: ")

File = Input.replace('/', '_')

if Input[0] == '/':

    Input = Input[1:]

if File[0] == '_':

    File = File[1:]

if File[-1] == '_':

    File = File[:-1]

SaveLocation = InputDir+Input

Filelocation = SaveLocation+File+".pkl"


with open(Filelocation, 'rb') as infile:
    result = pickle.load(infile)



print ("")

print("Here is what is Plotted \n")

print("Species, Time, Max_Density, Min_Temperature, Peak_Abundances, Mass_Fraction, Total_Mass_density, Total_Mass_temperature")

num_species, TIME, Density_max, Temperature_min, Peak_Abundances, MassFraction, TOTMASS_n, TOTMASS_T, A, B, C, D = result



#### Sorting out the ordering ####

if num_species == 3:

    fig1, ax = plt.subplots()

    ax.semilogy(TIME, Peak_Abundances[0, :], marker='None', linestyle='-', label="H2")

    ax.semilogy(TIME, Peak_Abundances[2, :], marker='None', linestyle='-.', label="CO")

    """
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

    ax.set_xlabel(r'$Time\ [Myr]$')

    ax.set_ylabel(r'$Peak\ Abundances$')

    Peak_Save = SaveLocation + "Peak Abundances.png"

    fig1.savefig(Peak_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


    fig2, ax = plt.subplots()


    ax.semilogy(TIME, MassFraction, marker='None', linestyle='-.', label="CO")

    """
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

    ax.set_xlabel(r'$Time\ [Myr]$')

    ax.set_ylabel(r'$Mass\ Fraction$')

    Fraction_Save = SaveLocation + "Mass Fraction.png"

    fig2.savefig(Fraction_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


if num_species == 9:


    fig1, ax = plt.subplots()

    ax.semilogy(TIME, Peak_Abundances[0, :], marker='None', linestyle='-', label="H2")

    ax.semilogy(TIME, Peak_Abundances[5, :], marker='None', linestyle='-.', label="CO")

    """
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

    ax.set_xlabel(r'$Time\ [Myr]$')

    ax.set_ylabel(r'$Peak\ Abundances$')

    Peak_Save = SaveLocation + "Peak Abundances.png"

    fig1.savefig(Peak_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


    fig2, ax = plt.subplots()


    ax.semilogy(TIME, MassFraction, marker='None', linestyle='-.', label="CO")

    """
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

    ax.set_xlabel(r'$Time\ [Myr]$')

    ax.set_ylabel(r'$Mass\ Fraction$')

    Fraction_Save = SaveLocation + "Mass Fraction.png"

    fig2.savefig(Fraction_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


T_Half = (TIME[-1] - TIME[0])/2


fig1, ax1 = plt.subplots()

ax1.semilogy(TIME, Density_max, color='b', marker='None', linestyle='-')

ax1.set_xlabel(r'$Time\ [Myr]$')

ax1.set_ylabel(r'$ Density\ \left[ cm^{-3} \right]$')

#ax1.tick_params('y', colors='b')

ax1.text(T_Half, 0.5e4, r'$max\ \rho$')


ax2 = ax1.twinx()

ax2.semilogy(TIME, Temperature_min, color='r', marker='None', linestyle='--')

ax2.set_ylabel('Temperature [K]')

#ax2.tick_params('y', colors='r')

ax2.text(T_Half, 10, r'$min\ T$')

fig1.tight_layout()

DandT_Save = SaveLocation + "Density and Temperature Plot.png"

fig1.savefig(DandT_Save, dpi=300, format='png', bbox_inches='tight')




fig2, ax = plt.subplots()

ax.semilogy(TIME, TOTMASS_n[0, :], marker='None', linestyle='-', label=r"$n\ >\ 1x10^2\ cm^{-3}$")

ax.semilogy(TIME, TOTMASS_n[1, :], marker='None', linestyle='-', label=r"$n\ >\ 1x10^3\ cm^{-3}$")

ax.semilogy(TIME, TOTMASS_n[2, :], marker='None', linestyle='-', label=r"$n\ >\ 1x10^4\ cm^{-3}$")

ax.semilogy(TIME, TOTMASS_n[3, :], marker='None', linestyle='-', label=r"$n\ >\ 1x10^5\ cm^{-3}$")

"""
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

ax.set_xlabel(r'$Time\ [Myr]$')

ax.set_ylabel(r'$Mass\ [M_{\odot}]$')

Density_Save = SaveLocation + "Mass evolution of Gas with Density.png"

fig2.savefig(Density_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


fig3, ax = plt.subplots()

ax.semilogy(TIME, TOTMASS_T[0, :], marker='None', linestyle='-', label=r"$T\ <\ %s K$" % A)

ax.semilogy(TIME, TOTMASS_T[1, :], marker='None', linestyle='-', label=r"$%s K\ < T\ <\ %s K$" % (A, B))

ax.semilogy(TIME, TOTMASS_T[2, :], marker='None', linestyle='-', label=r"$%s K\ <\ T\ <\ %s K$" % (C, D))

ax.semilogy(TIME, TOTMASS_T[3, :], marker='None', linestyle='-', label=r"$T\ >\ %s K$" % D)

"""
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))"""

lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

ax.set_xlabel(r'$Time\ [Myr]$')

ax.set_ylabel(r'$Mass\ [M_{\odot}]$')

Temperature_Save = SaveLocation + "Mass evolution of Gas with Temperature.png"

fig2.savefig(Temperature_Save, dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

