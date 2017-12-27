#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import os
import fnmatch
import numpy as np
import scipy as sc
import scipy.constants as sp
import astropy.constants as ap
import astropy.units as au
import math


""" Units"""


G = ap.G.to('cm3 / (g s2)')

M_cu = 1.991e33 * au.g

print "Mass unit : ", M_cu

dist_cu = 9.9999998e16 * au.cm

print "Distance unit : ", dist_cu

time_cu = 2.7436898e12 * au.s

print "Time unit : ", time_cu

v_cu = dist_cu / time_cu

vkms = 1e5*au.cm/au.km

vel_cu = v_cu/vkms

print "Velocity unit : ", vel_cu

rho_cu = M_cu/dist_cu**3

print "Density unit : ", rho_cu


K_E = 87288.56

m = 1e4*(ap.M_sun.to('g')/M_cu)

V = np.sqrt(2*K_E)/m

print V


POS = np.zeros((10, 1))

for i in range(len(POS[:, 0])):

    POS[i, 0] = np.linspace(0, 1, 20)