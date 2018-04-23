#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division
import os

import numpy as np
import matplotlib.pyplot as plt

n = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

nmin = np.min(n)

nmax = np.max(n)

nbins = 50

dlogn = (np.log10(nmax) - np.log10(nmin))/nbins


nbinmin = nmin

Massbin = np.zeros(nbins)

nbincent = np.zeros(nbins)

for i in range(0, nbins-1):

    nbinmax = 10**((i+1)*dlogn)

    nbincent[i] = 10**((i+1)*(dlogn/2))

    ParticleID = np.argwhere((n > nbinmin) & (n < nbinmax))

    nbinmin = nbinmax

    print Massbin








