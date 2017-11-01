#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import numpy as np
import scipy as sc
import scipy.constants as sp
import astropy.constants as ap
import math

G = ap.G.to('cm3 / (g s2)')

M_sol = ap.M_sun.to('g')

m_p = ap.m_p

dist_cu = 1e17  # ('cm')

a = G*M_sol

tu = np.sqrt((dist_cu**3)/a)

print tu



print 4*59.5996