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
import astropy.units as au
from astropy.units import imperial
import math

""" Units"""

G = ap.G.to('cm3 / (g s2)')

M_sol = ap.M_sun.to('g')

dist_cu = 1e17*au.cm

print("Distance unit : ", dist_cu)

a = G*M_sol

time_cu = np.sqrt((dist_cu**3)/a)

print("Time unit : ", time_cu)

v_u = dist_cu / time_cu

vkms = 1e5*au.cm/au.km

vel_cu = v_u/vkms

print("Velocity unit : ", vel_cu)


""" Calculator """


def CloudRadius(n, m_cloud):

    m_sol = ap.M_sun.to('g')

    m_p = ap.m_p.to('g')

    m_cloud_sol = m_cloud*m_sol

    r = ((3*m_cloud_sol)/(4*np.pi*1.4*m_p*n))**(1/3)

    distance_unit = 1e17  # cm

    radius = r/distance_unit

    return radius


def CloudVel(V, time_cu, dist_cu):

    Vu = dist_cu/time_cu

    Vel = (V*1e5)/Vu*(au.cm/au.s)

    return Vel


while True:

    parameter = raw_input("What to calculate or 0 to exit: ")

    if parameter == '0':

        break

    if parameter == 'radius':

        n = raw_input("What is the number density, n?: ")

        n_int = np.float64(n)

        m_c = raw_input("What is the mass of the cloud, m?: ")

        m_cint = np.float64(m_c)

        CR = CloudRadius(n_int, m_cint)

        print CR

        decision = raw_input("New calculation (y/n)?: ")

        if decision=='y' or decision=='yes':

            continue

        else:

            break

    if parameter == 'velocity':

        V = raw_input("What is the velocity (Kms-1), v?: ")

        Vint = np.float64(V)

        Vel = CloudVel(Vint, time_cu, dist_cu)

        print Vel

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'crossing-time-ext':

        V = raw_input("What is the velocity (Code_Units), v?: ")

        Vint = np.float64(V)

        L = raw_input("What is the radius of the cloud (Code_Units), r?: ")

        Lint = np.float64(L)

        t = raw_input("What is the seperation in terms of cloud radii (Code_Units), t?: ")

        tint = np.float64(t)

        T_Cross = (tint*Lint)/Vint

        print T_Cross

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'crossing-time-int':

        V_disp = raw_input("What is the velocity dispersion (Code_Units), v?: ")

        Vint_disp = np.float64(V_disp)

        R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

        Rint = np.float64(R)

        T_Cross = Rint/Vint_disp

        print T_Cross

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'optimal-distance':

        V = raw_input("What is the velocity (Code_Units), v?: ")

        Vint = np.float64(V)

        V_disp = raw_input("What is the velocity dispersion (Code_Units), v?: ")

        Vint_disp = np.float64(V_disp)

        R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

        Rint = np.float64(R)

        L = Vint*(Rint/Vint_disp)

        print L

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

