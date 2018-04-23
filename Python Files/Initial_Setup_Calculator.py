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


""" Units"""


G = ap.G.to('cm3 / (g s2)')

print("Gravitational const : ", G)

M_cu = 1.991e33 * au.g

print("Mass unit : ", M_cu)

dist_cu = 9.9999998e16 * au.cm

print("Distance unit : ", dist_cu)

time_cu = 2.7436898e12 * au.s

print("Time unit : ", time_cu)

v_cu = dist_cu / time_cu

vkms = 1e5*au.cm/au.km

vel_cu = v_cu/vkms

print("Velocity unit : ", vel_cu)

rho_cu = M_cu/dist_cu**3

print("Density unit : ", rho_cu)

ergs_cu = (M_cu*dist_cu**2)/time_cu**2

print("Energy unit (ergs) : ", ergs_cu)

print("Energy density : ", ergs_cu/dist_cu**3)

k_B = ap.k_B.cgs*((au.cm**2)*au.g/(au.s**2))/au.erg

print "Boltzmann Const : ", k_B

m_p = ap.m_p.to('g')

print "Mass of Proton : ", m_p

""" Calculator """


def CloudRadius(n, m_cloud):

    m_sol = ap.M_sun.to('g')

    m_p = ap.m_p.to('g')

    m_cloud_sol = m_cloud*m_sol

    r = ((3*m_cloud_sol)/(4*np.pi*1.4*m_p*n))**(1/3)

    radius = r/(dist_cu.value)

    return radius


while True:

    print(" ")

    print("Calculation options : ")

    print(" ")

    print("Radius, Velocity, Time, Crossing-Time")

    print(" ")

    print("Optimal-Distance, Boxsize, No.Density, Density, Second-Cloud-Distance, Freefall")

    print(" ")

    parameter = raw_input("What to calculate or 0 to exit: ")

    if parameter == '0':

        break

    if parameter == 'radius' or parameter == 'Radius':

        Decision = raw_input("Real or Code?: ")

        if Decision == 'code' or Decision == 'Code':

            n = raw_input("What is the number density, n?: ")

            n_int = np.float64(n)

            m_c = raw_input("What is the mass of the cloud, m?: ")

            m_cint = np.float64(m_c)

            Rint = CloudRadius(n_int, m_cint)

            print(Rint)

        elif Decision == 'real' or Decision == 'Real':

            n = raw_input("What is the number density, n?: ")

            n_int = np.float64(n)

            m_c = raw_input("What is the mass of the cloud, m?: ")

            m_cint = np.float64(m_c)

            Rint = CloudRadius(n_int, m_cint)*dist_cu

            cm_or_pc = raw_input("Cm or Parsecs (cm or pc): ")

            print(Rint.to(cm_or_pc))

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'velocity' or parameter == 'Velocity':

        Decision = raw_input("Real or Code?: ")

        if Decision == 'code' or Decision == 'Code':

            V = raw_input("What is the velocity (Kms-1), v?: ")

            Vel = np.float(V)

            Vint = Vel/vel_cu*(au.km/au.s)

            print(Vint)

        if Decision == 'real' or Decision == 'Real':

            V = raw_input("What is the velocity (Code), v?: ")

            Vel = np.float(V)

            Vint = Vel * vel_cu

            print(Vint)

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'time' or parameter == 'Time':

        Decision = raw_input("Real or Code?: ")

        if Decision == 'code' or Decision == 'Code':

            year_or_second = raw_input("MYears, Years or Seconds?: ")

            if year_or_second == 'MY' or year_or_second == 'my' or year_or_second == 'Myears' or year_or_second == 'MYears' or year_or_second == 'myears':

                t = raw_input("What is the time (Myear), t?: ")

                t_year = np.float64(t)*au.megayear.to('s')*au.s

                tint = t_year / time_cu

                print(tint)

            if year_or_second == 'Y' or year_or_second == 'y' or year_or_second == 'Years' or year_or_second == 'years' or year_or_second == 'year':

                t = raw_input("What is the time (year), t?: ")

                t_year = np.float64(t)*au.year.to('s')*au.s

                tint = t_year / time_cu

                print(tint)

            if year_or_second == 'S' or year_or_second == 's' or year_or_second == 'Seconds' or year_or_second == 'seconds':

                t = raw_input("What is the time (s), t?: ")

                tint = np.float(t)*au.s / time_cu

                print tint

        elif Decision == 'real' or Decision == 'Real':


            t = raw_input("What is the time (code), t?: ")

            t_int = np.float(t)

            year_or_second = raw_input("MYears, Years or Seconds?: ")

            if year_or_second == 'MY' or year_or_second == 'my' or year_or_second == 'Myears' or year_or_second == 'MYears' or year_or_second == 'myears':

                t_real = (t_int * time_cu)/(au.megayear.to('s')*(au.s / au.megayear))

            if year_or_second == 'Y' or year_or_second == 'y' or year_or_second == 'Years' or year_or_second == 'years' or year_or_second == 'year':

                t_real = (t_int * time_cu)/(au.year.to('s')*(au.s / au.year))

            if year_or_second == 'S' or year_or_second == 's' or year_or_second == 'Seconds' or year_or_second == 'seconds':

                t_real = t_int * time_cu

            tint = t_real

            print(tint)

        else :

            continue

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'crossing-time' or parameter == 'Crossing-Time' or parameter == 'crossing time':

        Decision = raw_input("External or Internal?: ")

        if Decision == 'Ext' or Decision == 'ext' or Decision == 'external':

            Mirror = raw_input("Are the velocities of the clouds mirrored (y/n)?: ")

            if Mirror == 'y' or Mirror == 'yes':

                try:

                    V1 = Vint
                    V2 = -Vint

                except NameError:

                    V = raw_input("What is the velocity (Code_Units), v?: ")

                    V1 = np.float(V)
                    V2 = -np.float(V)

            else:

                V_1 = raw_input("What is the velocity of Cloud 1 (Code_Units), v?: ")

                V_2 = raw_input("What is the velocity of Cloud 2 (Code_Units), v?: ")

                V1 = np.float(V_1)

                V2 = np.float(V_2)

            try:

                dint = L

            except NameError:

                d = raw_input("What is the seperation in terms of cloud radii (Code_Units), d?: ")

                dint = np.float(d)

            dV = V1 - V2

            T_Cross = dint/dV

            print(T_Cross)

        elif Decision == 'Int' or Decision == 'int' or Decision == 'internal':

            try:

                Vint_disp

            except NameError:

                V_disp = raw_input("What is the velocity dispersion (Code_Units), v?: ")

                Vint_disp = np.float(V_disp)

            try:

                Rint

            except NameError:

                R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

                Rint = np.float(R)

            T_Cross = Rint / Vint_disp

            print(T_Cross)

        else :

            continue

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'Freefall' or parameter == 'freefall-time' or parameter == 'freefall':

        try:

            Rho = rho

        except NameError:

            try:

                Rint

            except NameError:

                R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

                Rint = np.float(R)

            try:

                M_cint

            except NameError:

                M_c = raw_input("What is the mass of the cloud, m?: ")

                M_cint = np.float64(m_c)

            R_u = Rint * dist_cu

            M_u = M_cint * M_cu

            rho = M_u / ((4 / 3) * np.pi * R_u ** 3)


        Tff = np.sqrt((3*np.pi)/(32*G*rho))

        print(Tff)

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'optimal-distance' or parameter == 'Optimal-Distance' or parameter == 'optimal distance':

        Mirror = raw_input("Are the velocities of the clouds mirrored (y/n)?: ")

        if Mirror == 'y' or Mirror == 'yes':

            try:

                V1 = Vint
                V2 = -Vint

            except NameError:

                V = raw_input("What is the velocity (Code_Units), v?: ")

                V1 = np.float(V)
                V2 = -np.float(V)

        else:

            V_1 = raw_input("What is the velocity of Cloud 1 (Code_Units), v?: ")

            V_2 = raw_input("What is the velocity of Cloud 2 (Code_Units), v?: ")

            V1 = np.float(V_1)

            V2 = np.float(V_2)

        try:

            Vint_disp

        except NameError:

            V_disp = raw_input("What is the velocity dispersion (Code_Units), v?: ")

            Vint_disp = np.float(V_disp)

        try:

            Rint

        except NameError:

            R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

            Rint = np.float(R)

        Tc_ratio = raw_input("How much of crossing time (Percent or Fraction)?: ")

        Tcross_ratio = np.float(Tc_ratio)

        if Tcross_ratio > 1:

            Tcross_ratio = Tcross_ratio/100

        T_Cross = (Rint / Vint_disp)

        dV = V1 - V2

        L = dV*T_Cross*Tcross_ratio

        if L <= Rint:

            print("Less than Radius    ")

            code_or_pc = raw_input("Code or Parsecs (_ or pc): ")

            if code_or_pc == '' or code_or_pc == '_':

                print(Rint)

            else :

                Rint = Rint*dist_cu

                print(Rint.to(code_or_pc))

        else:

            code_or_pc = raw_input("Code or Parsecs (_ or pc): ")

            if code_or_pc == '' or code_or_pc == '_':

                print(L)

            else :

                L = L*dist_cu

                print(L.to(code_or_pc))

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'boxsize' or parameter == 'Boxsize':

        try:

            Rint

        except NameError:

            R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

            Rint = np.float(R)

        try:
            Dint = L

        except NameError:

            D = raw_input("What is the optimal distance between the two clouds (Code_Units), d?: ")

            Dint = np.float(D)

        boxsize = 4*Rint + Dint

        print("Boxsize =", boxsize)

        print("Y,Z =", boxsize/2)

        print("X1 =", Rint*2)

        print("X2 =", Rint*2 + Dint)

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'no.density' or parameter == 'No.Density':

        decision = raw_input("No.Density known ?: ")

        ABHE = 0.1

        conversion = ((1.0 + 4.0 * ABHE) * ap.m_p.to('g'))

        if decision == 'y' or decision == 'yes':

            decision1 = raw_input("Code or Real ?: ")

            if decision1 == 'Code' or decision1 == 'code':

                nd = raw_input("What is the number density (Code_Units), n?: ")

                N_int = np.float(nd)*(au.cm**-3)*(dist_cu**-3)

                print(N_int)

            else:

                nd = raw_input("What is the number density (cm-3), n?: ")

                N_int = np.float(nd)/((au.cm**-3)*(dist_cu**-3))

                print(N_int)

        else:

            d_int = raw_input("What is the number density (g/cm3)?: ")

            d_unit = au.g*(au.cm**-3)

            D_int = (np.float(d_int)*d_unit)/conversion

            print(D_int)

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'density' or parameter == 'Density':

        try:

            Rint

        except NameError:

            R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

            Rint = np.float(R)

        try:

            M_cint

        except NameError:

            M_c = raw_input("What is the mass of the cloud, m?: ")

            M_cint = np.float64(M_c)

        R_u = Rint*dist_cu

        M_u = M_cint*M_cu

        rho = M_u/((4/3)*np.pi*R_u**3)

        print(rho)

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'Jeans-Mass' or parameter == 'jeans mass':

        try:

            rho

        except NameError:

            try:

                Rint

            except NameError:

                R = raw_input("What is the radius of the cloud (Code_Units), r?: ")

                Rint = np.float(R)

            try:

                M_cint

            except NameError:

                M_c = raw_input("What is the mass of the cloud, m?: ")

                M_cint = np.float64(M_c)

            R_u = Rint*dist_cu

            M_u = M_cint*M_cu

            rho = M_u/((4/3)*np.pi*R_u**3)

        T = raw_input("What is the temperature of the cloud, K?: ")

        T_int = np.float64(T)*au.K

        Cs = ((k_B*T_int)/(1.4*m_p))**(1 / 2)

        Mj = (4*np.pi*Cs**3)/(3*((G**3)*rho)**(1/2))

        print(Mj/M_cu)

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break

    if parameter == 'Sound-Speed' or parameter == 'sound speed':

        T = raw_input("What is the temperature of the cloud, K?: ")

        T_int = np.float64(T)*au.K

        Cs = ((k_B*T_int)/(2.7*m_p))**(1/2)

        print(Cs.to('km/s'))

        decision = raw_input("New calculation (y/n)?: ")

        if decision == 'y' or decision == 'yes':

            continue

        else:

            break
