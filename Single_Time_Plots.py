#!/usr/bin/python
"""
Created on Wed Oct 04 10:57:34 2017

@author: Gerwyn
"""

from __future__ import division

import matplotlib.pyplot as plt

from Arepo_Read_Master import *

### THE MEAT OF THE MATTER ###

filecount = 1

scalecount = 1

thin_factor = 100

v_log = 1


while filecount > 0:

    inputdir = raw_input("Enter directory: ")

    print '\nAvailable files:\n'

    print os.listdir(str(inputdir))

    os_dir = os.listdir(str(inputdir))

    inputfile = raw_input("Enter filename: ")

    file_dir = inputdir+inputfile

    DATA = aread(file_dir)

    POS = DATA[0]

    VEL = DATA[1]

    ID = DATA[2]

    MASS = DATA[3]

    U = DATA[4]

    RHO = DATA[5]

    POT = DATA[6]

    DIVV = DATA[7]

    ACCE = DATA[8]

    DUST = DATA[9]

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

    POS = POS[:ngas, :]

    MASS = MASS[:ngas]

    ### TEMP ###

    ABHE = 0.1

    n = RHO * rho_cu / ((1.0 + 4.0 * ABHE) * m_p)

    Energy = (U * RHO * ergs_cu) / (dist_cu**3)

    ntot = (1.0 + ABHE - CHEM[0, :] + CHEM[1, :]) * n

    TEMP = 2.0 * Energy / (3.0 * ntot * k_B)

    # print np.min(RHO)
    # print np.max(RHO)
    #
    # print np.min(n)
    # print np.max(n)
    #
    # print np.min(ntot)
    # print np.max(ntot)
    #
    # print np.amin(U)
    # print np.amax(U)
    #
    # print np.min(energy)
    # print np.max(energy)
    #
    # print np.min(TEMP)
    # print np.max(TEMP)


    if scalecount == 1:

        x_min = np.min(POS[:, 0])

        x_max = np.max(POS[:, 0])

        y_min = np.min(POS[:, 1])

        y_max = np.max(POS[:, 1])

        z_min = np.min(POS[:, 2])

        z_max = np.max(POS[:, 2])

    while scalecount > 0:

        parameter = raw_input("Enter data tag to plot (default = RHO; options DUST, CHEM, MASS, TEMP, POS or 0 to change file or exit) : ")

        if parameter == '0':

            break

        if x_max - x_min > y_max - y_min:

            plt.figure(figsize=(6 * (x_max - x_min) / (y_max - y_min), 5))

        elif x_max - x_min < y_max - y_min:

            plt.figure(figsize=(6, 5 * (y_max - y_min) / (x_max - x_min)))

        else:

            plt.figure(figsize=(6, 5))

        if parameter == '' or parameter == 'rho' or parameter == 'RHO':

            parameter = 'RHO'

            v_min, v_max = np.min(RHO), np.max(RHO)

            POSRHO = np.hstack((POS, RHO.reshape(RHO.size, 1)))

            POSRHOX = POSRHO[np.argsort(POSRHO[:, 3])]

            if v_log == 1:

                plt.scatter(POSRHOX[:, 0], POSRHOX[:, 1], c=np.log10(POSRHOX[:, 3]),
                            vmin=np.log10(v_min), vmax=np.log10(v_max), marker='.', s=1, cmap='jet')

                plt.colorbar(label=r"$log(density/code-units)$")

                if nsink > 0:
                    plt.scatter(SINKPOS[:, 0], SINKPOS[:, 1], c='k', marker='*', s=3)

            else:

                plt.scatter(POSRHOX[:, 0], POSRHOX[:, 1], c=POSRHOX[:, 3],
                            vmin=v_min, vmax=v_max, marker='.', s=1, cmap='jet')

                plt.colorbar(label=r"$density/code-units$")

                if nsink > 0:
                    plt.scatter(SINKPOS[:, 0], SINKPOS[:, 1], c='k', marker='*', s=3)

            plt.xlim((x_min, x_max))

            plt.ylim((y_min, y_max))

            plt.title('Max. density along z axis @ Unit time:%5.2f' % time)

            plt.xlabel('x')

            plt.ylabel('y', rotation='horizontal')

            plt.show()

            xy_decision = raw_input("Change scale? : ")

            if xy_decision == 'y':

                scalecount = 2

                x_min = input("xmin: ")

                x_max = input("xmax: ")

                y_min = input("ymin: ")

                y_max = input("ymax: ")

                if parameter == 'POS':

                    z_min = input("zmin: ")

                    z_max = input("zmax: ")

                    thin_factor = input("Data thin factor: ")

                if parameter == 'RHO':

                    v_log = input("Colour scale (0 = linear, 1 = log10) : ")

                continue

            else:

                break

        elif parameter == 'dust' or parameter == 'DUST':

            v_min, v_max = np.min(DUST), np.max(DUST)

            POSDUST = np.hstack((POS, DUST.reshape(DUST.size, 1)))

            POSDUSTX = POSDUST[np.argsort(POSDUST[:, 3])]

            if v_log == 1:

                plt.scatter(POSDUSTX[:, 0], POSDUSTX[:, 1], c=np.log10(POSDUSTX[:, 3]),
                            vmin=np.log10(v_min), vmax=np.log10(v_max), marker='.', s=1, cmap='jet')

                plt.colorbar(label=r"$log(Temp \/ /K)$")

                if nsink > 0:
                    plt.scatter(SINKPOS[:, 0], SINKPOS[:, 1], c='k', marker='*', s=3)

            else:

                plt.scatter(POSDUSTX[:, 0], POSDUSTX[:, 1], c=POSDUSTX[:, 3],
                            vmin=v_min, vmax=v_max, marker='.', s=1, cmap='jet')

                plt.colorbar(label=r"$Temp \/ /K$")

                if nsink > 0:
                    plt.scatter(SINKPOS[:, 0], SINKPOS[:, 1], c='k', marker='*', s=3)

            plt.xlim((x_min, x_max))

            plt.ylim((y_min, y_max))

            plt.title('Max. T along z axis @ Unit time:%5.2f' % time)

            plt.xlabel('x')

            plt.ylabel('y', rotation='horizontal')

            plt.show()

            xy_decision = raw_input("Change scale? : ")

            if xy_decision == 'y':

                scalecount = 2

                x_min = input("xmin: ")

                x_max = input("xmax: ")

                y_min = input("ymin: ")

                y_max = input("ymax: ")

                v_log = input("Colour scale (0 = linear, 1 = log10) : ")

                continue

            else:

                break

        elif parameter == 'chem' or parameter == 'CHEM':

            print "Options for chem species", np.arange(nspecies)

            ### 0 for H2    1 for H+    2 for CO ###

            X = raw_input("Choose which chem species to plot (See above for options) : ")

            if v_log == 1:

                plt.scatter(np.log10(n.value), np.log10(TEMP.value), c=np.log10(CHEM[np.int(X), :] + 1e-4), marker='.', s=0.5,
                            cmap='jet_r')

                plt.colorbar(label=r"$log(chem)$")

            else:

                plt.scatter(n.value, TEMP.value, c=(CHEM[X, :] + 1e-4), marker='.', s=0.5,
                            cmap='jet_r')

                plt.colorbar(label=r"$chem$")

            plt.xlabel(r'$n \/ (cm^{-3})$')

            plt.ylabel(r'$Temp \/ (K)$')

            plt.show()

            xy_decision = raw_input("Change to linear scale? : ")

            if xy_decision == 'y':

                v_log = 0

                continue

            else:

                break

        elif parameter == 'mass' or parameter == 'MASS':

            if v_log == 1:

                plt.scatter(np.log10(n.value), np.log10(TEMP.value), c=MASS, marker='.', s=0.5, cmap='jet_r')

                plt.colorbar(label=r"$Mass$")

            else:

                plt.scatter(n.value, TEMP.value, c=MASS, marker='.', s=0.5, cmap='jet_r')

                plt.colorbar(label=r"$Mass$")

            plt.xlabel(r'$n \/ (cm^{-3})$')

            plt.ylabel(r'$Temp \/ (K)$')

            plt.show()

            xy_decision = raw_input("Change to linear scale? : ")

            if xy_decision == 'y':

                v_log = 0

                continue

            else:

                break

            # elif parameter == 'POS':
            #
            #     POS_count = 0
            #
            #     for i in range(int(len(POS)/thin_factor)):
            #
            #         i *= thin_factor
            #
            #         if POS[i,0]>=x_min and POS[i,0]<=x_max and POS[i,1]>=y_min and POS[i,1]<=y_max
            #                                        and POS[i,2]>=z_min and POS[i,2]<=z_max:
            #
            #             if POS_count == 0:
            #
            #                 POSPLOT = POS[i,:]
            #
            #                 POS_count = 1
            #
            #             else:
            #
            #                 POSPLOT = np.vstack((POSPLOT,POS[i,:]))
            #
            #
            #     APEX_DISTANCE = np.sqrt(POSPLOT[:,0]**2+POSPLOT[:,1]**2)-np.sqrt(x_min**2+y_min**2)
            #
            #     MARKER_SIZE = 100*(1-(APEX_DISTANCE/np.max(APEX_DISTANCE))**2)
            #
            #
            #     fig = plt.figure(figsize=(10,10))
            #
            #     ax = Axes3D(fig)
            #
            #     ax.scatter(POSPLOT[:,0],POSPLOT[:,1], POSPLOT[:,2], marker='.', s=MARKER_SIZE, depthshade=True)
            #
            #     ax.view_init(azim=225, elev=10)
            #
            #     ax.set_zlabel('z')

        else:

            print '\n       NOT YET DEFINED'

            continue

