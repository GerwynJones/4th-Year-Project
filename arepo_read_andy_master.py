#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

k_B =1.3806581e-16
mp = 1.6726575e-24

def aread(filename):
        def readtag (fobj):
                tag = np.fromfile(fobj,dtype="|S4",count=1)
                return tag

        def readlen (fobj):
                skip = np.fromfile(fobj,dtype=np.int32,count=1)
                # recordlength is length in bits
                recordunit = np.fromfile(fobj,dtype=np.int32,count=1)
                numbits = np.fromfile(fobj,dtype=np.int32,count=1)
                recordlength = (numbits/recordunit)
                return recordlength
        
        print "Reading file:", filename

        bytesleft= 256 - 6*4 - 6*8 - 2*8 - 2*4 -6*4 - 2*4 - 4*8 - 2*4 - 6*2 - 3*4 - 4

        with open(filename,"rb") as file:
                ##### Header
                dummy = np.fromfile(file,dtype=np.int32,count=1)
                tag = np.fromfile(file,dtype="|S4",count=1)
                dummy = np.fromfile(file,dtype=np.int32,count=3)
                npart = np.fromfile(file, dtype=np.int32, count=6)
                massarr = np.fromfile(file, dtype=np.double,count=6)
                time = np.fromfile(file,dtype=np.double,count=1)##[0]
                redshift = np.fromfile(file,dtype=np.double,count=1)[0]
                flag_sfr = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_feedback = np.fromfile(file,dtype=np.int32,count=1)[0]
                npartTotal = np.fromfile(file,dtype=np.int32,count=6)
                flag_cooling = np.fromfile(file,dtype=np.int32,count=1)[0]
                num_files = np.fromfile(file,dtype=np.int32,count=1)[0]
                boxsize = np.fromfile(file,dtype=np.double,count=1)##[0]
                cos1 = np.fromfile(file,dtype=np.double,count=1)[0]
                cos2 = np.fromfile(file,dtype=np.double,count=1)[0]
                hubble_param = np.fromfile(file,dtype=np.double,count=1)[0]
                flag_stellarage = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_metals = np.fromfile(file,dtype=np.int32,count=1)[0]
                npartHighword = np.fromfile(file,dtype='int16',count=6)
                flag_entropy = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_dp = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_1pt = np.fromfile(file,dtype=np.int32,count=1)[0]
                scalefactor = np.fromfile(file,dtype='float32',count=1)[0]
                dum = np.fromfile(file,dtype='int16',count=bytesleft/2)
                pad = np.fromfile(file,dtype='int16',count=4)

                print "npart array:", npart
                print "mass array:", massarr
                print "time in codeunits:", time
                print "Total npart:", npartTotal

                print "Header finished\n"

                N = int(sum(npart))
                ngas = npart[0]
                nsink = npart[5]
                tag = readtag(file)
                recordlength = readlen(file)
                print "Found array with tag ",tag, "and record length", recordlength

                pos = np.zeros((N,3),dtype=np.double)
                vel = np.zeros((N,3),dtype=np.double)
                partid = np.zeros(N,dtype=np.int32)
                mass = np.zeros(N,dtype=np.double)
                u = np.zeros(ngas,dtype=np.double)
                rho = np.zeros(ngas,dtype=np.double)
                potential = np.zeros(N,dtype=np.double)
                divv = np.zeros(ngas,dtype=np.double)
                accel = np.zeros((N,3),dtype=np.double)
                tdust = np.zeros(ngas,dtype=np.double)
                dtime = np.zeros(ngas,dtype=np.double)
                bfield = np.zeros((ngas,3),dtype=np.double)
                divb = np.zeros(ngas,dtype=np.double)
                softening = np.zeros(ngas,dtype=np.double)
                chem = np.zeros((ngas,3),dtype=np.double)   
                
                while len(recordlength > 0):
                        #sys.exit("Stop")
                        ## Case statement
                        if(tag=="POS "):
                                print "Reading positions"
                                posdum = np.fromfile(file,dtype=np.double,count=3*N)
                                print "Read positions, min and max: ",np.amin(posdum),np.amax(posdum)
                                pos = np.reshape(posdum,(-1,3))
                        elif(tag=="VEL "):
                                print "Reading velocities"
                                veldum = np.fromfile(file,dtype=np.double,count=3*N)
                                print "Read velocities, min and max: ",np.amin(veldum),np.amax(veldum)
                                vel = np.reshape(veldum,(-1,3))
                        elif(tag=="ID  "):
                                print "Reading particle IDs"
                                partid = np.fromfile(file,dtype=np.int32,count=N)
                                print "Read IDs, min and max:", np.amin(partid),np.amax(partid)
                        elif(tag=="MASS"):
                                print "Reading particle masses"
                                mass = np.fromfile(file,dtype=np.double,count=N)
                                print "Read masses, min and max: ",np.amin(mass),np.amax(mass)
                        elif(tag=="U   "):
                                print "Reading u"
                                u = np.fromfile(file,dtype=np.double,count=ngas)
                                print "Read u, min and max: ",np.amin(u),np.amax(u)
                        elif(tag=="RHO "):
                                print "Reading densities"
                                rho = np.fromfile(file,dtype=np.double,count=ngas)
                                print "Read densities, min and max: ",np.amin(rho),np.amax(rho)
                        elif(tag=="POT "):
                                print "Reading potentials"
                                potential = np.fromfile(file,dtype=np.double,count=N)
                                print "Read potentials, min and max: ",np.amin(potential),np.amax(potential)
                        elif(tag=="DIVV"):
                                print "Reading velocity divergence"
                                divv = np.fromfile(file,dtype=np.double,count=ngas)
                                print "Read velocity divergence, min and max:",np.amin(divv),np.amax(divv)
                                u = np.fromfile(file,dtype=np.double,count=ngas)
                                print "Read u, min and max: ",np.amin(u),np.amax(u)
                        elif(tag=="ACCE"):
                                print "Reading accelerations"
                                acceldum = np.fromfile(file,dtype=np.double,count=N*3)
                                print "Read accelerations, min and max: ",np.amin(accel),np.amax(accel)
                                accel = np.reshape(acceldum,(-1,3))
                        elif(tag=="DUST"):
                                print "Reading dust temperatures"
                                tdust = np.fromfile(file,dtype=np.double,count=ngas)
                                print "Read dust temperatures, min and max: ", np.amin(tdust),np.amax(tdust)
                        elif(tag=="TSTP"):
                                print "Reading timesteps"
                                dtime = np.fromfile(file,dtype=np.double,count=N)
                                print "Read timesteps, min and max:",np.amin(dtime),np.amax(dtime)
                        elif(tag=="BFLD"):
                                print "Reading magnetic field"
                                bfielddum = np.fromfile(file,dtype=np.double,count=3*ngas)
                                print "Read magnetic field, min and max: ",np.amin(bfield),np.amax(bfield)
                                bfield = np.reshape(bfielddum,(-1,3))
                        elif(tag=="DIVB"):
                                print "Reading magnetic field divergence"
                                divb = np.fromfile(file,dtype=np.double,count=ngas)
                                print "Read magnetic field divergence, min and max: ",np.amin(divb),np.amax(divb)
                        elif(tag=="SOFT"):
                                print "Reading softening"
                                softening = np.fromfile(file,dtype=np.double,count=N)
                                print "Read softenings, min and max: ",np.amin(softening),np.amax(softening)
                        elif(tag=="CHEM"):
                                print "Reading chemistry"
                                chemdum = np.fromfile(file,dtype=np.double,count=3*ngas)
                                print "Read chemistry, min and max: ",np.amin(chem),np.amax(chem)
                                chem = np.reshape(chemdum,(-1,3))
                        else:
                                print "Skipping through property",tag," with record length", recordlength
                                dummy = np.fromfile(file,dtype=np.double,count=recordlength)[0]
                        
                        pad = np.fromfile(file,dtype='int16',count=4)
                        tag = readtag(file)
                        recordlength = readlen(file)
                        if(len(recordlength) > 0):
                                print "Found array with tag ",tag, "and record length", recordlength
        
        if(nsink>0):
            print "Reading sinks"
            idsink = np.linspace(0,nsink-1,nsink,dtype='int16') + ngas
            sinkpos = pos[idsink,:]
            sinkvel = vel[idsink,:]
            sinkmass = mass[idsink]
            sinkid = partid[idsink]
        else:
            sinkpos, sinkvel, sinkmass, sinkid = 0, 0, 0, 0
     
        print "Finished reading file:", filename,"\n" 
        
	return pos, vel, partid, mass, u, rho, potential, divv, accel, tdust, dtime, bfield, divb, softening, chem, time, nsink, sinkpos, sinkvel, sinkmass, sinkid
    


###THE MEAT OF THE MATTER

filecount = 1
scalecount = 1
thin_factor = 100
v_log = 1

while filecount>0:
    print '\nAvailable files:\n'
    print os.listdir('.')
    inputfile = raw_input("Enter filename: ")
    
    DATA = aread(inputfile)
    
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
    SINKPOS = DATA[17]
    SINKVEL = DATA[18]
    SINKMASS = DATA[19]
    SINKID = DATA[20]
    
    POS = POS[:len(POS)-nsink,:]
    
    if scalecount==1:
        x_min = np.min(POS[:,0])
        x_max = np.max(POS[:,0])
        y_min = np.min(POS[:,1])
        y_max = np.max(POS[:,1])
        z_min = np.min(POS[:,2])
        z_max = np.max(POS[:,2])
    
    while scalecount >0:
        parameter = raw_input("Enter data tag to plot (default = RHO; options DUST, IMF or 0 to change file or exit): ")
        if parameter == '0':
            break
        
        if x_max - x_min > y_max - y_min:
            plt.figure(figsize=(6*(x_max - x_min)/(y_max - y_min),5))
        elif x_max - x_min < y_max - y_min:
            plt.figure(figsize=(6,5*(y_max - y_min)/(x_max - x_min)))
        else:
            plt.figure(figsize=(6,5))        
                    
        if parameter == '' or parameter == 'RHO':
            parameter = 'RHO'
            v_min, v_max = np.min(RHO), np.max(RHO)
            POSRHO = np.hstack((POS,RHO.reshape(RHO.size,1)))
            POSRHOX = POSRHO[np.argsort(POSRHO[:,3])]
            if v_log == 1:
                plt.scatter(POSRHOX[:,0],POSRHOX[:,1], c=np.log10(POSRHOX[:,3]), vmin = np.log10(v_min), vmax = np.log10(v_max), marker='.', s=1, cmap = 'jet')
                plt.colorbar(label="log(density/code-units)")
                if(nsink>0):
                    plt.scatter(SINKPOS[:,0],SINKPOS[:,1], c='k', marker='*', s=10)
            else:
                plt.scatter(POSRHOX[:,0],POSRHOX[:,1], c=POSRHOX[:,3], vmin = v_min, vmax = v_max, marker='.', s=1, cmap = 'jet')
                plt.colorbar(label="density/code-units")
                if(nsink>0):
                    plt.scatter(SINKPOS[:,0],SINKPOS[:,1], c='k', marker='*', s=10)
            plt.xlim((x_min,x_max))
            plt.ylim((y_min,y_max))
            plt.title('Max. density along z axis @ Unit time:%5.2f' %(time))
            plt.xlabel('x')
            plt.ylabel('y',rotation='horizontal')
            plt.show()

            xy_decision = raw_input("Change scale?: ")
            
            if xy_decision == 'y':
                scalecount = 2
                x_min = input("xmin: ")
                x_max = input("xmax: ")
                y_min = input("ymin: ")
                y_max = input("ymax: ")
                if parameter=='POS':
                    z_min = input("zmin: ")
                    z_max = input("zmax: ")
                    thin_factor = input("Data thin factor: ")
                if parameter=='RHO':
                    v_log = input("Colour scale (0 = linear, 1 = log10): ")
                continue
            else:
                break


        elif parameter == 'DUST':
            v_min, v_max = np.min(DUST), np.max(DUST)
            POSDUST = np.hstack((POS,DUST.reshape(DUST.size,1)))
            POSDUSTX = POSDUST[np.argsort(POSDUST[:,3])]
            if v_log == 1:
                plt.scatter(POSDUSTX[:,0],POSDUSTX[:,1], c=np.log10(POSDUSTX[:,3]), vmin = np.log10(v_min), vmax = np.log10(v_max), marker='.', s=1, cmap = 'jet')
                plt.colorbar(label="log(Temp./K)")
                if(nsink>0):
                    plt.scatter(SINKPOS[:,0],SINKPOS[:,1], c='k', marker='*', s=10)
            else:
                plt.scatter(POSDUSTX[:,0],POSDUSTX[:,1], c=POSDUSTX[:,3], vmin = v_min, vmax = v_max, marker='.', s=1, cmap = 'jet')
                plt.colorbar(label="Temp./K")
                if(nsink>0):
                    plt.scatter(SINKPOS[:,0],SINKPOS[:,1], c='k', marker='*', s=10)
            plt.xlim((x_min,x_max))
            plt.ylim((y_min,y_max))
            plt.title('Max. T along z axis @ Unit time:%5.2f' %(time))
            plt.xlabel('x')
            plt.ylabel('y',rotation='horizontal')
            plt.show()

            xy_decision = raw_input("Change scale?: ")
            
            if xy_decision == 'y':
                scalecount = 2
                x_min = input("xmin: ")
                x_max = input("xmax: ")
                y_min = input("ymin: ")
                y_max = input("ymax: ")
                if parameter=='POS':
                    z_min = input("zmin: ")
                    z_max = input("zmax: ")
                    thin_factor = input("Data thin factor: ")
                if parameter=='DUST':
                    v_log = input("Colour scale (0 = linear, 1 = log10): ")
                continue
            else:
                break





       
#        elif parameter == 'POS':
#            POS_count = 0
#            for i in range(int(len(POS)/thin_factor)):
#                i *= thin_factor
#                if POS[i,0]>=x_min and POS[i,0]<=x_max and POS[i,1]>=y_min and POS[i,1]<=y_max and POS[i,2]>=z_min and POS[i,2]<=z_max:
#                    if POS_count == 0:
#                        POSPLOT = POS[i,:]
#                        POS_count = 1
#                    else:
#                        POSPLOT = np.vstack((POSPLOT,POS[i,:]))
#            
#            APEX_DISTANCE = np.sqrt(POSPLOT[:,0]**2+POSPLOT[:,1]**2)-np.sqrt(x_min**2+y_min**2)
#            MARKER_SIZE = 100*(1-(APEX_DISTANCE/np.max(APEX_DISTANCE))**2)
#            
#            fig = plt.figure(figsize=(10,10))
#            ax = Axes3D(fig)
#            ax.scatter(POSPLOT[:,0],POSPLOT[:,1], POSPLOT[:,2], marker='.', s=MARKER_SIZE, depthshade=True)
#            ax.view_init(azim=225, elev=10)
#            ax.set_zlabel('z')
        
        elif parameter == 'IMF':
            plt.hist(SINKMASS, bins = np.logspace(-2,3,50))
            plt.xscale('log')
            plt.title('Sink mass distribution @ Unit time:%5.2f\n Total sinks:%d' %(time, nsink))
            plt.xlabel('Sink mass')
            plt.ylabel('Number')
            plt.show()
            continue
                    
        else:
            print '\n       NOT YET DEFINED'
            continue

    file_decision=raw_input("Make a new plot?: ")
    if file_decision=='y':
        continue
    else:
        break
