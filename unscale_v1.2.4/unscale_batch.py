import numpy as np
import scipy.interpolate
import math
import matplotlib.pyplot as plt
import subprocess
import sys
from astropy.io import fits
from subprocess import Popen, PIPE
#%matplotlib inline

deftype = int(sys.argv[1])
if(deftype==1):
    our_file = fits.open('Trf_Johannsen_a13.fits')
elif(deftype==2):
    our_file = fits.open('Trf_Johannsen_a22.fits')
elif(deftype==3):
    our_file = fits.open('Trf_Johannsen_e3.fits') 
elif(deftype==11):
    our_file = fits.open('Trf_KRZ_d1.fits') # THIS IS FOR KRZ d1  
elif(deftype==12):
    our_file = fits.open('Trf_KRZ_d2.fits') # THIS IS FOR KRZ d1  
elif(deftype==13):
    our_file = fits.open('Trf_KRZ_d3.fits') # THIS IS FOR KRZ d1  
elif(deftype==14):
    our_file = fits.open('Trf_KRZ_d4.fits') # THIS IS FOR KRZ d1  
elif(deftype==15):
    our_file = fits.open('Trf_KRZ_d5.fits') # THIS IS FOR KRZ d1  
elif(deftype==16):
    our_file = fits.open('Trf_KRZ_d6.fits') # THIS IS FOR KRZ d1  
else:
    print("deftype doesn't match with any deformation parameter")


spindpgrid = our_file[1].data
#spin, dpscaled, redchi = np.loadtxt('output_relxill_raw.dat', unpack=True)

filein = sys.argv[2]
fileout = sys.argv[3]

chi, redchi, spind, spin, defind, defpar= np.loadtxt(filein, unpack=True)

b = open(filein, "r").readlines()
aind = int(max(spind))+1
dpind = int(max(defind))+1
redchi = np.zeros([aind,dpind])
spin = np.zeros([aind,dpind])
dpscaled = np.zeros([aind,dpind])
dpcorrect = np.zeros([aind,dpind])
for line in b:
    spin[int(line.split()[2]),int(line.split()[4])] = line.split()[3]
    dpscaled[int(line.split()[2]),int(line.split()[4])] = line.split()[5]
    redchi[int(line.split()[2]),int(line.split()[4])] = line.split()[1]

#print dpscaled[0][0]

for i in np.arange(aind):
    for j in np.arange(29):
        if((spin[i][0]<spindpgrid[j+1][0]) and (spin[i][0]>spindpgrid[j][0])):
            #print spin[i][0],spindpgrid[j][0], spindpgrid[j+1][0]
            deflim = 0.0
            ifaca = (spin[i][0]-spindpgrid[j][0])/(spindpgrid[j+1][0]-spindpgrid[j][0])
            for k in np.arange(dpind):
                if(dpscaled[i][k] < 0):
                    deflim = spindpgrid[j][1][0] + ifaca*(spindpgrid[j+1][1][0]-spindpgrid[j][1][0])
                    #print(spindpgrid[j][0])
                    dpcorrect[i][k] = -dpscaled[i][k]*deflim
                else:
                    deflim = spindpgrid[j][1][29] + ifaca*(spindpgrid[j+1][1][29]-spindpgrid[j][1][29])
                    dpcorrect[i][k] = dpscaled[i][k]*deflim
                
f = open(fileout,"w")
for i in np.arange(aind):
    for j in np.arange(dpind):
        f.write(" " + str(spin[i][j]) + "\t" + str(dpcorrect[i][j]) + "\t" + str(redchi[i][j]) + "\n" )
        #f.write(" " + str(spin[i*dpind+j]) + "\t" + str(dpcorrect[i*dpind+j]) + "\t" + str(redchi[i*dpind+j]) + "\n" )
    f.write("\n")
    
f.close()
