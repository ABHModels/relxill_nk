import numpy as np
from astropy.io import fits
from subprocess import Popen, PIPE
import sys

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
else:
    print("deftype doesn't match with any deformation parameter")
#our_file = fits.open('Trf_ordered_e_0.00e+00.a22_0.00e+00.a52_0.00e+00.fits')
spindpgrid = our_file[1].data
dpcorrect = 0

spin = float(sys.argv[2])
dpscaled = float(sys.argv[3])

for j in np.arange(29):
    if((spin<spindpgrid[j+1][0]) and (spin>=spindpgrid[j][0])):
            #print spin[i][0],spindpgrid[j][0], spindpgrid[j+1][0]
            #for k in np.arange(dpind):
            deflim=0.0
            ifaca = (spin-spindpgrid[j][0])/(spindpgrid[j+1][0]-spindpgrid[j][0])
            if(dpscaled < 0):
                deflim = spindpgrid[j][1][0] + ifaca*(spindpgrid[j+1][1][0]-spindpgrid[j][1][0])
                dpcorrect = -dpscaled*deflim
            else:
                deflim = spindpgrid[j][1][29] + ifaca*(spindpgrid[j+1][1][29]-spindpgrid[j][1][29])
                dpcorrect = dpscaled*deflim
                        
print(dpcorrect)
