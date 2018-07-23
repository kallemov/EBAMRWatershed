import sys
import numpy as np
#import matplotlib.pyplot as plt
from scipy import ndimage, misc

if (len(sys.argv)<3):
    sys.exit("Usage: averageDEM   input_file_name   average_size   type(0-uniform, 1-Gaussian)")
try:
    demfile = open(sys.argv[1],'r')
except:
    sys.exit("Cannot read the file ", sys.argv[1])
smoothN=int(sys.argv[2])

try:
    method = int(sys.argv[3])
except:
    method=0

##demfile = open('copper_stamp_zero.dem','r')
toplines = filter(lambda x: str(x).isdigit(),demfile.readline().split())
res = list(map(int,toplines))
spacing=[]
for x in demfile.readline().split():
    try:
        float(x)
    except:
        continue
    else:
        spacing.append(float(x))

print(res, spacing)
demvalues=np.empty((res[0],res[1]))
counter=0

for line in demfile:
    for x in line.split():
        demvalues[int(counter%res[0])][int(counter/res[0])]=float(x)
        counter += 1

#plt.subplot(211)
#plt.imshow(demvalues)


## need to choose mode from reflect, mirror,nearest (skipping unfit ones)
if method==0:
    smoothed = ndimage.uniform_filter(demvalues, size=smoothN, mode='nearest')
else:
    smoothed = ndimage.gaussian_filter(demvalues, sigma=smoothN, mode='nearest')

smoothedDEMfile = open(sys.argv[1].replace(".dem","")+"_smoothed_"+str(smoothN)+".dem",'w')
smoothedDEMfile.write("dem_ncells "+str(res[0])+' '+str(res[1])+'\n')
smoothedDEMfile.write("dem_grid_spacing "+str(spacing[0])+' '+str(spacing[1])+'\n')
np.savetxt(smoothedDEMfile,np.transpose(smoothed),fmt='%.5f')

#plt.subplot(212)
#plt.imshow(smoothed)
#plt.show()

