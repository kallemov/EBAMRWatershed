import sys
import numpy as np
#import matplotlib.pyplot as plt
from scipy import ndimage, misc

if (len(sys.argv)<3):
    sys.exit("Usage: coarseDEM   input_file_name   Xres   Yres")
try:
    demfile = open(sys.argv[1],'r')
except:
    sys.exit("Cannot read the file ", sys.argv[1])
Xres=int(sys.argv[2])
Yres=int(sys.argv[3])


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
print("performing coarsening")

coarsed = misc.imresize(demvalues,(Xres,Yres),interp='bicubic', mode='F')

coarsedDEMfile = open(sys.argv[1].replace(".dem","")+"_coarsed_"+str(Xres)+"_"+str(Yres)+".dem",'w')
coarsedDEMfile.write("dem_ncells "+str(Xres)+' '+str(Yres)+'\n')
coarsedDEMfile.write("dem_grid_spacing "+str(spacing[0]*res[0]/Xres)+' '+str(spacing[1]*res[1]/Yres)+'\n')

#plt.subplot(212)
#plt.imshow(coarsed)
#plt.show()
#print(coarsed)
np.savetxt(coarsedDEMfile,np.transpose(coarsed),fmt='%.5f')

