import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, misc

if (len(sys.argv)<2):
    sys.exit("need input file!")
try:
    vertexfile = open(sys.argv[1],'r')
except:
    sys.exit("Cannot read the file ", sys.argv[1])
import re
figc=0
while True:
    line=vertexfile.readline()
    if not line:break
    if "=>" not in line:
        continue
    templine=re.findall(r"[\d]+",line)
    origin=map(int,templine[0:3])
    if (int(templine[3])>0): 
        sys.exit("no support for multicells")
    nump=(templine[4])
    numpair=int(templine[5])
    #print origin, nump, numpair
    points=[]
    pairs=[]
    condition= True ##(origin[0]==46 and origin[1]==23)
    if (condition): 
        plt.figure(figc/9+1)
        ##plt.subplot(3,3,int(figc%9)+1)# figsize=(10,10))
        plt.subplot(1,2,int(figc%9)+1)# figsize=(10,10))
        plt.title(str(origin)) # +'  '+ str(nump)+'-'+str(numpair))
        plt.axis([origin[0]-0.2, origin[0]+1.2, origin[1]-0.2, origin[1]+1.2])
        ax=[origin[0],origin[0],origin[0]+1,origin[0]+1,origin[0]]
        ay=[origin[1],origin[1]+1,origin[1]+1,origin[1], origin[1]]
        plt.plot(ax,ay,'green',linestyle='dashed')
        plt.hold(True)

    #print origin, numpair

    for i in range(0,numpair):
        nextline=vertexfile.readline()
        templine=re.split('=|,|\)+|\(+', nextline)
        #print templine
        point1=map(float,templine[2:5])
        if point1 not in points:
            points.append(point1)
        point2=map(float,templine[7:10])
        if point2 not in points:
            points.append(point2)
        pairs.append((point1,point2))
        normal=map(float,templine[19:21])
        #lengths.append(float(templine[21]))
        grad1=map(float,templine[25:27])
        grad2=map(float,templine[30:32])
        #print templine[33:38]
        dotGrad=int(templine[34])


        x=[point1[0],point2[0]]
        y=[point1[1],point2[1]]
        #print '\t', x,y

        if (point1[2]==point2[2] and point1[2]==origin[2]):
            col='red'
        elif (point1[2]==point2[2] and point1[2]==origin[2]+1):
            col='blue'
        else:
            col='black'
            
        ls='-'
        if dotGrad<0:
            continue
            ls='--'
        elif dotGrad==0:
            ls=':'
            col='purple'
        if (condition):            
            plt.plot(x,y,color=col, marker='o', linestyle=ls)
        ##plt.scatter(x,y,color='yellow',marker='o')
       
        if dotGrad<=0:
            continue
        xv = 0.5*(point1[0]+point2[0])
        yv = 0.5*(point1[1]+point2[1])
        if (condition):            
            plt.quiver(xv,yv,normal[0],normal[1], scale=None)
           # plt.quiver(point1[0],point1[1],grad1[0],grad1[1], scale=None, color='brown')
           # plt.quiver(point2[0],point2[1],grad2[0],grad2[1], scale=None, color='brown')
    #plt.draw()    
    #plt.show()
    if (condition): 
        figc +=1
#    print len(pairs), len(points)
plt.show()

