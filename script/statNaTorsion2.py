#!/usr/bin/env python
docstring='''
statNaTorsion.py NaTorsion.raw.gz
    read torsion angle data, draw Torsion statistics
'''
import sys, os
import gzip

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import numpy as np

def statNaTorsion2(rawfile):
    data=[]
    pucker_list=[]
    N_list=[]
    angle_list=[]
    fp=gzip.open(rawfile)
    for line in fp.read().splitlines():
        if line.startswith("####"):
            continue
        if line.startswith("N c resi"):
            angle_list=line[9:].strip().split()
            continue
        N_list.append(line[0])
        #pucker=float(line[9:17])
        pucker=float(line[26:33])
        if pucker>0:
            pucker_list.append(1)
        elif pucker<0 and pucker>=-180:
            pucker_list.append(-1)
        else:
            pucker_list.append(0)
        data.append([float(t) for t in line[9:].split()])
    fp.close()
    data=np.array(data)
    pucker_array=np.array(pucker_list,dtype=int)
    N_array=np.array(N_list,dtype=str)
    
    #### set parameters for matplotlib ####
    fontsize=8
    matplotlib.rcParams["legend.borderaxespad"]=0.1
    matplotlib.rcParams["legend.borderpad"]    =0.1
    matplotlib.rcParams["legend.columnspacing"]=1
    matplotlib.rcParams["legend.handlelength"] =1
    matplotlib.rcParams["legend.handletextpad"]=0.1 
    matplotlib.rcParams["legend.labelspacing"] =0.1
    xticks=np.arange(-180,181,60)

    hpp='''#ifndef pseudo_torsion_hpp
#define pseudo_torsion_hpp 1
const int pseudo_torsion_cdf[360][360]={
'''
    
    #### pseudo torsion ####
    plt.figure(figsize=(7.87,7.87))
    plt.subplot(2,2,1)
    valid_residues=(data[:,7]>=-180)*(data[:,8]>=-180)
    eta_data=data[valid_residues,7]
    theta_data=data[valid_residues,8]
    data_fine=np.zeros((360,360))
    for eta,theta in zip(eta_data.tolist(),theta_data.tolist()):
        i=int(theta+180)
        j=int(eta+180)
        if i==data_fine.shape[0]:
            i=0
        if j==data_fine.shape[1]:
            j=0
        data_fine[i][j]+=1
    vmax=data_fine.max(axis=0).mean()
    plt.imshow(data_fine,vmin=0,vmax=vmax,cmap="jet",interpolation='nearest')
    plt.axis([-0.5,360-0.5,-0.5,360-0.5])
    plt.xticks(xticks+180-0.5,xticks,fontsize=fontsize)
    plt.yticks(xticks+180-0.5,xticks,fontsize=fontsize)
    plt.xlabel("eta",fontsize=fontsize)
    plt.ylabel("theta",fontsize=fontsize)
    #data_sum=1.*data_fine.sum()
    for i in range(data_fine.shape[0]):
        hpp+='{'
        for j in range(data_fine.shape[1]):
            #hpp+="%.5f,"%((data_fine[:i,:].sum()+ \
                        #data_fine[i,:(j+1)].sum())/data_sum)).lstrip('0')
            hpp+="%6d,"%(data_fine[:i,:].sum()+data_fine[i,:(j+1)].sum())
        hpp=hpp[:-1]+'},\n'
    hpp=hpp[:-2]+'''}
;
#endif
'''
    
    plt.subplot(2,2,2)
    step=30
    data_coarse=np.zeros((int(360/step),int(360/step)))
    for eta in range(-180,180,step):
        j=int((eta+180)/step)
        if j==data_coarse.shape[1]:
            j=0
        for theta in range(-180,180,step):
            i=int((theta+180)/step)
            if i==data_coarse.shape[0]:
                i=0
            data_coarse[i][j]=sum((eta_data>=  eta)*(  eta_data<  eta+step)*
                              (theta_data>=theta)*(theta_data<theta+step))
    vmax=data_coarse.max(axis=0).mean()
    plt.imshow(data_coarse,vmin=0,vmax=vmax,cmap="jet",interpolation='nearest')
    plt.axis([-0.5,360/step-0.5,-0.5,360/step-0.5])
    
    plt.xticks((xticks+180)/step-0.5,xticks,fontsize=fontsize)
    plt.yticks((xticks+180)/step-0.5,xticks,fontsize=fontsize)
    plt.xlabel("eta",fontsize=fontsize)
    plt.ylabel("theta",fontsize=fontsize)

    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    print("statNaTorsion2.png")
    plt.savefig("statNaTorsion2.png",dpi=300)
    fp=open("pseudo_torsion_prob.hpp",'w')
    fp.write(hpp)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()

    statNaTorsion2(sys.argv[1])
