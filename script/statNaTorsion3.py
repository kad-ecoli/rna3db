#!/usr/bin/env python
docstring='''
statNaTorsion3.py NaTorsion2.raw.gz
    read torsion angle data, draw Torsion statistics
'''
import sys, os
import gzip

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import numpy as np

def statNaTorsion(rawfile):
    data=[]
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
        data.append([float(t) for t in line[9:].split()])
    fp.close()
    data=np.array(data)
    N_array=np.array(N_list,dtype=str)
    
    #### set parameters for matplotlib ####
    fontsize=8
    matplotlib.rcParams["legend.borderaxespad"]=0.1
    matplotlib.rcParams["legend.borderpad"]    =0.1
    matplotlib.rcParams["legend.columnspacing"]=1
    matplotlib.rcParams["legend.handlelength"] =1
    matplotlib.rcParams["legend.handletextpad"]=0.1 
    matplotlib.rcParams["legend.labelspacing"] =0.1
    xticks=range(-180,181,60)
    
    xlabel_dict={
        "eta"    :"eta (C4'[-1]-P-C4'-P[+1])"  ,
        "theta"  :"theta (P-C4'-P[+1]-C4'[+1])",
        "PPC4C1" :"P[+1]-P-C4'-C1'"            ,
        "PPC4N"  :"P[+1]-P-C4'-N"              ,
        "C4PC4C1":"C4'[-1]-P-C4'-C1'"          ,
        "C4PC4N" :"C4'[-1]-P-C4'-N"            ,

        "PC4'"   :"P-C4'"      ,
        "C4'Pp"  :"C4'-P[+1]"  ,
        "C4'C1'" :"C4'-C1'"    ,
        "C4'N"   :"C4'-N"      ,
        "PPp"    :"P-P[+1]"    ,
        "C4'mC4'":"C4'[-1]-C4'",

        "C4'mPC4'": "C4'[-1]-P-C4'",
        "PC4'Pp"  : "P-C4'-P[+1]"  ,
        "PC4'C1'" : "P-C4'-C1'"    ,
        "PC4'N"   : "P-C4'-N"      ,
    }

    #### pseudobond length ####
    plt.figure(figsize=(7.87,7.87))

    for a in range(len(data[0])):
        if a<6:
            width=1
            bins=np.arange(-180,180+width,width)
            xticks=range(-180,181,60)
        elif 6<=a and a<10:
            width=0.1
            bins=np.arange(2,6+width,width)
            xticks=range(2,7,2)
        elif 10<=a and a<12:
            width=0.1
            bins=np.arange(2,8+width,width)
            xticks=range(2,9,2)
        elif 12<=a:
            width=1
            bins=np.arange(0,180+width,width)
            xticks=range(0,181,30)

        xlabel=xlabel_dict[angle_list[a]]
        plt.subplot(4,4,a+1)
        
        bond_data=data[:,a]
        mu=bond_data.mean()
        sigma=bond_data.std()
        if a<6:
            bond_data2=np.array(bond_data)
            idx=bond_data<0
            bond_data2[idx]=bond_data[idx]+360
            sigma2=bond_data2.std()
            mu2=bond_data2.mean()
            if sigma2<sigma:
                mu=mu2
                sigma=sigma2
                if mu>180:
                    mu-=360
            for it in range(10):
                bond_data2=np.array(bond_data)
                idx=bond_data<mu-180
                if sum(idx):
                    bond_data2[idx]=bond_data[idx]+360
                idx=bond_data>mu+180
                if sum(idx):
                    bond_data2[idx]=bond_data[idx]-360
                mu=bond_data2.mean()
                sigma=bond_data2.std()
                if mu>180:
                    mu-=360
                elif mu<-180:
                    mu+=360
        label=r"$\mu$"+"=%.3f\n"%mu+ \
            r"$\sigma$"+"=%.3f"%sigma
        n=plt.hist(bond_data,bins=bins,
            edgecolor="none",facecolor="grey",
            width=width, align="left", label=label)[0]
        ymax=n.max()
        loc="upper right"
        if a in [0,4,5,10,11,12,13,14,15]:
            loc="upper left"
        plt.legend(loc=loc,handlelength=0,fontsize=fontsize)
        plt.axis([min(xticks),max(xticks),0,ymax*1.1])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("NaTorsion2.png",dpi=300)
    return

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()

    statNaTorsion(sys.argv[1])
