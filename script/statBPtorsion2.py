#!/usr/bin/env python
docstring='''
statBPtorsion2.py BPtorsion2.raw.gz
    read torsion angle data, draw Torsion statistics
'''
import sys, os
import gzip

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import numpy as np

#import scipy.stats

def get_mean_std(angle_data,is_torsion=False):
    mu=angle_data.mean()
    sigma=angle_data.std()
    if is_torsion:
        angle_data2=np.array(angle_data)
        idx=angle_data<0
        angle_data2[idx]+=360
        sigma2=angle_data2.std()
        mu2=angle_data2.mean()
        if sigma2<sigma:
            mu=mu2
            sigma=sigma2
            if mu>180:
                mu-=360
        for it in range(10):
            angle_data2=np.array(angle_data)
            idx=angle_data<mu-180
            if sum(idx):
                angle_data2[idx]+=360
            idx=angle_data>mu+180
            if sum(idx):
                angle_data2[idx]-=360
            mu=angle_data2.mean()
            sigma=angle_data2.std()
            if mu>180:
                mu-=360
            elif mu<-180:
                mu+=360
    return mu,sigma

def statBPtorsion2(rawfile):
    data=[]
    bptype_list=[]
    angle_list=[]
    fp=gzip.open(rawfile)
    for line in fp.read().splitlines():
        if line.startswith("####"):
            continue
        if line.startswith("N c resi"):
            angle_list=line[32:].strip().split()
            continue
        name=line[20:32].strip()
        bptype_list.append(name)
        data.append([float(t) for t in line[32:].split()])
    fp.close()
    data=np.array(data)
    bptype_array=np.array(bptype_list,dtype=str)

    BPtorsion2_txt="#short\tALLmean\tALLstd\tWCmean\tWCstd\tGUmean\tGUstd\tfull\n"
    
    #### set parameters for matplotlib ####
    fontsize=8
    matplotlib.rcParams["legend.borderaxespad"]=0.1
    matplotlib.rcParams["legend.borderpad"]    =0.1
    matplotlib.rcParams["legend.columnspacing"]=1
    matplotlib.rcParams["legend.handlelength"] =0.5
    matplotlib.rcParams["legend.handletextpad"]=0.1 
    matplotlib.rcParams["legend.labelspacing"] =0.1
    
    xlabel_dict={
        "PCCP":"P[i]-C4'[i]-C4'[j]-P[j]",
        "CNNC":"C4'[i]-N[i]-N[j]-C4'[j]",
        "PNNP":"P[i]-N[i]-N[j]-P[j]"    ,

        "PP"  :"P[i]-P[j]"              ,
        "CC"  :"C4'[i]-C4'[j]"          ,
        "NN"  :"N[i]-N[j]"              ,

        "aPC" :"<P[i]C4'[i],P[j]C4'[j]>",
        "aCN" :"<C4'[i]N[i],C4'[j]N[j]>",
        "aPN" :"<P[i]N[i],P[j]N[j]>"    ,
    }
    # WC          19-XIX    cWW  cW-W
    # Wobble      28-XXVIII cWW  cW-W
    name_list=["WC","Wobble"]

    #### dihedral ####
    plt.figure(figsize=(7.87,7.87))
    bins=np.arange(-180,181,1)
    xticks=range(-180,181,60)
    for a in range(len(data[0])):
        print(angle_list[a])
        if a<3:
            width=1
            bins=np.arange(-180,181,width)
            xticks=range(-180,181,60)
        elif 3<=a and a<6:
            width=0.1
            bins=np.arange(8,24+width,width)
            xticks=range(8,25,4)
        else:
            width=1
            bins=np.arange(0,181,width)
            xticks=range(0,181,30)
        
        plt.subplot(3,3,a+1)
        ymax=0
        angle_data=data[:,a]
        if a<3 or a>6:
            angle_data=angle_data[angle_data>=-180]
        else:
            angle_data=angle_data[angle_data>0]
        mu,sigma=get_mean_std(angle_data,a<3)
        label="all(%.3f,%.3f)"%(mu,sigma)
        BPtorsion2_txt+="%s\t%.3f\t%.3f\t"%(angle_list[a],mu,sigma)
        n=plt.hist(angle_data,bins=bins,
            edgecolor="none",facecolor="white",
            width=1, align="left", label=label)[0]
        for l,name in enumerate(name_list):
            angle_data=data[bptype_array==name,a]
            if a<3 or a>6:
                angle_data=angle_data[angle_data>=-180]
            else:
                angle_data=angle_data[angle_data>0]
            mu,sigma=get_mean_std(angle_data,a<3)
            label="%s(%.3f,%.3f)"%(name.replace("Wobble","g/u"), mu,sigma)
            BPtorsion2_txt+="%.3f\t%.3f\t"%(mu,sigma)
            n=plt.hist(angle_data,bins=bins,
                edgecolor="none",facecolor=["lightgrey","grey"][l],
                width=1, align="left", label=label)[0]
            ymax=max([ymax,n.max()*1.3])
        BPtorsion2_txt+=xlabel_dict[angle_list[a]]+'\n'
        plt.axis([min(bins),max(bins),0,ymax])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel_dict[angle_list[a]],fontsize=fontsize)
        plt.legend(loc="upper center",fontsize=fontsize)
    
    fp=open("BPtorsion2.txt",'w')
    fp.write(BPtorsion2_txt)
    fp.close()
    plt.tight_layout(pad=0.1,h_pad=0,w_pad=0)
    plt.savefig("BPtorsion2.png",dpi=300)
    return

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()

    statBPtorsion2(sys.argv[1])
