#!/usr/bin/env python
docstring='''
statBPtorsion.py BPtorsion.raw.gz
    read torsion angle data, draw Torsion statistics
'''
import sys, os
import gzip

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import numpy as np

#import scipy.stats

def statBPtorsion(rawfile):
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

    BPtorsion_txt="#short\tALLmean\tALLstd\tWCmean\tWCstd\tGUmean\tGUstd\tfull\n"
    
    #### set parameters for matplotlib ####
    fontsize=7
    matplotlib.rcParams["legend.borderaxespad"]=0.1
    matplotlib.rcParams["legend.borderpad"]    =0.1
    matplotlib.rcParams["legend.columnspacing"]=1
    matplotlib.rcParams["legend.handlelength"] =0.5
    matplotlib.rcParams["legend.handletextpad"]=0.1 
    matplotlib.rcParams["legend.labelspacing"] =0.1
    
    xlabel_dict={
        "Pm"    : "P[i-1]-P[i]-P[j]-P[j-1]"        ,
        "Pp"    : "P[i+1]-P[i]-P[j]-P[j+1]"        ,
        "O5m"   : "O5'[i-1]-O5'[i]-O5'[j]-O5'[j-1]",
        "O5p"   : "O5'[i+1]-O5'[i]-O5'[j]-O5'[j+1]",
        "C5m"   : "C5'[i-1]-C5'[i]-C5'[j]-C5'[j-1]",
        "C5p"   : "C5'[i+1]-C5'[i]-C5'[j]-C5'[j+1]",
        "C4m"   : "C4'[i-1]-C4'[i]-C4'[j]-C4'[j-1]",
        "C4p"   : "C4'[i+1]-C4'[i]-C4'[j]-C4'[j+1]",
        "C3m"   : "C3'[i-1]-C3'[i]-C3'[j]-C3'[j-1]",
        "C3p"   : "C3'[i+1]-C3'[i]-C3'[j]-C3'[j+1]",
        "C2m"   : "C2'[i-1]-C2'[i]-C2'[j]-C2'[j-1]",
        "C2p"   : "C2'[i+1]-C2'[i]-C2'[j]-C2'[j+1]",
        "C1m"   : "C1'[i-1]-C1'[i]-C1'[j]-C1'[j-1]",
        "C1p"   : "C1'[i+1]-C1'[i]-C1'[j]-C1'[j+1]",
        "O4m"   : "O4'[i-1]-O4'[i]-O4'[j]-O4'[j-1]",
        "O4p"   : "O4'[i+1]-O4'[i]-O4'[j]-O4'[j+1]",
        "O3m"   : "O3'[i-1]-O3'[i]-O3'[j]-O3'[j-1]",
        "O3p"   : "O3'[i+1]-O3'[i]-O3'[j]-O3'[j+1]",
        "P44P"  : "P[i]-C4'[i]-C4'[j]-P[j]"        ,
        "C4114C": "C4'[i]-C1'[i]-C1'[j]-C4'[j]"    ,
        "PppP"  : "P[i+1]-P[i]-P[j]-P[j+1]"        ,

        "PP"    : "P[i]-P[j]"    ,
        "O5O5"  : "O5'[i]-O5'[j]",
        "C5C5"  : "C5'[i]-C5'[j]",
        "C4C4"  : "C4'[i]-C4'[j]",
        "C3C3"  : "C3'[i]-C3'[j]",
        "C2C2"  : "C2'[i]-C2'[j]",
        "C1C1"  : "C1'[i]-C1'[j]",
        "O4O4"  : "O4'[i]-O4'[j]",
        "O3O3"  : "O3'[i]-O3'[j]",
        "NN"    : "N[i]-N[j]"    ,

        "aPm"   : "<P[i-1]P[i],P[j+1]P[j]>"        ,
        "aPp"   : "<P[i+1]P[i],P[j-1]P[j]>"        ,
        "aO5m"  : "<O5'[i-1]O5'[i],O5'[j+1]O5'[j]>",
        "aO5p"  : "<O5'[i+1]O5'[i],O5'[j-1]O5'[j]>",
        "aC5m"  : "<C5'[i-1]C5'[i],C5'[j+1]C5'[j]>",
        "aC5p"  : "<C5'[i+1]C5'[i],C5'[j-1]C5'[j]>",
        "aC4m"  : "<C4'[i-1]C4'[i],C4'[j+1]C4'[j]>",
        "aC4p"  : "<C4'[i+1]C4'[i],C4'[j-1]C4'[j]>",
        "aC3m"  : "<C3'[i-1]C3'[i],C3'[j+1]C3'[j]>",
        "aC3p"  : "<C3'[i+1]C3'[i],C3'[j-1]C3'[j]>",
        "aC2m"  : "<C2'[i-1]C2'[i],C2'[j+1]C2'[j]>",
        "aC2p"  : "<C2'[i+1]C2'[i],C2'[j-1]C2'[j]>",
        "aC1m"  : "<C1'[i-1]C1'[i],C1'[j+1]C1'[j]>",
        "aC1p"  : "<C1'[i+1]C1'[i],C1'[j-1]C1'[j]>",
        "aO4m"  : "<O4'[i-1]O4'[i],O4'[j+1]O4'[j]>",
        "aO4p"  : "<O4'[i+1]O4'[i],O4'[j-1]O4'[j]>",
        "aO3m"  : "<O3'[i-1]O3'[i],O3'[j+1]O3'[j]>",
        "aO3p"  : "<O3'[i+1]O3'[i],O3'[j-1]O3'[j]>",
        "aPC"   : "<P[i]C4'[i],P[j]C4'[j]>"        ,
        "aCC"   : "<C4'[i]C1'[i],C4'[j]C1'[j]>"    ,
    }
    # WC          19-XIX    cWW  cW-W
    # Wobble      28-XXVIII cWW  cW-W
    name_list=["WC","Wobble"]

    #### dihedral ####
    plt.figure(figsize=(8,16))
    bins=np.arange(-180,181,1)
    xticks=range(-180,181,60)
    for a in range(20):
        print(angle_list[a])
        if a % 2 ==0:
            plt.subplot(10,5,5*a/2+1)
        else:
            plt.subplot(10,5,5*(a-1)/2+2)
        ymax=0
        angle_data=data[:,a]
        angle_data=angle_data[angle_data>=-180]
        mean=angle_data.mean()
        std =angle_data.std()
        angle_data2=np.array(angle_data)
        angle_data2[angle_data<0]+=360
        mean2=angle_data2.mean()
        std2 =angle_data2.std()
        if std2<std:
            std=std2
            mean=mean2
            if mean>180:
                mean-=360
        label="all(%.2f,%.2f)"%(mean,std)
        BPtorsion_txt+="%s\t%.3f\t%.3f\t"%(angle_list[a],mean,std)
        n=plt.hist(angle_data,bins=bins,
            edgecolor="none",facecolor="white",
            width=1, align="left", label=label)[0]
        for l,name in enumerate(name_list):
            angle_data=data[bptype_array==name,a]
            angle_data=angle_data[angle_data>=-180]
            mean=angle_data.mean()
            std =angle_data.std()
            angle_data2=np.array(angle_data)
            angle_data2[angle_data<0]+=360
            mean2=angle_data2.mean()
            std2 =angle_data2.std()
            if std2<std:
                std=std2
                mean=mean2
                if mean>180:
                    mean-=360
            label="%s(%.2f,%.2f)"%(name.replace("Wobble","g/u"), mean,std)
            BPtorsion_txt+="%.3f\t%.3f\t"%(mean,std)
            n=plt.hist(angle_data,bins=bins,
                edgecolor="none",facecolor=["lightgrey","grey"][l],
                width=1, align="left", label=label)[0]
            ymax=max([ymax,n.max()*1.5])
        BPtorsion_txt+=xlabel_dict[angle_list[a]]+'\n'
        plt.axis([min(bins),max(bins),0,ymax])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel_dict[angle_list[a]],fontsize=fontsize)
        plt.legend(loc="upper center",fontsize=fontsize)
        #if a in range(0,21,2):
            #plt.ylabel("Counts",fontsize=fontsize)
    
    #### distance ####
    bins=np.arange(8,24.00001,0.1)
    xticks=range(8,25,4)
    for a in range(20,30):
        plt.subplot(10,5,5*(a-20)+3)
        print(angle_list[a])
        ymax=0
        angle_data=data[:,a]
        angle_data=angle_data[angle_data>0]
        #mean, std = scipy.stats.norm.fit(angle_data, 
            #loc=angle_data.mean(),
            #scale=angle_data.std())
        mean=angle_data.mean()
        std=angle_data.std()
        label="all(%.3f,%.3f)"%(mean,std)
        BPtorsion_txt+="%s\t%.3f\t%.3f\t"%(angle_list[a],mean,std)
        n=plt.hist(angle_data,bins=bins,
            edgecolor="none",facecolor="white",
            width=1, align="left", label=label)[0]
        for l,name in enumerate(name_list):
            angle_data=data[bptype_array==name,a]
            #mean, std = scipy.stats.norm.fit(angle_data,
                #loc=angle_data.mean(),
                #scale=angle_data.std())
            mean=angle_data.mean()
            std=angle_data.std()
            label="%s(%.3f,%.3f)"%(name.replace("Wobble","g/u"),mean,std)
            BPtorsion_txt+="%.3f\t%.3f\t"%(mean,std)
            n=plt.hist(angle_data,bins=bins,
                edgecolor="none",facecolor=["lightgrey","grey"][l],
                width=1, align="left", label=label)[0]
            ymax=max([ymax,n.max()*1.5])
        BPtorsion_txt+=xlabel_dict[angle_list[a]]+'\n'
        plt.axis([min(bins),max(bins),0,ymax])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel_dict[angle_list[a]],fontsize=fontsize)
        plt.legend(loc="upper center",fontsize=fontsize)
    
    #### angle ####
    bins=np.arange(0,181,1)
    xticks=range(0,181,30)
    for a in range(30,50):
        print(angle_list[a])
        if a % 2 ==0:
            plt.subplot(10,5,5*(a-30)/2+4)
        else:
            plt.subplot(10,5,5*(a-31)/2+5)
        ymax=0
        angle_data=data[:,a]
        angle_data=angle_data[angle_data>=0]
        #fa,fb,mean,std=scipy.stats.truncnorm.fit(angle_data,
            #loc=angle_data.mean(),
            #scale=angle_data.std(),
            #fa=0,
            #fb=180,
            #)
        mean=angle_data.mean()
        std =angle_data.std()
        label="all(%.2f,%.2f)"%(mean,std)
        BPtorsion_txt+="%s\t%.3f\t%.3f\t"%(angle_list[a],mean,std)
        n=plt.hist(angle_data,bins=bins,
            edgecolor="none",facecolor="white",
            width=1, align="left", label=label)[0]
        for l,name in enumerate(name_list):
            angle_data=data[bptype_array==name,a]
            angle_data=angle_data[angle_data>=0]
            #fa,fb,mean,std=scipy.stats.truncnorm.fit(angle_data,
                #loc=angle_data.mean(),
                #scale=angle_data.std(),
                #fa=0,
                #fb=180,
                #)
            mean=angle_data.mean()
            std =angle_data.std()
            label="%s(%.2f,%.2f)"%(name.replace("Wobble","g/u"),mean,std)
            BPtorsion_txt+="%.3f\t%.3f\t"%(mean,std)
            n=plt.hist(angle_data,bins=bins,
                edgecolor="none",facecolor=["lightgrey","grey"][l],
                width=1, align="left", label=label)[0]
            ymax=max([ymax,n.max()*1.5])
        BPtorsion_txt+=xlabel_dict[angle_list[a]]+'\n'
        plt.axis([min(bins),max(bins),0,ymax])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel_dict[angle_list[a]],fontsize=fontsize)
        plt.legend(loc="upper center",fontsize=fontsize)

    fp=open("BPtorsion.txt",'w')
    fp.write(BPtorsion_txt)
    fp.close()
    plt.tight_layout(pad=0.05,h_pad=-0.5,w_pad=-0.5)
    plt.savefig("BPtorsion.png",dpi=300)
    return

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()

    statBPtorsion(sys.argv[1])
