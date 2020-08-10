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

def statNaTorsion(rawfile):
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
    xticks=range(-180,181,60)
    
    xlabel_dict={
        "v0":"v0 (C4'-O4'-C1'-C2')",
        "v1":"v1 (O4'-C1'-C2'-C3')",
        "v2":"v2 (C1'-C2'-C3'-C4')",
        "v3":"v3 (C2'-C3'-C4'-O4')",
        "v4":"v4 (C3'-C4'-O4'-C1')",
        "v5":"v5 (C5'-C4'-C3'-C2')",
        "v6":"v6 (C4'-C3'-C2'-O2')",

        "PmP"    :"P[-1]-P"    ,
        "O5'mO5'":"O5'[-1]-O5'",
        "C5'mC5'":"C5'[-1]-C5'",
        "C4'mC4'":"C4'[-1]-C4'",
        "C3'mC3'":"C3'[-1]-C3'",
        "O3'mO3'":"O3'[-1]-O3'",
        "C2'mC2'":"C2'[-1]-C2'",
        "C1'mC1'":"C1'[-1]-C1'",
        "O4'mO4'":"O4'[-1]-O4'",
        "PC4'"   :"P-C4'"      ,
        "C4'Pp"  :"C4'-P[+1]"  ,
        "PC1'"   :"P-C1'"      ,
        "C1'Pp"  :"C1'-P[+1]"  ,
        "O3'mO5'":"O3'[-1]-O5'",
        "PC5'"   :"P-C5'"      ,
        "C5'C3'" :"C5'-C3'"    ,
        "C4'O3'" :"C4'-O3'"    ,
        "C3'Pp"  :"C3'-P[+1]"  ,
        "C4'C2'" :"C4'-C2'"    ,
        "O4'C3'" :"O4'-C3'"    ,
        "C4'C1'" :"C4'-C1'"    ,
        "C3'C1'" :"C3'-C1'"    ,
        "O4'N"   :"O4'-N"      ,

        "O3'mP"  :"O3'[-1]-P",
        "PO5'"   :"P-O5'"    ,
        "O5'C5'" :"O5'-C5'"  ,
        "C5'C4'" :"C5'-C4'"  ,
        "C4'C3'" :"C4'-C3'"  ,
        "C3'O3'" :"C3'-O3'"  ,
        "C3'C2'" :"C3'-C2'"  ,
        "C2'C1'" :"C2'-C1'"  ,
        "C4'O4'" :"C4'-O4'"  ,
        "O4'C1'" :"O4'-C1'"  ,
        "C2'O2'" :"C2'-O2'"  ,
        "C1'N"   :"C1'-N"    ,
        "NC"     :"N-C"      ,

        "C4'mPC4'" : "C4'[-1]-P-C4'",
        "PC4'Pp"   : "P-C4'-P[+1]"  ,
        "C1'mPC1'" : "C1'[-1]-P-C1'",
        "PC1'Pp"   : "P-C1'-P[+1]"  ,
        "PC4'C1'"  : "P-C4'-C1'"    ,
        
        "O3'mPO5'" : "O3'[-1]-P-O5'",
        "PO5'C5'"  : "P-O5'-C5'"    ,
        "O5'C5'C4'": "O5'-C5'-C4'"  ,
        "C5'C4'C3'": "C5'-C4'-C3'"  ,
        "C4'C3'O3'": "C4'-C3'-O3'"  ,
        "C3'O3'Pp" : "C3'-O3'-P[+1]",
        "C4'C3'C2'": "C4'-C3'-C2'"  ,
        "O3'C3'C2'": "O3'-C3'-C2'"  ,
        "C3'C2'C1'": "C3'-C2'-C1'"  ,
        "C3'C4'O4'": "C3'-C4'-O4'"  ,
        "C5'C4'O4'": "C5'-C4'-O4'"  ,
        "C4'O4'C1'": "C4'-O4'-C1'"  ,
        "O4'C1'C2'": "O4'-C1'-C2'"  ,
        "O2'C2'C3'": "O2'-C2'-C3'"  ,
        "O4'C1'N"  : "O4'-C1'-N"    ,
        "C1'NC"    : "C1'-N-C"      ,
    }

    #### pucker ####
    plt.figure(figsize=(7.87,2.62))
    plt.subplot(1,3,1)
    angle_data=data[:,2]
    angle_data=angle_data[angle_data>=-180]
    n1=plt.hist(angle_data[angle_data>0],bins=np.arange(-180,181,1),
        edgecolor="none",facecolor="blue",width=1,align="left",
        label="C3' endo (%d%%)"%(.5+100.*sum(angle_data>0)/len(angle_data)))[0]
    n2=plt.hist(angle_data[angle_data<0],bins=np.arange(-180,181,1),
        edgecolor="none",facecolor="red",width=1,align="left",
        label="C2' endo (%d%%)"%(.5+100.*sum(angle_data<0)/len(angle_data)))[0]
    plt.axis([-180,180,0,max((n1.max(),n2.max()))*1.1])
    plt.legend(loc="upper left",fontsize=fontsize)
    plt.xticks(xticks,xticks,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylabel("Counts",fontsize=fontsize)
    plt.xlabel("v2 (C1'-C2'-C3'-C4')",fontsize=fontsize)
    for a,na in enumerate("aucg"):
        plt.subplot(2,3,a+2+(a>=2))
        angle_data=data[N_array==na,2]
        angle_data=angle_data[angle_data>=-180]
        n1=plt.hist(angle_data[angle_data>0],bins=np.arange(-180,181,1),
            edgecolor="none",facecolor="blue",width=1,align="left",
            label="C3' endo (%d%%)"%(.5+100.*sum(angle_data>0)/len(angle_data)))[0]
        n2=plt.hist(angle_data[angle_data<0],bins=np.arange(-180,181,1),
            edgecolor="none",facecolor="red",width=1,align="left",
            label="C2' endo (%d%%)"%(.5+100.*sum(angle_data<0)/len(angle_data)))[0]
        plt.axis([-180,180,0,max((n1.max(),n2.max()))*1.1])
        plt.legend(loc="upper left",fontsize=fontsize)
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel("v2 for "+na,fontsize=fontsize)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("pucker.png",dpi=300)

    #### pseudo torsion ####
    plt.figure(figsize=(7.87,7.87))
    for row in range(2):
        for col in range(2):
            for l,label in enumerate(["C3' endo","C2' endo"]):
                plt.subplot(4,4,col+3+8*row+4*l)
                angle_data=data[pucker_array==1-l*2,col+7+row*2]
                n=plt.hist(angle_data,bins=np.arange(-180,181,1),
                    edgecolor="none",facecolor=["blue","red"][l],
                    width=1, align="left", label=label)[0]
                plt.axis([-180,180,0,n.max()*1.1])
                plt.xticks(xticks,xticks,fontsize=fontsize)
                plt.yticks(fontsize=fontsize)
                plt.xlabel(angle_list[col+7+row*2],fontsize=fontsize)
                if col==0:
                    plt.ylabel("Counts (%s)"%label,fontsize=fontsize)

        plt.subplot(2,2,1+2*row)
        x_data=data[pucker_array==-1,7+row*2]
        y_data=data[pucker_array==-1,8+row*2]
        plt.plot(x_data,y_data,'.',color="red",markersize=0.05)
        x_data=data[pucker_array==1,7+row*2]
        y_data=data[pucker_array==1,8+row*2]
        plt.plot(x_data,y_data,'.',color="blue",markersize=0.05)
        plt.axis([-180,180,-180,180])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(xticks,xticks,fontsize=fontsize)
        plt.xlabel(angle_list[7+row*2],fontsize=fontsize)
        plt.ylabel(angle_list[8+row*2],fontsize=fontsize)
        
    plt.tight_layout(pad=0.1,h_pad=-0.3,w_pad=-1.5)
    plt.savefig("pseudotorsion.png",dpi=300)

    #### sugar torsion ####
    plt.figure(figsize=(7.87,7.87))
    for a in range(0,7):
        angle=angle_list[a]
        for l,label in enumerate(["C3' endo","C2' endo"]):
            plt.subplot(3,3,a+1)
            angle_data=data[pucker_array==1-l*2,a]
            n=plt.hist(angle_data,bins=np.arange(-180,181,1),
                edgecolor="none",facecolor=["blue","red"][l],
                width=1, align="left", label=label)[0]
            if l==0:
                plt.axis([-180,180,0,n.max()*1.1])
                plt.xticks(xticks,xticks,fontsize=fontsize)
                plt.yticks(fontsize=fontsize)
                plt.xlabel(xlabel_dict[angle_list[a]],fontsize=fontsize)
            if a in [0,3,6]:
                plt.ylabel("Counts",fontsize=fontsize)
            plt.legend(loc="upper right",fontsize=fontsize)

    
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("sugartorsion.png",dpi=300)

    #### backbone torsion ####
    plt.figure(figsize=(7.87,7.87))
    for a in range(11,18):
        angle=angle_list[a]
        for l,label in enumerate(["C3' endo","C2' endo"]):
            if a<=14:
                plt.subplot(4,4,(a-10)+l*4)
            else:
                plt.subplot(4,4,(a-6)+l*4)
            angle_data=data[pucker_array==1-l*2,a]
            n=plt.hist(angle_data,bins=np.arange(-180,181,1),
                edgecolor="none",facecolor=["blue","red"][l],
                width=1, align="left", label=label)[0]
            plt.axis([-180,180,0,n.max()*1.1])
            plt.xticks(xticks,xticks,fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.xlabel(angle_list[a],fontsize=fontsize)
            if a==17:
                plt.xlabel(angle_list[a]+" (C1'-P-C4'-P[+1])",fontsize=fontsize)
            if a in [11,15]:
                plt.ylabel("Counts (%s)"%label,fontsize=fontsize)
        
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("backbonetorsion.png",dpi=300)
    
    #### sidechain torsion ####
    plt.figure(figsize=(7.87,7.87))
    a=18
    angle=angle_list[a]
    for l,label in enumerate(["C3' endo","C2' endo"]):
        for col in range(3):
            plt.subplot(4,4,1+col+4*l)
            if col==0:
                angle_data=data[pucker_array==1-l*2,a]
                xlabel="All base "+angle_list[a]
                plt.ylabel("Counts (%s)"%label,fontsize=fontsize)
            elif col==1:
                angle_data=data[(pucker_array==1-l*2)*(
                    (N_array=='a')+(N_array=='g')),a]
                xlabel="Base a/g "+angle_list[a]
            elif col==2:
                angle_data=data[(pucker_array==1-l*2)*(
                    (N_array=='c')+(N_array=='u')),a]
                xlabel="Base c/u "+angle_list[a]
            n=plt.hist(angle_data,bins=np.arange(-180,181,1),
                edgecolor="none",facecolor=["blue","red"][l],
                width=1, align="left", label=label)[0]
            plt.axis([-180,180,0,n.max()*1.1])
            plt.xticks(xticks,xticks,fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.xlabel(xlabel,fontsize=fontsize)

        for col,base in enumerate("agcu"):
            plt.subplot(4,4,9+col+4*l)
            angle_data=data[(pucker_array==1-l*2)*(N_array==base),a]
            n=plt.hist(angle_data,bins=np.arange(-180,181,1),
                edgecolor="none",facecolor=["blue","red"][l],
                width=1, align="left", label=label)[0]
            plt.axis([-180,180,0,n.max()*1.1])
            plt.xticks(xticks,xticks,fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.xlabel("Base %s %s"%(base,angle_list[a]),fontsize=fontsize)
            if col==0:
                plt.ylabel("Counts (%s)"%label,fontsize=fontsize)
        
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("sidechaintorsion.png",dpi=300)
    
    #### pseudobond length ####
    plt.figure(figsize=(7.87,7.87))
    bins=np.arange(0,12,0.1)
    xticks=range(0,13,2)
    width=0.1
    for a in range(19,32):
        xlabel=xlabel_dict[angle_list[a]]
        plt.subplot(4,4,a-18)
        bond_data=data[:,a]
        bond_data=bond_data[bond_data>0]
        bond_data=bond_data[bond_data<12]
        label=r"$\mu$"+"=%.3f\n"%bond_data.mean()+ \
            r"$\sigma$"+"=%.3f"%bond_data.std()
        n=plt.hist(bond_data,bins=bins,
            edgecolor="none",facecolor="grey",
            width=width, align="left", label=label)[0]
        plt.axis([0,max(xticks),0,n.max()*1.1])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.legend(loc="upper right",handlelength=0,fontsize=fontsize)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("pseudobond.png",dpi=300)

    plt.figure(figsize=(7.87,5.90))
    bins=np.arange(0,4,0.001)
    xticks=range(0,5,1)
    width=0.001
    for a in range(32,42):
        xlabel=xlabel_dict[angle_list[a]]
        plt.subplot(3,4,a-31)
        bond_data=data[:,a]
        bond_data=bond_data[bond_data>0]
        bond_data=bond_data[bond_data<4]
        label=r"$\mu$"+"=%.3f\n"%bond_data.mean()+ \
            r"$\sigma$"+"=%.3f"%bond_data.std()
        n=plt.hist(bond_data,bins=bins,
            edgecolor="none",facecolor="grey",
            width=width, align="left", label=label)[0]
        plt.axis([0,max(xticks),0,n.max()*1.1])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.legend(loc="upper left",handlelength=0,fontsize=fontsize)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("pseudobond2.png",dpi=300)
    
    #### covalent bond length ####
    plt.figure(figsize=(7.87,7.87))
    for a in range(42,55):
        xlabel=xlabel_dict[angle_list[a]]
        plt.subplot(4,4,a-41)
        bond_data=data[:,a]
        bond_data=bond_data[bond_data>0]
        bond_data=bond_data[bond_data<4]
        label=r"$\mu$"+"=%.3f\n"%bond_data.mean()+ \
            r"$\sigma$"+"=%.3f"%bond_data.std()
        n=plt.hist(bond_data,bins=bins,
            edgecolor="none",facecolor="grey",
            width=width, align="left", label=label)[0]
        plt.axis([0,max(xticks),0,n.max()*1.1])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.legend(loc="upper right",handlelength=0,fontsize=fontsize)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("bond.png",dpi=300)

    #### pseudo angle ####
    plt.figure(figsize=(7.87,5.25))
    bins=np.arange(0,181,1)
    xticks=range(0,181,30)
    width=1
    for a in range(55,60):
        xlabel=xlabel_dict[angle_list[a]]
        plt.subplot(2,3,a-54)
        angle_data=data[:,a]
        angle_data=angle_data[angle_data>=0]
        label=r"$\mu$"+"=%.2f\n"%angle_data.mean()+ \
            r"$\sigma$"+"=%.2f"%angle_data.std()
        n=plt.hist(angle_data,bins=bins,
            edgecolor="none",facecolor="grey",
            width=width, align="left", label=label)[0]
        plt.axis([0,max(xticks),0,n.max()*1.1])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.legend(loc="upper left",handlelength=0,fontsize=fontsize)
    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("pseudoangle.png",dpi=300)

    #### angle ####
    plt.figure(figsize=(7.87,7.87))
    bins=np.arange(0,181,0.01)
    width=0.01
    for a in range(60,75):
        xlabel=xlabel_dict[angle_list[a]]
        plt.subplot(4,4,a-59)
        angle_data=data[:,a]
        angle_data=angle_data[angle_data>=0]
        label=r"$\mu$"+"=%.2f\n"%angle_data.mean()+ \
            r"$\sigma$"+"=%.2f"%angle_data.std()
        n=plt.hist(angle_data,bins=bins,
            edgecolor="none",facecolor="grey",
            width=width, align="left", label=label)[0]
        plt.axis([0,max(xticks),0,n.max()*1.1])
        plt.xticks(xticks,xticks,fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.legend(loc="upper left",handlelength=0,fontsize=fontsize)

    a=75
    xlabel=xlabel_dict[angle_list[a]]
    plt.subplot(4,4,a-59)
    angle_data=data[(N_array=='a')+(N_array=='g'),a]
    angle_data=angle_data[angle_data>=0]
    label="a/g:"+r"$\mu$"+"=%.2f;"%angle_data.mean()+ \
            r"$\sigma$"+"=%.2f"%angle_data.std()
    n1=plt.hist(angle_data,bins=bins,
        edgecolor="none",facecolor="blue",
        width=width, align="left", label=label)[0]
    angle_data=data[(N_array=='c')+(N_array=='u'),a]
    angle_data=angle_data[angle_data>=0]
    label="c/u:"+r"$\mu$"+"=%.2f;"%angle_data.mean()+ \
            r"$\sigma$"+"=%.2f"%angle_data.std()
    n2=plt.hist(angle_data,bins=bins,
        edgecolor="none",facecolor="red",
        width=width, align="left", label=label)[0]
    plt.axis([0,max(xticks),0,max([n1.max(),n2.max()])*1.2])
    plt.xticks(xticks,xticks,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    plt.legend(loc="upper left",fontsize=fontsize)

    plt.tight_layout(pad=0.1,h_pad=0.1,w_pad=0.1)
    plt.savefig("bondangle.png",dpi=300)
    return

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()

    statNaTorsion(sys.argv[1])
