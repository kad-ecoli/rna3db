#!/usr/bin/env python
docstring='''
statLvsDSSR.py ../cull/all_c1.0_s1.0/list.atomic ../cull/pdb_atom.sort.4.0_c0.8_s0.8 ../cull/all_c1.0_s1.0/DSSR/
    Read target list ../cull/all_c1.0_s1.0/list.atomic
    fasta file ../cull/pdb_atom.sort.4.0_c0.8_s0.8, and
    DSSR output folder ../cull/all_c1.0_s1.0/DSSR/
    Draw the relation between length and number of base pairs.
'''

import sys, os
import textwrap
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import numpy as np

def fasta2len(inputfasta,target_set,dssrfolder):
    len_dict=dict()
    na_dict=dict(a=0,c=0,g=0,u=0)
    na_pair_dict={
        "a-a":0, "a-c":0, "a-g":0, "a-u":0, 
        "c-a":0, "c-c":0, "c-g":0, "c-u":0, 
        "g-a":0, "g-c":0, "g-g":0, "g-u":0, 
        "u-a":0, "u-c":0, "u-g":0, "u-u":0, 
    }
    target_list=[]
    fp=open(inputfasta,'r')
    for block in ('\n'+fp.read()).split('\n>'):
        lines=block.splitlines()
        if len(lines)<2:
            continue
        target=lines[0].split()[0]
        if not target in target_set:
            print("skip unlisted target %s"%target)
        filename=os.path.join(dssrfolder,target+".dssr")
        if not os.path.isfile(filename):
            print("Warning! %s missing"%filename)
            continue
        sequence=''.join(lines[1:])
        for na in na_dict:
            na_dict[na]+=sequence.count(na)
        for i in range(len(sequence)):
            for j in range(i+1,len(sequence)):
                na_pair_dict[sequence[i]+'-'+sequence[j]]+=1
        len_dict[target]=len(sequence)
        target_list.append(target)
    fp.close()
    return target_list,len_dict,na_dict,na_pair_dict

def dssr2num(target_list,dssrfolder):
    dssr2num_dict=dict()
    na_bp_dict={
        "a-a":0, "a-c":0, "a-g":0, "a-u":0, 
        "c-a":0, "c-c":0, "c-g":0, "c-u":0, 
        "g-a":0, "g-c":0, "g-g":0, "g-u":0, 
        "u-a":0, "u-c":0, "u-g":0, "u-u":0, 
    }
    for target in target_list:
        dssr2num_dict[target]=[0,0]
        filename=os.path.join(dssrfolder,target+".dssr")
        #if not os.path.isfile(filename):
            #print("Warning! %s missing"%filename)
            #continue
        fp=open(filename,'r')
        lines=fp.read().splitlines()
        fp.close()
        if len(lines)<3:
            print("Warning! %s empty"%filename)
            continue
        dssr2num_dict[target]=[int(lines[2].split()[2]),0]
        for line in lines[4:]:
            items=line.strip().split()
            name=items[4]
            if not name in ["WC","Wobble"]:
                continue
            key=items[3].lower().replace('+','-')
            na_bp_dict[key]+=1
            dssr2num_dict[target][1]+=1
    return dssr2num_dict,na_bp_dict

def FitLvsDSSR(target_list,len_dict,dssr2num_dict,Lcut):
    x=np.array([len_dict[target] for target in target_list \
        if len_dict[target]<=Lcut])
    y=np.array([dssr2num_dict[target][1] for target in target_list \
        if len_dict[target]<=Lcut])
    a1,a0=np.linalg.lstsq( np.vstack([x,np.ones(len(x))]).T, y, rcond=None)[0]
    a=np.linalg.lstsq( np.vstack([x,np.zeros(len(x))]).T, y, rcond=None)[0][0]
    return a1,a0,a,x,y

def DrawLvsDSSR(target_list,len_dict,dssr2num_dict,prefix="LvsDSSR"):
    fontsize=8
    #Lcut_list=[100,350,1500,max(len_dict.values())]
    Lcut_list=[75,150,300,600,1200,max(len_dict.values())]
    plt.figure(figsize=(7.87,5.25))
    txt="#Lcut\ta1\ta0\ta\n"
    a1_list=[]
    a0_list=[]
    a_list =[]
    for l,Lcut in enumerate(Lcut_list):
        a1,a0,a,x,y=FitLvsDSSR(target_list,len_dict,dssr2num_dict,Lcut)
        txt+="%d\t%.2f\t%.2f\t%.2f\n"%(Lcut,a1,a0,a)
        a1_list.append(a1)
        a0_list.append(a0)
        a_list.append(a)
        plt.subplot(2,3,l+1)
        plt.plot(x,y,'o',markerfacecolor="lightgrey",clip_on=False)
        plt.plot([0,Lcut],[a0,a1*Lcut+a0],'k-',
            label="y = %.2f * x %.2f"%(a1,a0))
        plt.xlabel("Chain length ("+r"$\leq$"+"%d nt)"%Lcut,fontsize=fontsize)
        plt.ylabel("Base pairs per chain",fontsize=fontsize)
        plt.axis([0,Lcut*1.01,0,max(y)*1.01])
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(loc="upper left",fontsize=fontsize,handlelength=0,
            handletextpad=0.1,borderaxespad=0.1)
    xticks=range(0,Lcut,1000)
    plt.xticks(xticks,xticks,fontsize=fontsize)
    txt+="#mean\t%.2f\t%.2f\t%.2f\n"%(
        sum(a1_list)/len(a1_list),
        sum(a0_list)/len(a0_list),
        sum(a_list)/len(a_list))
    txt+="##bp=\ta1*L+a0\tor\ta*L\n"
    fp=open("%s.txt"%prefix,'w')
    fp.write(txt)
    fp.close()
    plt.tight_layout(pad=0.1,w_pad=0.1,h_pad=0.1)
    for outfmt in ["png"]:
        plt.savefig("%s.%s"%(prefix,outfmt),dpi=300)
    return

def DrawNaDSSR(na_dict,na_pair_dict,na_bp_dict,prefix="NAinDSSR"):
    fontsize=8
    base_type_list=list("aucg")
   
    txt=''
    plt.figure(figsize=(7.87,6))
    ax=plt.subplot(2,2,1)
    ax.tick_params(axis='x',length=0)
    x=np.arange(len(base_type_list))
    y=[na_dict[base] for base in base_type_list]
    plt.bar(x,y,color="lightgrey",align="center")
    plt.xticks(x,base_type_list,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylabel("number of nucleotide in sequence",fontsize=fontsize)
    txt+='#### number of nucleotide in sequence ####\n'
    #total=1.*sum(na_dict.values())
    for b,base in enumerate(base_type_list):
        plt.text(b,na_dict[base],"%d"%(na_dict[base]
            ),fontsize=fontsize,ha="center",va="bottom")
        txt+="%s\t%d\n"%(base,na_dict[base])
    
    color_list=["white","lightgrey","grey","darkgrey"]
    ax=plt.subplot(2,2,2)
    ax.tick_params(axis='x',length=0)
    width=0.2
    plt.ylabel("Number of nucleotide combination in sequence",fontsize=fontsize)
    txt+="#### Number of nucleotide combination in sequence ####\n"
    for b,base2 in enumerate(base_type_list):
        y=[]
        x=[]
        for a,base1 in enumerate(base_type_list):
            x.append(a+width*(b-1.5))
            y.append(na_pair_dict[base1+'-'+base2])
            plt.text(a+width*(b-1.5)*1.1,y[-1],str(y[-1]),rotation=90,
                fontsize=fontsize,ha="center",va="bottom")
            txt+="%s-%s\t%d\n"%(base1,base2,y[-1])
        plt.bar(x,y,width=width,color=color_list[b],align="center",label="*-"+base2)
    plt.xticks(range(len(base_type_list)),[a+"-*" for a in base_type_list],fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.axis([-0.5,len(base_type_list)-0.5,0,max(na_pair_dict.values())*1.2])
    plt.legend(loc="upper left",fontsize=fontsize,ncol=len(base_type_list), columnspacing=0.5,
        handletextpad=0.1,borderpad=0.1,handlelength=1,labelspacing=0.1,borderaxespad=0.1)
    
    plt.subplot(2,2,3)
    ax.tick_params(axis='x',length=0)
    x=np.arange(len(base_type_list))
    y=[]
    y1=[]
    y2=[]
    width=0.3
    plt.ylabel("number of nucleotide in base pair",fontsize=fontsize)
    txt+="#### number of nucleotide in base pair ####\n"
    txt+="#base\tupsteam\tdownstream\ttotal\n"
    #total=1.*sum(na_bp_dict.values())
    for b,base in enumerate(base_type_list):
        y1.append(sum([na_bp_dict[key] for key in na_bp_dict if key[0]==base]))
        y2.append(sum([na_bp_dict[key] for key in na_bp_dict if key[-1]==base]))
        y.append(y1[-1]+y2[-1])
        txt+="%s\t%d\t%s\t%d\n"%(base,y1[-1],y2[-1],y[-1])
        plt.text(b-width,y[-1], str(y[-1]), fontsize=fontsize,ha="center",va="bottom",rotation=90)
        plt.text(b,      y1[-1],str(y1[-1]),fontsize=fontsize,ha="center",va="bottom",rotation=90)
        plt.text(b+width,y2[-1],str(y2[-1]),fontsize=fontsize,ha="center",va="bottom",rotation=90)
    plt.bar(x-width,y, width=width,color=color_list[0],align="center",label="Either base of the pair")
    plt.bar(x      ,y1,width=width,color=color_list[1],align="center",label="Upstream base in pair")
    plt.bar(x+width,y2,width=width,color=color_list[2],align="center",label="Downstream base in pair")
    plt.legend(loc="upper left",fontsize=fontsize, columnspacing=0.5,
        handletextpad=0.1,borderpad=0.1,handlelength=1,labelspacing=0.1,borderaxespad=0.1)
    plt.xticks(x,base_type_list,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.axis((-0.5,len(base_type_list)-0.5,1,max(y)*1.2))

    ax=plt.subplot(2,2,4)
    ax.tick_params(axis='x',length=0)
    width=0.9
    plt.ylabel("Number of nucleotide combination in base pair",fontsize=fontsize)
    txt+="#### Number of nucleotide combination in base pair ####\n"
    y=[]
    xticks=[]
    for b,base2 in enumerate(base_type_list):
        for a,base1 in enumerate(base_type_list):
            key=base1+'-'+base2
            if na_bp_dict[key]==0:
                continue
            xticks.append(key)
            y.append(na_bp_dict[base1+'-'+base2])
            plt.text(len(y)-1,y[-1],str(y[-1]),fontsize=fontsize,ha="center",va="bottom")
    plt.bar(range(len(y)),y,width=width,color="lightgrey",align="center")
    plt.xticks(range(len(xticks)),xticks,fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.axis([-0.5,len(y)-0.5,0,max(y)*1.2])
    #plt.legend(loc="upper left",fontsize=fontsize,ncol=len(base_type_list),columnspacing=0.5, 
        #handletextpad=0.1,borderpad=0.1,handlelength=1,labelspacing=0.1,borderaxespad=0.1)

    plt.tight_layout(pad=0.1,w_pad=0.1,h_pad=0.1)
    for outfmt in ["png"]:
        plt.savefig("%s.%s"%(prefix,outfmt),dpi=300)
    fp=open(prefix+".txt",'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)<=3:
        sys.stderr.write(docstring)
        exit()

    listfile  =sys.argv[1]
    inputfasta=sys.argv[2]
    dssrfolder=sys.argv[3]
    
    fp=open(listfile,'r')
    target_set=set([line.split()[0] for line in fp.read().splitlines()])
    fp.close()

    target_list,len_dict,na_dict,na_pair_dict=fasta2len(
        inputfasta,target_set,dssrfolder)
    print("%d targets"%len(target_list))
    print("%d nucleotides"%sum(len_dict.values()))
    
    dssr2num_dict,na_bp_dict=dssr2num(target_list,dssrfolder)
    
    print("Drawing LvCanonicalPair.png")
    DrawLvsDSSR(target_list,len_dict,dssr2num_dict,prefix="LvsCanonicalPair")
    print("Drawing NAinCanonicalPair.png")
    DrawNaDSSR(na_dict,na_pair_dict,na_bp_dict,prefix="NAinCanonicalPair")
