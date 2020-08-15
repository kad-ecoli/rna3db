#!/usr/bin/python
fp=open("BPtorsion.txt",'r')
lines=fp.read().splitlines()
fp.close()

txt='''/* parameters for RNA base pairing */
#ifndef BPstat_HPP
#define BPstat_HPP 1
'''

for line in lines[1:]:
    line=line.split()
    short=line[0]
    mu=short+"_mu"
    mu="const float "+mu+' '*(9-len(mu))+"="+' '*(7-len(line[1]))+line[1]
    sd=short+"_sd"
    sd="const float "+sd+' '*(9-len(sd))+"="+' '*(6-len(line[2]))+line[2]
    txt+="%s; %s; // %s\n"%(mu,sd,line[-1])
txt+='''#endif
'''

fp=open("BPstat.hpp",'w')
fp.write(txt)
fp.close()
