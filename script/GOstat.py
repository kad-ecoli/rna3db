#!/usr/bin/env python
docstring='''
GOstat.py GO/goa/UNIPROT/goa_uniprot_all.gaf.gz obo/go/go-basic.obo stat.txt
    Read go annotations and go obo definition. Write statistics.
'''
import os, sys
import gzip
from obo2csv import parse_obo_file
import urllib2


guide_go_evidence_codes=dict(
    EXP="Inferred from Experiment (Experimental)",
    IDA="Inferred from Direct Assay (Experimental)",
    IPI="Inferred from Physical Interaction (Experimental)",
    IMP="Inferred from Mutant Phenotype (Experimental)",
    IGI="Inferred from Genetic Interaction (Experimental)",
    IEP="Inferred from Expression Pattern (Experimental)",

    HTP="Inferred from High Throughput Experiment (high throughput)",
    HDA="Inferred from High Throughput Direct Assay (high throughput)",
    HMP="Inferred from High Throughput Mutant Phenotype (high throughput)",
    HGI="Inferred from High Throughput Genetic Interaction (high throughput)",
    HEP="Inferred from High Throughput Expression Pattern (high throughput)",

    IBA="Inferred from Biological aspect of Ancestor (Phylogenetically-inferred)",
    IBD="Inferred from Biological aspect of Descendant (Phylogenetically-inferred)",
    IKR="Inferred from Key Residues (Phylogenetically-inferred)",
    IRD="Inferred from Rapid Divergence (Phylogenetically-inferred)",

    ISS="Inferred from Sequence or structural Similarity (Computational analysis)",
    ISO="Inferred from Sequence Orthology (Computational analysis)",
    ISA="Inferred from Sequence Alignment (Computational analysis)",
    ISM="Inferred from Sequence Model (Computational analysis)",
    IGC="Inferred from Genomic Context (Computational analysis)",
    RCA="Inferred from Reviewed Computational Analysis (Computational analysis)",

    TAS="Traceable Author Statement (Author statement)",
    NAS="Non-traceable Author Statement (Author statement)",

    IC="Inferred by Curator (Curator statement)",
    ND="No biological Data available (Curator statement)",

    IEA="Inferred from Electronic Annotation (Electronic annotation)",
)

go_aspect_dict=dict(
    F="molecular_function",
    P="biological_process",
    C="cellular_component",
)

topnum=100

def GOstat(gaf_file,obo_file,stat_file):
    obo_dict=parse_obo_file(obo_file)
    fp=gzip.open(gaf_file,'r')
    lines=[line.split('\t') for line in fp.read().splitlines() \
        if not line.startswith('!')]
    fp.close()

    DB_Object_Type_dict    =dict()
    EXP_DB_Object_Type_dict=dict()
    ann_dict        =dict(F=dict(),P=dict(),C=dict())
    EXP_ann_dict    =dict(F=dict(),P=dict(),C=dict())
    DB_Object_dict  =dict(F=dict(),P=dict(),C=dict())
    EXP_Object_dict =dict(F=dict(),P=dict(),C=dict())
    GO_ID_dict      =dict(F=dict(),P=dict(),C=dict())
    EXP_GO_ID_dict  =dict(F=dict(),P=dict(),C=dict())
    GO_REF_dict     =dict()
    evidence_dict   =dict()
    withfrom_dict   =dict()
    taxon_dict      =dict()
    assigned_by_dict=dict()
    for items in lines:
        if items[6]=="ND" or "NOT" in items[3]:
            continue # skip negative annotation
        DB_Object =items[1]
        GO_ID     =items[4]
        GO_REF    =items[5]
        evidence  =items[6]
        withfrom  =items[7]
        aspect    =items[8]
        DB_Object_Type=items[11]
        taxon     =items[12].split('|')[0][len("taxon:"):]
        assigned  =items[14]
        annotation=DB_Object+'\t'+GO_ID

        #if GO_ID in ["GO:0003735","GO:0030533","GO:0006412","GO:0005840"]:
            #continue
        if evidence in ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", 
                        "HTP", "HDA",        "HMP", "HGI", "HEP", "TAS"]:
            if not GO_ID in EXP_GO_ID_dict[aspect]:
                EXP_GO_ID_dict[aspect][GO_ID]=[]
            EXP_GO_ID_dict[aspect][GO_ID].append(DB_Object)
            if not DB_Object in EXP_Object_dict[aspect]:
                EXP_Object_dict[aspect][DB_Object]=0
            EXP_Object_dict[aspect][DB_Object]+=1
            if not annotation in EXP_ann_dict[aspect]:
                EXP_ann_dict[aspect][annotation]=0
            EXP_ann_dict[aspect][annotation]+=1
            if not DB_Object_Type in EXP_DB_Object_Type_dict:
                EXP_DB_Object_Type_dict[DB_Object_Type]=0
            EXP_DB_Object_Type_dict[DB_Object_Type]+=1

        if not DB_Object in DB_Object_dict[aspect]:
            DB_Object_dict[aspect][DB_Object]=0
        DB_Object_dict[aspect][DB_Object]+=1

        if not annotation in ann_dict[aspect]:
            ann_dict[aspect][annotation]=0
        ann_dict[aspect][annotation]+=1

        if not GO_ID in GO_ID_dict[aspect]:
            GO_ID_dict[aspect][GO_ID]=[]
        GO_ID_dict[aspect][GO_ID].append(DB_Object)

        if not GO_REF in GO_REF_dict:
            GO_REF_dict[GO_REF]=0
        GO_REF_dict[GO_REF]+=1

        if not evidence in evidence_dict:
            evidence_dict[evidence]=0
        evidence_dict[evidence]+=1

        if not withfrom in withfrom_dict:
            withfrom_dict[withfrom]=0
        withfrom_dict[withfrom]+=1
        
        if not DB_Object_Type in DB_Object_Type_dict:
            DB_Object_Type_dict[DB_Object_Type]=0
        DB_Object_Type_dict[DB_Object_Type]+=1

        if not taxon in taxon_dict:
            taxon_dict[taxon]=0
        taxon_dict[taxon]+=1

        if not assigned in assigned_by_dict:
            assigned_by_dict[assigned]=0
        assigned_by_dict[assigned]+=1

    stat_txt ="#### number of all child annotations with any code and ####\n"
    stat_txt+="#### number of experimental child annotations with     ####\n"
    stat_txt+="#### (EXP|IDA|IPI|IMP|IGI|IEP|HTP|HDA|HMP|HGI|HEP|TAS) ####\n"
    stat_txt+="#aspect\tcount_all\tcount_all_nonredundant\tcount_EXP\tcount_EXP_nonredundant\tname\n"
    data=[(sum(ann_dict[aspect].values()),len(ann_dict[aspect]),
           sum(EXP_ann_dict[aspect].values()),len(EXP_ann_dict[aspect]),
           aspect) for aspect in "FPC"]
    tot_nrnum=0
    tot_rnum=0
    tot_expnrnum=0
    tot_exprnum=0
    for rnum,nrnum,exprnum,expnrnum,aspect in data:
        tot_nrnum   +=nrnum
        tot_rnum    +=rnum
        tot_expnrnum+=expnrnum
        tot_exprnum +=exprnum
        stat_txt+="%s\t%d\t%d\t%d\t%d\t%s\n"%(aspect,nrnum,rnum,
            expnrnum,exprnum,go_aspect_dict[aspect])
    stat_txt+="A\t%d\t%d\t%d\t%d\tall_three_aspects\n"%(
        tot_rnum,tot_nrnum,tot_exprnum,tot_expnrnum)
    
    stat_txt+="#### all %d directly annotated GO terms ####\n"%(
        sum([len(GO_ID_dict[aspect]) for aspect in "FPC"]))
    for aspect in "FPC":
        stat_txt+="## top %d GO terms for %s among %d RNAs ##\n"%(
            topnum,aspect,len(DB_Object_dict[aspect]))
        stat_txt+="#GO term\tcount\tfrequency\tname\n"
        for count,GO_ID in sorted([(len(set(GO_ID_dict[aspect][GO_ID])),GO_ID
            ) for GO_ID in GO_ID_dict[aspect]],reverse=True)[:topnum]:
            name=obo_dict.short(GO_ID).strip().split(' ! ')[1]
            stat_txt+="%s\t%d\t%.7f\t%s\n"%(GO_ID,count,
                1.*count/len(DB_Object_dict[aspect]),name)

    stat_txt+="#### %d directly annotated experimental GO terms with  ####\n"%(
        sum([len(GO_ID_dict[aspect]) for aspect in "FPC"]))
    stat_txt+="#### (EXP|IDA|IPI|IMP|IGI|IEP|HTP|HDA|HMP|HGI|HEP|TAS) ####\n"
    for aspect in "FPC":
        stat_txt+="## %d experimental GO terms for %s among %d RNAs ##\n"%(
            len(EXP_GO_ID_dict[aspect]),aspect,len(EXP_Object_dict[aspect]))
        stat_txt+="#GO term\tcount\tfrequency\tcount_all\tfrequency_all\tname\n"
        for count,GO_ID in sorted([(len(set(EXP_GO_ID_dict[aspect][GO_ID])),GO_ID
            ) for GO_ID in EXP_GO_ID_dict[aspect]],reverse=True):
            name=obo_dict.short(GO_ID).strip().split(' ! ')[1]
            count_all=len(set(GO_ID_dict[aspect][GO_ID]))
            stat_txt+="%s\t%d\t%.7f\t%d\t%.7f\t%s\n"%(GO_ID,count,
                1.*count/len(EXP_Object_dict[aspect]),count_all,
                1.*count_all/len(DB_Object_dict[aspect]), name)

    stat_txt+="#### %d evidence codes ####\n"%len(evidence_dict)
    stat_txt+="#evidence\tcount\tfrequency\tname\n"
    data=sorted([(evidence_dict[evidence],evidence) for evidence in \
        evidence_dict],reverse=True)
    tot_rnum=sum(evidence_dict.values())
    for count,evidence in data:
        stat_txt+="%s\t%d\t%.7f\t%s\n"%(evidence,count,
            1.*count/tot_rnum,guide_go_evidence_codes[evidence])
    
    stat_txt+="#### top %d references ####\n"%topnum
    stat_txt+="#reference\tcount\tfrequency\tname\n"
    data=sorted([(GO_REF_dict[reference],reference) for reference in \
        GO_REF_dict],reverse=True)[:topnum]
    tot_rnum=sum(GO_REF_dict.values())
    for count,reference in data:
        name=reference
        if reference.startswith("GO_REF:"):
            for line in urllib2.urlopen(
                "https://raw.githubusercontent.com/geneontology/go-site/master/metadata/gorefs/goref-%s.md"%(reference.split(':')[1])):
                if line.startswith("## "):
                    name=line[len("## "):].strip()
                    break
        elif reference.startswith("PMID:"):
            for line in urllib2.urlopen(
                "https://pubmed.ncbi.nlm.nih.gov/%s/"%(reference.split(':')[1])):
                if "<title>" in line:
                    name=line.split('>')[1].split('<')[0]
                    break
        stat_txt+="%s\t%d\t%.7f\t%s\n"%(reference,count,1.*count/tot_rnum,name)

    stat_txt+="#### top %d with/from ####\n"%topnum
    stat_txt+="#with/from\tcount\tfrequency\tname\n"
    data=sorted([(withfrom_dict[withfrom],withfrom) for withfrom in \
        withfrom_dict],reverse=True)[:topnum]
    tot_rnum=sum(withfrom_dict.values())
    for count,withfrom in data:
        name=withfrom
        name_list=[]
        for name in withfrom.split('|'):
            if name.startswith("GO:"):
                name_list.append(obo_dict.short(withfrom).strip().split(' ! ')[1])
            elif name.startswith("Rfam:"):
                for line in urllib2.urlopen(
                    "http://rfam.xfam.org/family/%s/"%(name.split(':')[1])):
                    if "<title>" in line:
                        name_list.append(line.split('>')[1].split('<')[0].strip())
                        break
        name=' | '.join(name_list)
        stat_txt+="%s\t%d\t%.7f\t%s\n"%(withfrom,count,1.*count/tot_rnum,name)
    
    stat_txt+="#### top %d taxa ####\n"%topnum
    stat_txt+="#taxon\tcount\tfrequency\tname\n"
    data=sorted([(taxon_dict[taxon],taxon) for taxon in \
        taxon_dict],reverse=True)[:topnum]
    tot_rnum=sum(taxon_dict.values())
    for count,taxon in data:
        name=""
        for line in urllib2.urlopen(
            "https://www.uniprot.org/taxonomy/%s.rdf"%taxon):
            if "<scientificName>" in line:
                name=line.split('>')[1].split('<')[0]
                break
        stat_txt+="%s\t%d\t%.7f\t%s\n"%(taxon,count,1.*count/tot_rnum,name)
    
    stat_txt+="#### assigned by ####\n"
    stat_txt+="#assigned\tcount\tfrequency\n"
    data=sorted([(assigned_by_dict[assigned],assigned) for assigned in \
        assigned_by_dict],reverse=True)
    tot_rnum=sum(assigned_by_dict.values())
    for count,assigned in data:
        stat_txt+="%s\t%d\t%.7f\n"%(assigned,count,1.*count/tot_rnum)
    
    stat_txt+="#### DB_Object_Type ####\n"
    stat_txt+="#type\tcount_all\tfrequency_all\tcount_EXP\tfrequency_all\n"
    data=sorted([(DB_Object_Type_dict[DB_Object_Type],DB_Object_Type
        ) for DB_Object_Type in DB_Object_Type_dict],reverse=True)
    tot_rnum=sum(DB_Object_Type_dict.values())
    tot_exprnum=sum(EXP_DB_Object_Type_dict.values())
    for count,DB_Object_Type in data:
        count_EXP=0
        if DB_Object_Type in EXP_DB_Object_Type_dict:
            count_EXP=EXP_DB_Object_Type_dict[DB_Object_Type]
        stat_txt+="%s\t%d\t%.7f\t%d\t%.7f\n"%(DB_Object_Type,
            count,1.*count/tot_rnum,count_EXP,1.*count_EXP/tot_exprnum)
        
    fp=open(stat_file,'w')
    fp.write(stat_txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=4:
        sys.stderr.write(docstring)
        exit()
    GOstat(sys.argv[1],sys.argv[2],sys.argv[3])
