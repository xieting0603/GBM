import sys
import os
import os.path
import numpy as np
import scipy.stats as stats
import random
genes={}
e={}
EXTEND=int(sys.argv[1])#######the distance flanking loop anchor to be consider

kind=sys.argv[3]########Protein coding gene or all genes
rootdir='.'
rsemdir='/scratch/users/txie/RNA-seq/3rsem_out/'###the place to put the rsem RNA results############
cutoff=sys.argv[2]
f=open("./neoloop95_"+kind+"_ext"+str(EXTEND//1000)+'kb'+'.'+str(cutoff)+'.tpm','w')
tpm_neoloops=[]
tpm_loops=[]
tpm_nonloops=[]
tpm_loopsLa=[]
tpm_neoloopsLa=[]
tpm_nonloopsLa=[]
locs_loop={}
locs_neoloop={}
locs_nonloop={}
sw={}
def get_filenames(file_dir):
    L=[]
    list = os.listdir(file_dir) 
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]) and not 'R' in list[i]:
            L.append(os.path.join(list[i]))
    return L
#####Function is for getting the list of loop and neoloop on reassemble SVs########
def get_location(filename,label):
    L_neoloop={}
    L_loop={}
    L_neoloop[label]=[]
    L_loop[label]=[]
    for i in open(rootdir+'/'+filename):
        t=i.strip().split()
        L_loop[label].append((t[0],int(t[1]),int(t[2])))
        L_loop[label].append((t[3],int(t[4]),int(t[5])))
        if '000,1' in t[-1]:
            L_neoloop[label].append((t[0],int(t[1]),int(t[2]),label))
            L_neoloop[label].append((t[3],int(t[4]),int(t[5]),label))
    L_loop[label]=set(L_loop[label])
    L_neoloop[label]=set(L_neoloop[label])
    return ((list(L_neoloop[label]),list(L_loop[label])))
for filename in get_filenames(rsemdir):
    fil=filename.split('.')[0]
    if  filename.endswith('.genes.results') and filename.startswith("GSC") :
        for i in open(rsemdir+filename):
            t=i.strip().split()
            if not t[0]=='gene_id':
                e[(fil,t[0].split('_')[0])]=np.log2(float(t[-2])+1)
#####Function is for getting the list of nonloops on reassemble SVs########
def get_SVlocation(filename,label):
    L_nonloop={}
    L_nonloop[label]=[]
    for i in open(rootdir+'/'+filename):
        t=i.strip().split()
        if t[0].startswith('C'):
            col1=t[1].split(',')
            col2=t[2].split(',')
            col3=t[3].split(',')
            start1=min(col1[2],col2[1])
            end1=max(col1[2],col2[1])
            start2=min(col1[5],col3[1])
            end2=max(col1[5],col3[1])
            L_nonloop[label].append(('chr'+col1[1],int(start1),int(end1)))
            L_nonloop[label].append(('chr'+col1[4],int(start2),int(end2)))
    return list(set(L_nonloop[label]))

for filename in get_filenames(rootdir):
    fil=filename.split('_')[0]
    if  filename.endswith('.assemblies.txt') and filename.startswith("GSC"):
        locs_nonloop[fil]=get_SVlocation(filename,fil)

if kind=='pc':
    for i in open(f"/scratch/users/txie/RNA-seq/0ref/v29.gff"):
        t=i.strip().split('\t')
        if t[-1]=='+':
            loc=int(t[3])
        if t[-1]=='-':
            loc=int(t[4])
        if not t[2] in genes.keys():
            genes[t[2]]=[]
        genes[t[2]].append((loc,t[0]))
        sw[t[0]]=t[1]
if kind=='all':
    for i in open(f"/scratch/users/txie/RNA-seq/0ref/v29_all.gff"):
        t=i.strip().split('\t')
        if t[-1]=='+':
            loc=int(t[3])
        if t[-1]=='-':
            loc=int(t[4])
        if not t[2] in genes.keys():
            genes[t[2]]=[]
        genes[t[2]].append((loc,t[0]))
        sw[t[0]]=t[1]
f.write("geneID"+'\t'+"sample"+'\t'+"type"+'\t'+"chr1"+"\t"+'loc1'+"\t"+"loc2"+"\t"+"log2(tpm+1)"+'\n')

for filename in get_filenames(rootdir):
    fil=filename.split('.')[0]
    if  filename.endswith('95.Mneoloop') and filename.startswith("GSC") and not 'R' in filename:
        locs_neoloop[fil]=get_location(filename,fil)[0]
        locs_loop[fil]=get_location(filename,fil)[1]

#################get the rna levels of three list neoloop associated genes, loop associated genes and nonloop associated genes and all the genes are on the SVs#########
for label in locs_loop.keys():
    locs=locs_neoloop[label]
    for loc in locs:
        if loc[2]-loc[1]>10000:
            extend=0
        else:
            extend=EXTEND
        for gene in genes[loc[0]]:
            if gene[0]>loc[1]-extend and gene[0]<loc[2]+extend:
                if e[(label,gene[1])]>float(cutoff) and not (label,gene[1]) in tpm_neoloopsLa:
                    f.write(gene[1]+'\t'+label+'\t'+'neoloop'+'\t'+loc[0]+'\t'+str(loc[1]-extend)+'\t'+str(loc[2]+extend)+'\t'+str(e[(label,gene[1])])+'\n')
                    tpm_neoloopsLa.append((label,gene[1]))
                    tpm_neoloops.append(e[(label,gene[1])])
    for loc in locs_loop[label]:
        if loc[2]-loc[1]>10000:
            extend=0
        else:
            extend=EXTEND
        for gene in genes[loc[0]]:
            if  gene[0]>loc[1]-extend and gene[0]<loc[2]+extend:
                if not (label,gene[1]) in tpm_neoloopsLa and not (label,gene[1]) in tpm_loopsLa and e[(label,gene[1])]>float(cutoff):
                    f.write(gene[1]+'\t'+label+'\t'+'loop'+'\t'+loc[0]+'\t'+str(loc[1]-extend)+'\t'+str(loc[2]+extend)+'\t'+str(e[(label,gene[1])])+'\n')
                    tpm_loops.append(e[(label,gene[1])])
                    tpm_loopsLa.append((label,gene[1]))
    for loc in locs_nonloop[label]:
        for gene in genes[loc[0]]:
            if  gene[0]>loc[1] and gene[0]<loc[2]:
                if not (label,gene[1]) in tpm_neoloopsLa and not (label,gene[1]) in tpm_loopsLa and not (label,gene[1]) in tpm_nonloopsLa  and e[(label,gene[1])]>float(cutoff):
                    f.write(gene[1]+'\t'+label+'\t'+'nonloop'+'\t'+loc[0]+'\t'+str(loc[1])+'\t'+str(loc[2])+'\t'+str(e[(label,gene[1])])+'\n')
                    tpm_nonloops.append(e[(label,gene[1])])
                    tpm_nonloopsLa.append((label,gene[1]))
print(stats.mannwhitneyu(tpm_loops,tpm_neoloops)[1])
print(stats.mannwhitneyu(tpm_loops,tpm_nonloops)[1])
print(stats.mannwhitneyu(tpm_neoloops,tpm_nonloops)[1])

print(np.median(tpm_neoloops))
print(np.median(tpm_loops))
print(np.median(tpm_nonloops))
