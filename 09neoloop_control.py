import sys
import scipy.stats as stats

control1_1={}
control2_1={}
control1_2={}
control2_2={}
list1=[]
list2=[]
list4=[]
#################The script is for comparing the genes between GSCs samples and astrocyte samples.
rsemdir='/scratch/users/txie/RNA-seq/3rsem_out/'
f=open(sys.argv[1].replace("neoloop95",'Control95'),'w')####The input is the output of neoloop_rna_rsem.py####
f.write("geneID"+'\t'+"log2(tpm+1)"+'\t'+'sample'+'\t'+"type"+'\n')
cutoff=int(sys.argv[1].split('.')[-2])
print(cutoff)
import numpy as np
###########extract tpm of astrs##############
for i in open(f'/scratch/users/txie/RNA-seq/3rsem_out/SRR5442971.genes.results'):
    t=i.strip().split()
    control1_1[t[0].split('_')[0]]=t[5]
for i in open(f'/scratch/users/txie/RNA-seq/3rsem_out/SRR5442972.genes.results'):
    t=i.strip().split()
    control1_2[t[0].split('_')[0]]=t[5]
for i in open(f'/scratch/users/txie/RNA-seq/3rsem_out/SRR5442951.genes.results'):
    t=i.strip().split()
    control2_1[t[0].split('_')[0]]=t[5]
for i in open(f'/scratch/users/txie/RNA-seq/3rsem_out/SRR5442952.genes.results'):
    t=i.strip().split()
    control2_2[t[0].split('_')[0]]=t[5]
for i in open(sys.argv[1]):
    t=i.strip().split()
    if t[2]=='neoloop':
        tpm1=np.log2(float(control1_1[t[0]])/4+float(control1_2[t[0]])/4+float(control2_1[t[0]])/4+float(control2_2[t[0]])/4+1)
        if tpm1>cutoff:
            list1.append(float(t[-1]))
            list2.append(tpm1)

            f.write(t[0]+'\t'+str(t[-1])+'\t'+t[1]+'\t'+'GSCs'+'\n')
            f.write(t[0]+'\t'+str(tpm1)+'\t'+t[1]+'\t'+'Astr'+'\n')
print(stats.wilcoxon(list1,list2)[1])
