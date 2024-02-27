import sys
gene={}
###This file is to classify EE/EP/PP/NN neoloops###########
cutoff=sys.argv[1]#############cutoff of peaks############
for i in open('/scratch/users/txie/RNA-seq/0ref/v29_all.gff'):
    t=i.strip().split()
    if not t[2] in gene.keys():
        gene[t[2]]=[]
    if t[-1]=='+':
        loc=t[-3]
    else:
        loc=t[-2]
    gene[t[2]].append(int(loc))
################d
def peaks_overlap(key,d,loc1,loc2,ext,value):
    ori='N'
    for i in d[key]:
      if min(i[1],loc2+ext)-max(i[0],loc1-ext)>0:
          ori=value
          break
    return ori

def gene_overlap(key,d,loc1,loc2,ext,value):
    ori='N'
    for i in d[key]:
      if min(i+1,loc2+ext)-max(i-1,loc1-ext)>0:
          ori=value
          break
    return ori
list=['GSC275_H3K27AC','GSC275B_H3K27AC','GSC486_H3K27AC','GSC428_H3K27AC','GSC402_H3K27AC','GSC412_H3K27AC','GSC452C_H3K27AC_2','GSC394B_H3K27AC','GSC148_H3K27AC','GSC450_H3K27AC_2']
print(''+'\t'+'EE'+'\t'+'EP'+'\t'+'PP'+'\t'+'n')
for name in list:
    nameloop=name.split('_')[0]
    fout=open(f'./{nameloop}.95.Mneoloop_{cutoff}.EP','w')
    peak={}
    Nneo=0
    EPneo=0
    EEneo=0
    PPneo=0
    summneo=0
    N=0
    EP=0
    EE=0
    PP=0
    summ=0
    for i in open (f'/home/uni08/txie/analysis/CP/07Chips/cuttag/{name}_results/peakCalling/{name}_seacr_top{cutoff}.peaks.stringent.bed'):
        t=i.strip().split()
        if not t[0] in peak.keys():
            peak[t[0]]=[]
        peak[t[0]].append((int(t[1]),int(t[2])))
    for i in open(f'./{nameloop}.95.Mneoloop'):
        if not i.startswith('#'):
            t=i.strip().split()
            chr1=t[0]
            chr2=t[3]
            a1=''
            a2=''
            p1=''
            p2=''
            if int(t[2])-int(t[1])>10000:
                ext=0
            else:
                ext=5000
            a1=peaks_overlap(chr1,peak,int(t[1]),int(t[2]),ext,'E')
            a2=peaks_overlap(chr2,peak,int(t[4]),int(t[5]),ext,'E')
            p1=gene_overlap(chr1,gene,int(t[1]),int(t[2]),ext,'P')
            p2=gene_overlap(chr2,gene,int(t[4]),int(t[5]),ext,'P')
            if (t[-1]).endswith('00,1'):
                summneo+=1
                if a1=='E' and a2=='E':
                    EEneo+=1
                    fout.write(i.strip()+'\t'+'EE'+'\n')
                if p1=='P' and p2=='P':
                    PPneo+=1
                    fout.write(i.strip()+'\t'+'PP'+'\n')
                if a1=='E' and p2=='P' :
                    EPneo+=1
                    fout.write(i.strip()+'\t'+'EP'+'\n')
                if p1=='P' and a2=='E':
                    EPneo+=1
                    fout.write(i.strip()+'\t'+'PE'+'\n')
                if (a1=='N' and p1=='N') or (a2=='N' and p2=='N'):
                    fout.write(i.strip()+'\t'+'--'+'\n')
                    Nneo+=1
            if (t[-1]).endswith('00,0'):
                summ+=1
                if a1=='E' and a2=='E':
                    EE+=1
                    fout.write(i.strip()+'\t'+'EE'+'\n')
                if p1=='P' and p2=='P':
                    PP+=1
                    fout.write(i.strip()+'\t'+'PP'+'\n')
                if (a1=='E' and p2=='P' ):
                    EP+=1
                    fout.write(i.strip()+'\t'+'EP'+'\n')
                if p1=='P' and a2=='E':
                    EP+=1
                    fout.write(i.strip()+'\t'+'PE'+'\n')
                if (a1=='N' and p1=='N') or (a2=='N' and p2=='N'):
                    fout.write(i.strip()+'\t'+'--'+'\n')
                    N+=1
    print(nameloop+'\t'+str(EEneo)+'\t'+str(EPneo)+'\t'+str(PPneo)+'\t'+str(Nneo)+'\t'+str(summneo))
    print(nameloop+'\t'+str(EE)+'\t'+str(EP)+'\t'+str(PP)+'\t'+str(N)+'\t'+str(summ))
