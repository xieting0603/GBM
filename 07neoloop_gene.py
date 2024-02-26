import sys
import os
import os.path
d={}

#################33This script is for calling all genes associated with neoloops





rootdir='/home/uni08/txie/analysis/GSC/analysis/SV'#####The path of SVs files for all samples
extend=int(sys.argv[1])##########the distance flanking loop anchor to be consider
kind=sys.argv[2]##########Protein coding gene or all genes
f=open('/home/uni08/txie/analysis/GSC/analysis/SV/Mneoloop95_ext'+str(extend//1000)+'k.'+kind+'.gene','w')
sw={}
def get_filenames(file_dir):
    L=[]
    list = os.listdir(file_dir) #列出文件夹下所有的目录与文件
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]):#判断是否是文件
            L.append(os.path.join(list[i]))
    return L

def get_location(filename,label):
    listname=[]
    for i in open(rootdir+'/'+filename):
        t=i.strip().split()
        if ("000,1") in t[-1]:
            if not (t[0],int(t[1]),int(t[2]),label)  in listname:
                listname.append((t[0],int(t[1]),int(t[2]),label))
            if not (t[3],int(t[4]),int(t[5]),label) in listname:
                listname.append((t[3],int(t[4]),int(t[5]),label))
    return list(set(listname))
if kind=='pc':
    for i in open(f"/scratch/users/txie/RNA-seq/0ref/v29.gff"):
        t=i.strip().split()
        if t[-1]=='+':
            d[t[0]]=(t[2],int(t[3]))
        if t[-1]=='-':
            d[t[0]]=(t[2],int(t[4]))
        sw[t[0]]=t[1]
if kind=='all':
    for i in open(f"/scratch/users/txie/RNA-seq/0ref/v29_all.gff"):
        t=i.strip().split()
        if t[-1]=='+':
            d[t[0]]=(t[2],int(t[3]))
        if t[-1]=='-':
            d[t[0]]=(t[2],int(t[4]))
        sw[t[0]]=t[1]

for filename in get_filenames(rootdir):
    if  filename.endswith('95.Mneoloop') and not 'R' in filename and not 'astr' in filename:
        fil=filename.split('.')[0]
        locs=get_location(filename,fil)
        lis=[]
        for gene in d.keys():
            for loc in locs:
               if abs(loc[2]-loc[1])>10000:
                   if d[gene][0]==loc[0] and d[gene][1]>loc[1] and d[gene][1]<loc[2]:
                       f.write(gene+'\t'+sw[gene]+'\t'+fil+'\t'+loc[0]+'\t'+str(loc[1])+'\t'+str(loc[2])+'\n')
                       break
               else:
                   if d[gene][0]==loc[0] and d[gene][1]>loc[1]-extend and d[gene][1]<loc[2]+extend:
                       f.write(gene+'\t'+sw[gene]+'\t'+fil+'\t'+loc[0]+'\t'+str(loc[1]-extend)+'\t'+str(loc[2]+extend)+'\n')
                       break

        
