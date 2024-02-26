import sys
import os
import os.path
##################This script is for calculate the Jaccard-index of SV breakpoints between samples#
mismatch=int(sys.argv[1])
binsizes=5000#################acceptable mismatch bins for loop anchor########

rootdir='.'
f=open(f'./SVPairwise_'+'mis'+str(sys.argv[1])+'_withoutS_index.txt','w')

def get_filenames(file_dir):
    L=[]
    list = os.listdir(rootdir) #列出文件夹下所有的目录与文件
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]):#判断是否是文件
            L.append(os.path.join(list[i]))
    return L
def Ratio(filname1,filename2,B,M):
    d={}
    dd={}
    m=0
    n=0
    C1=0
    C2=0
    loc1=0
    loc2=0
    strand1=0
    strand2=0
    for i in open(rootdir+'/'+filname1):
        t=i.strip().split()
        if t[0].startswith('C'):
            m=m+1
            start=t[1].strip().split(',')
            SV_type=start[0]
            C1=start[1]
            C2=start[4]
            strand1=start[3]
            strand2=start[6]
            loc1=int(start[2])
            loc2=int(start[5])
            if loc1>loc2:
                loc2=int(start[2])
                loc2=int(start[5])
            if C1>C2:
                C1=start[4]
                C2=start[1]
            if not (C1,C2) in d.keys():
                d[(C1,C2)]=[]
            d[(C1,C2)].append((loc1,loc2,strand1,strand2))
    for i in open(rootdir+'/'+filename2):
        t=i.strip().split()
        if t[0].startswith('C'):
            m=m+1
            start=t[1].strip().split(',')
            SV_type=start[0]
            C1=start[1]
            C2=start[4]
            strand1=start[3]
            strand2=start[6]
            loc1=int(start[2])
            loc2=int(start[5])
            if loc1>loc2:
                loc2=int(start[2])
                loc2=int(start[5])
            if C1>C2:
                C1=start[4]
                C2=start[1]
            if not (C1,C2) in dd.keys():
                dd[(C1,C2)]=[]
            dd[(C1,C2)].append((loc1,loc2,strand1,strand2))
    for i in d.keys():
        if i in dd:
            for j in d[i]:
                for k in dd[i]:
                    if abs(j[0]-k[0])<M*B and abs(j[1]-k[1])<M*B:
                        n=n+1
                        dd[i].remove(k)
                        break
                        
    ratio=n/(m-n)
#    return (str(n)+'_'+str(m))
    return str(ratio)

for filename in get_filenames(rootdir):
    if  filename.endswith('_CNV.assemblies.txt') :
        fil=filename.split('_')[0]
        f.write('\t'+fil)
f.write('\n')
for filename1 in get_filenames(rootdir):
    if  filename1.endswith('_CNV.assemblies.txt') :
        fil=filename1.split('_')[0]
        f.write(fil)
        for filename2 in get_filenames(rootdir):
            if  filename2.endswith('_CNV.assemblies.txt'):
                f.write('\t'+Ratio(filename2,filename1,binsizes,mismatch))
        f.write('\n')
