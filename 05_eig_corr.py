import os
import sys
import os.path
import numpy as np


binsizes=int(sys.argv[1])
rootdir='/home/uni08/txie/analysis/CP/04Compartment/cooltools'
f=open(f'/home/uni08/txie/analysis/CP/04Compartment/summary_results/ASTRY_Eig_{binsizes//1000}k.pcc','w')
def get_filenames(file_dir):
    L=[]
    list = os.listdir(rootdir) #列出文件夹下所有的目录与文件
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]):#判断是否是文件
            L.append(os.path.join(list[i]))
    return L
def overlap(filename1,filename2):
    d1={}
    d2={}
    m=0
    x=[]
    y=[]
    for i in open(rootdir+'/'+filename1):
        t=i.strip().split()
        if len(t)==11:
            if not t[-3]=='E1':
                d1[(t[0],t[1])]=float(t[-3])
    for i in open(rootdir+'/'+filename2):
        t=i.strip().split()
        if len(t)==11:
            if not t[-3]=='E1' :
                d2[(t[0],t[1])]=float(t[-3])
    listname=set(d1.keys())&set(d2.keys())
    for i in listname:
        x.append(d1[i])
        y.append(d2[i])
    pccs=np.corrcoef(x,y)[0][1]
    return str(pccs)
for filename in get_filenames(rootdir):
    if  filename.endswith(f'{binsizes//1000}kb.cis.vecs.tsv') :
        fil=filename.split('_')[0]
        f.write('\t\t'+fil)
f.write('\n')             
for filename1 in get_filenames(rootdir):
    if  filename1.endswith(f'{binsizes//1000}kb.cis.vecs.tsv') :
        fil=filename1.split('_')[0]
        f.write(fil)
        for filename2 in get_filenames(rootdir):
            if  filename2.endswith(f'{binsizes//1000}kb.cis.vecs.tsv'):
                f.write('\t'+overlap(filename1,filename2))
        f.write('\n')
