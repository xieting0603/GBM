import sys
import os
import os.path
###############################################To convert neoloop CNV segment to circos format CNV## new format CNV outputs are located in the same folder with old##########
inputdir=sys.argv[1]

##########################'/home/uni08/txie/analysis/Pan/08trans/cnv'



def get_filenames(file_dir):
    L=[]
    list = os.listdir(inputdir) #列出文件夹下所有的目录与文件
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]):#判断是否是文件
            L.append(os.path.join(list[i]))
    return L

for filename in get_filenames(inputdir):
    if  filename.endswith('_25kb.2') :
        print(filename)
        fil=filename.split('.')[0].strip('_25kb')
        f_de=open(f'{inputdir}/{fil}_25kb.3','w')
        for i in open(f'{inputdir}/{filename}'):
            t=i.strip().split('\t')
            f_de.write('hs'+t[0].strip('chr')+'\t'+t[1]+'\t'+t[2]+'\t'+str(int(t[3])-2)+'\n')
