import sys
import os
import os.path


#############################################to convert EagleC SV to circos format SV### new format SV outputs are located in the same folder with old####################
inputdir='/home/uni08/txie/analysis/Pan/analysis/SV'

##############################################################
def get_filenames(file_dir):
    L=[]
    list = os.listdir(inputdir) 
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]):
            L.append(os.path.join(list[i]))
    return L

def start_end(A,B):
    if int(A) < int(B):
        Str = A+'\t'+B
    else:
        Str = B+'\t'+A
    return Str
    
def col(C):
    if C=='translocation':
#       color='set3-4-qual-3'
        color='(124,111,176)'
    if C=='duplication':
#       color='set3-4-qual-4'
        color='(255,2,43)'
    if C=='deletion':
#       color='set3-4-qual-1'
        color='(160,217,246)'
    if C=='inversion':
        color='(191,191,191)'
#       color='493657'
    return color

for filename in get_filenames(inputdir):
    if  filename.endswith('.CNN_SVs.5K_combined.merge') :
        fil=filename.split('.')[0]
        f_de=open(f'{inputdir}/{fil}.forcircos.type','w')
        for i in open(f'{inputdir}/{filename}'):
            t=i.strip().split('\t')
            f_de.write('hs'+t[0].strip('chr')+'\t'+str(t[3])+'\t'+str(int(t[3])+1)+'\t'+'hs'+t[1].strip('chr')+'\t'+str(t[4])+'\t'+str(int(t[4])+1)+'\tcolor='+col(t[-1])+'\n')
    