import sys
import os
import os.path
pos=[]
chrsize={}
coolpath=sys.argv[1]
svpath=sys.argv[2]########location with SV/CNV created by EagelC files###########
outpath=sys.argv[3] 


 ######output folder################3
def get_filenames(file_dir):
    L=[]
    list = os.listdir(file_dir) #列出文件夹下所有的目录与文件
    for i in range(0,len(list)):
        if os.path.isfile(file_dir+'/'+list[i]):#判断是否是文件
            L.append(os.path.join(list[i]))
    return L
    
for i in open('/scratch1/users/txie/data/hg38/hg38.chrom.sizes'):
    t=i.strip().split()
    chrsize[t[0]]=int(t[1])

for filename in get_filenames(svpath):
    listchrpairs=[]
    N=0
    if  ('combine')  in filename and filename.endswith('txt'):
        fil=filename.split('.')[0]
        print(fil)
        for i in open(svpath+'/'+filename):
            t=i.strip().split()
            if not t[0].startswith('chrom'):
                N=N+1
                if not t[0]==t[2] :
                    os.system('plot-interSVs --cool-uri '+coolpath+'/'+fil+\
                    '-Arima-allReps-filtered.mcool::/resolutions/1000000    --full-sv-file \
                    '+svpath+'/'+filename+' --output-figure-name  \
                    '+outpath+'/'+fil+'.'+str(N)+'.pdf  -C '+t[0]+' '+t[2]+' --balance-type Raw --dpi 1000')
                elif abs(int(t[3])-int(t[1])<10000000) :#####skip very large intraSVs, as it is very time and memory consuming
                    extend=int(t[3])-int(t[1])###### decide the flanking size along SVs
                    if extend<300000:
                        extend=300000
                    elif extend>2000000:
                        extend=1000000
                    show_start=int(t[1])-extend
                    show_end=int(t[3])+extend
                    if show_start<0:
                        show_start=0
                    if show_end>chrsize[t[0]]:
                        show_end=chrsize[t[0]]
                    os.system('plot-intraSVs  --cnv-max-value  4 --cool-uri '+coolpath+'/'+fil+\
                    '-Arima-allReps-filtered.mcool::/resolutions/5000 --full-sv-file '
                    +svpath+'/'+filename+' --output-figure-name  '+outpath+'/'+fil+'.'\
                    +str(N)+'.pdf   --coordinates-to-display '+t[1]+' '+t[3]+' --balance-type CNV --dpi 1000  --cnv-file '\
                    +svpath+'/'+fil+'_5kb.bw   --region '+t[0]+':'+str(show_start)+'-'+str(show_end))
                    
                    
