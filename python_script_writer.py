import os
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-csv", "--csv", help="path to csv-file")
parser.add_argument("-d", "--work_dir", help="your work directory")
parser.add_argument("-bwp", "--BigWig_pipeline", help="path to bw_pipeline")
parser.add_argument("-t", "--threads", default=1, help="the maximal number of threads")
parser.add_argument("-g", "--reference", help="path to reference genome")
parser.add_argument("-q", "--MAPQ", default=30, help="required MAPQ")
args = parser.parse_args()
df=pd.read_csv(f"{args.csv}")
work_dir=args.work_dir
aligner=args.BigWig_pipeline
threads=args.threads
reference=args.reference
MAPQ=args.MAPQ





struct={}
for i, item in enumerate(df['Sample (е-эмбрион. к-биоптат)']):
    if item in struct:
        link_R1=df['R1'][i]
        link_R2=df['R2'][i]
        struct[item]['count']+=1
        n=struct[item]['count']
        struct[item][f'n{n}_R1']=df['R1'][i]
        struct[item][f'n{n}_R2']=df['R2'][i]

    else:
        link_R1=df['R1'][i]
        link_R2=df['R2'][i]
        struct[item]={}
        struct[item]['count']=1
        struct[item]['n1_R1']=df['R1'][i]
        struct[item]['n1_R2']=df['R2'][i]
data_dir=work_dir+'/data'
try:
    os.mkdir(path=data_dir)
except:
    print('Data directory also exist. Please check files whithin.')
script_dir=work_dir+'/scripts'
try:
    os.mkdir(path=script_dir)
except:
    print('Scripts directory also exist. Please check files whithin.')
bw_dir=data_dir+'/bw'
try:
    os.mkdir(path=bw_dir)
except:
    print('BW directory also exist. Please check files whithin.')

#a way to trick combination of bash and f-string syntaxis
id="{id}"
r='{r}'
for i in struct:
    f=open(script_dir+'/'+i+'.sh', 'w+')
    if struct[i]['count']==1:
        f.write(f'''#!/bin/bash
            
mkdir {data_dir}/{i}
cd {data_dir}/{i}
mkdir fastq
cd {data_dir}/{i}/fastq
wget {struct[i]['n1_R1']}
wget {struct[i]['n1_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"
mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done

bash {aligner} {data_dir}/{i}/fastq/ {reference} {data_dir}/bw {threads} {MAPQ} {i}

                ''')
    elif struct[i]['count']==2:
        f.write(f'''#!/bin/bash
            
mkdir {data_dir}/{i}
mkdir {data_dir}/{i}/{i}_n1
cd {data_dir}/{i}/{i}_n1
mkdir fastq
cd {data_dir}/{i}/{i}_n1/fastq
wget {struct[i]['n1_R1']}
wget {struct[i]['n1_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"
mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done

mkdir {data_dir}/{i}/{i}_n2
cd {data_dir}/{i}/{i}_n2
mkdir fastq
cd {data_dir}/{i}/{i}_n2/fastq
wget {struct[i]['n2_R1']}
wget {struct[i]['n2_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"

mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done


cd
mkdir {data_dir}/{i}/fastq
zcat {data_dir}/{i}/{i}_n1/fastq/*R1*.fastq.gz {data_dir}/{i}/{i}_n2/fastq/*R1*.fastq.gz|gzip > {data_dir}/{i}/fastq/{i}_R1.fastq.gz
zcat {data_dir}/{i}/{i}_n1/fastq/*R2*.fastq.gz {data_dir}/{i}/{i}_n2/fastq/*R2*.fastq.gz|gzip > {data_dir}/{i}/fastq/{i}_R2.fastq.gz

bash {aligner} {data_dir}/{i}/fastq/ {reference} {data_dir}/bw {threads} {MAPQ} {i}

                               ''')
    elif struct[i]['count']==4:
        f.write(f'''#!/bin/bash
            
mkdir {data_dir}/{i}
mkdir {data_dir}/{i}/{i}_n1
cd {data_dir}/{i}/{i}_n1
mkdir fastq
cd {data_dir}/{i}/{i}_n1/fastq
wget {struct[i]['n1_R1']}
wget {struct[i]['n1_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"
mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done        
mkdir {data_dir}/{i}/{i}_n2
cd {data_dir}/{i}/{i}_n2
mkdir fastq
cd {data_dir}/{i}/{i}_n2/fastq
wget {struct[i]['n2_R1']}
wget {struct[i]['n2_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"
mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done

mkdir {data_dir}/{i}
mkdir {data_dir}/{i}/{i}_n3
cd {data_dir}/{i}/{i}_n3
mkdir fastq
cd {data_dir}/{i}/{i}_n3/fastq
wget {struct[i]['n3_R1']}
wget {struct[i]['n3_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"
mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done

mkdir {data_dir}/{i}
mkdir {data_dir}/{i}/{i}_n4
cd {data_dir}/{i}/{i}_n4
mkdir fastq
cd {data_dir}/{i}/{i}_n4/fastq
wget {struct[i]['n4_R1']}
wget {struct[i]['n4_R2']}
r=1
for file in *
do
if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
then
id="$( cut -d '.' -f1 <<< "$file" )"
new_name="${id}R${r}.fastq.gz"
mv "$file" "$new_name"
fi
if [ "$r" = 1 ]
then
r=2
else
r=1
fi
done

cd
mkdir {data_dir}/{i}/fastq
zcat {data_dir}/{i}/{i}_n1/fastq/*R1*.fastq.gz {data_dir}/{i}/{i}_n2/fastq/*R1*.fastq.gz {data_dir}/{i}/{i}_n3/fastq/*R1*.fastq.gz {data_dir}/{i}/{i}_n4/fastq/*R1*.fastq.gz|gzip > {data_dir}/{i}/fastq/{i}_R1.fastq.gz
zcat {data_dir}/{i}/{i}_n1/fastq/*R2*.fastq.gz {data_dir}/{i}/{i}_n2/fastq/*R2*.fastq.gz {data_dir}/{i}/{i}_n3/fastq/*R2*.fastq.gz {data_dir}/{i}/{i}_n4/fastq/*R2*.fastq.gz|gzip > {data_dir}/{i}/fastq/{i}_R2.fastq.gz

bash {aligner} {data_dir}/{i}/fastq/ {reference} {data_dir}/bw {threads} {MAPQ} {i}
                               ''')
    else:
        print("There is a problem with count of files pairs. Please check, that there are 1, 2 or 4 file pairs in table")