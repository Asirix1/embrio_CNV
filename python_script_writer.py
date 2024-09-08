import os
import pandas as pd
import argparse

def create_directory(path):
    """Create directory if it does not exist."""
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        print(f"{path} directory already exists. Please check the files within.")

def write_script(f, aligner, data_dir, reference, threads, MAPQ, sample_name, sample_data):
    """Write the appropriate bash script based on the count of paired files."""
    f.write(f'''#!/bin/bash
mkdir -p {data_dir}/{sample_name}
cd {data_dir}/{sample_name}
mkdir -p fastq

''')
    
    if sample_data['count']==1:
        f.write(f'''# Processing {sample_name}
cd {data_dir}/{sample_name}/fastq
wget {sample_data[f'n1_R1']}
wget {sample_data[f'n1_R2']}

r=1
for file in *
do
    if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
    then
        id="$( cut -d '.' -f1 <<< "$file" )"
        new_name="${{id}}R${{r}}.fastq.gz"
        mv "$file" "$new_name"
    fi
    r=$((3-r)) # Toggle between 1 and 2
done

cd

''')

    if sample_data['count'] > 1:
        for n in range(1, sample_data['count'] + 1):
            f.write(f'''# Processing {sample_name}_n{n}
cd {data_dir}/{sample_name}
mkdir -p {sample_name}_n{n}/fastq
cd {data_dir}/{sample_name}/{sample_name}_n{n}/fastq
wget {sample_data[f'n{n}_R1']}
wget {sample_data[f'n{n}_R2']}
r=1
for file in *
do
    if [[ "$file" != *R1*.fastq.gz && "$file" != *R2*.fastq.gz ]]
    then
        id="$( cut -d '.' -f1 <<< "$file" )"
        new_name="${{id}}R${{r}}.fastq.gz"
        mv "$file" "$new_name"
    fi
    r=$((3-r)) # Toggle between 1 and 2
done

cd
zcat {" ".join([f"{data_dir}/{sample_name}/{sample_name}_n{n}/fastq/*R1*.fastq.gz" for n in range(1, sample_data['count'] + 1)])} | gzip > {data_dir}/{sample_name}/fastq/{sample_name}_R1.fastq.gz
zcat {" ".join([f"{data_dir}/{sample_name}/{sample_name}_n{n}/fastq/*R2*.fastq.gz" for n in range(1, sample_data['count'] + 1)])} | gzip > {data_dir}/{sample_name}/fastq/{sample_name}_R2.fastq.gz
''')

    f.write(f'bash {aligner} {data_dir}/{sample_name}/fastq/ {reference} {data_dir}/bw {threads} {MAPQ} {sample_name}\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-csv", "--csv", help="path to csv-file", required=True)
    parser.add_argument("-d", "--work_dir", help="your work directory", required=True)
    parser.add_argument("-bwp", "--BigWig_pipeline", help="path to bw_pipeline", required=True)
    parser.add_argument("-t", "--threads", type=int, default=1, help="the maximal number of threads")
    parser.add_argument("-g", "--reference", help="path to reference genome", required=True)
    parser.add_argument("-q", "--MAPQ", type=int, default=30, help="required MAPQ")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    work_dir = args.work_dir
    aligner = args.BigWig_pipeline
    threads = args.threads
    reference = args.reference
    MAPQ = args.MAPQ

    struct = {}
    for i, item in enumerate(df['Sample (е-эмбрион. к-биоптат)']):
        if item not in struct:
            struct[item] = {'count': 0}
        
        struct[item]['count'] += 1
        count = struct[item]['count']
        struct[item][f'n{count}_R1'] = df['R1'][i]
        struct[item][f'n{count}_R2'] = df['R2'][i]

    data_dir = os.path.join(work_dir, 'data')
    create_directory(data_dir)
    
    script_dir = os.path.join(work_dir, 'scripts')
    create_directory(script_dir)

    bw_dir = os.path.join(data_dir, 'bw')
    create_directory(bw_dir)

    for sample_name, sample_data in struct.items():
        script_path = os.path.join(script_dir, f'{sample_name}.sh')
        with open(script_path, 'w+') as f:
            if sample_data['count'] not in [1, 2, 4]:
                print(f"Sample {sample_name} has {sample_data['count']} file pairs. Please ensure there are 1, 2, or 4 file pairs.")
            else:
                write_script(f, aligner, data_dir, reference, threads, MAPQ, sample_name, sample_data)

if __name__ == "__main__":
    main()
