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
wget {sample_data[f'n1_R1']} -O {sample_name}_R1.fastq.gz
wget {sample_data[f'n1_R2']} -O {sample_name}_R2.fastq.gz

cd

''')

    if sample_data['count'] > 1:
      for n in range(1, sample_data['count'] + 1):
            f.write(f'''# Processing {sample_name}_n{n}
cd {data_dir}/{sample_name}
mkdir -p {sample_name}_n{n}/fastq
cd {data_dir}/{sample_name}/{sample_name}_n{n}/fastq
wget {sample_data[f'n{n}_R1']} -O {sample_name}_n{n}_R1.fastq.gz
wget {sample_data[f'n{n}_R2']} -O {sample_name}_n{n}_R2.fastq.gz 
''')
            if n==sample_data['count']:
                f.write(f'''
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
    parser.add_argument("-ch", "--checking", type=str, default=True, help="Checking the correctness of links.")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    checking=args.checking
    work_dir = args.work_dir
    aligner = args.BigWig_pipeline
    threads = args.threads
    reference = args.reference
    MAPQ = args.MAPQ

    struct = {}
    for i, item in enumerate(df['Sample']):
        if item not in struct:
            struct[item] = {'count': 0}
        
        struct[item]['count'] += 1
        count = struct[item]['count']
        struct[item][f'n{count}_R1'] = df['R1'][i]
        struct[item][f'n{count}_R2'] = df['R2'][i]
        if checking:
            assert struct[item][f'n{count}_R1']!=struct[item][f'n{count}_R2'], f'Tere are the same R1 and R2 links for {item} with {df["Sequenced Read Pairs"][i]} Sequenced Read Pairs.'
            assert struct[item][f'n{count}_R1'].split('/')[-1].translate(str.maketrans('', '', '12'))==struct[item][f'n{count}_R2'].split('/')[-1].translate(str.maketrans('', '', '12')), f'Tere is wrong R1 or R2 link for {item} with {df["Sequenced Read Pairs"][i]} Sequenced Read Pairs.'

    data_dir = os.path.join(work_dir, 'data')
    create_directory(data_dir)
    
    script_dir = os.path.join(work_dir, 'scripts')
    create_directory(script_dir)

    bw_dir = os.path.join(data_dir, 'bw')
    create_directory(bw_dir)

    for sample_name, sample_data in struct.items():
        script_path = os.path.join(script_dir, f'{sample_name}.sh')
        with open(script_path, 'w+') as f:
            write_script(f, aligner, data_dir, reference, threads, MAPQ, sample_name, sample_data)

if __name__ == "__main__":
    main()
