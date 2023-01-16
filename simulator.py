import argparse
import time
parser = argparse.ArgumentParser(description="Simulate reads with kde model",add_help=True)
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('--kde_model_path', type=str, help="Path of KDE model",required=True)
requiredNamed.add_argument('--insertion_rate', type=float,default=0.01, help="Insertion rate",required=False)
requiredNamed.add_argument('--deletion_rate', type=float,default=0.01, help="Deletion rate",required=False)
requiredNamed.add_argument('--substitution_rate', type=float,default=0.01, help="Substitution rate",required=False)
requiredNamed.add_argument('--coord_error_in_5_end', type=float,default=0.05, help="Coord randomness in 5' end",required=False)
requiredNamed.add_argument('--coord_error_in_3_end', type=float,default=0.05, help="Coord randomness in 3' end",required=False)
requiredNamed.add_argument('-expr','--expression_profile', type=str, help="Expression profile",required=True)
requiredNamed.add_argument('-cdna','--reference_transcriptome', type=str, help="Reference transcriptome",required=True)
# requiredNamed.add_argument('-gtf','--reference_annotation', type=str, help="Reference annnotation",required=True)
# requiredNamed.add_argument('-prot','--library_protocol', type=str, default='direct_RNA',help="Library protocol",required=False)
requiredNamed.add_argument('--num_reads', type=int, help="The number of simulated reads",required=True)
requiredNamed.add_argument('-t','--threads', type=int, help="Threads",required=False)
requiredNamed.add_argument('-o','--output', type=str, help="The path of output simulated reads",required=True)
import dill as pickle
import numpy as np
# import sys
import os
import subprocess
import shutil

# ref_file_path = '/fs/project/PCON0009/Au-scratch2/haoran/reference/genecode/gencode.v38.annotation.gtf'
# READ_LEN=0
# READ_JUNC_MIN_MAP_LEN= 1
# short_read_alignment_file_path = None

# threads = 3
# [gene_exons_dict,gene_points_dict, gene_isoforms_dict,_,
#             _, gene_regions_dict, gene_isoforms_length_dict,_,_,_] = parse_annotation(ref_file_path,threads)
import re
import pysam
from Bio.Seq import Seq
from pathlib import Path
import pandas as pd
import concurrent.futures
def extract_error_rate(err_sub,err_ins,err_del):
    error_type = ["no","sub","ins","del"]
    error_prob = [(1.0-(err_sub + err_ins + err_del)),err_sub,err_ins,err_del]
    return error_type,error_prob
def mutate_read(read_seq,error_type,error_prob):
    read_seq = read_seq.upper()
    dic_error_no = {"A":"A","C":"C","G":"G","T":"T"}
    dic_error_sub = {"A":np.random.choice(["C","G","T"]),"C":np.random.choice(["A","G","T"]),"G":np.random.choice(["A","C","T"]),"T":np.random.choice(["A","C","G"])}
    dic_error_ins = {"A":"A"+np.random.choice(["A","C","G","T"]),"C":"C"+np.random.choice(["A","C","G","T"]),"G":"G"+np.random.choice(["A","C","G","T"]),"T":"T"+np.random.choice(["A","C","G","T"])}
    dic_error_del = {"A":"","C":"","G":"","T":""}

    dic_error = {"no":dic_error_no,"sub":dic_error_sub,"ins":dic_error_ins,"del":dic_error_del}
    read_seq_new = ""
    for nuc in read_seq:
        nuc_new = dic_error[np.random.choice(error_type,p=error_prob)][nuc]
        read_seq_new += nuc_new
    return read_seq_new
import random
def get_start_end_kde(kde, num):
    return kde.sample(num)
def select_nearest_kde3d(sampled_3d, isoform_len):
    fc = sampled_3d[:, 0]
    dist = np.abs(fc - isoform_len)
    indices = np.where(dist == dist.min())[0]
    idx = np.random.choice(indices, 1)[0]
    return int(sampled_3d[idx][1]),int(sampled_3d[idx][2])
def randomize_offset(ref_5_end_offset,ref_3_end_offset,error_rate_5_end,error_rate_3_end,isoform_len):
    MIN_READ_LEN = 15
    new_ref_5_end_offset = ref_5_end_offset + int(np.round((random.random() - 0.5) * error_rate_5_end * isoform_len))
    new_ref_3_end_offset = ref_3_end_offset + int(np.round((random.random() - 0.5) * error_rate_3_end * isoform_len))
    if new_ref_5_end_offset < 0:
        new_ref_5_end_offset = ref_5_end_offset
    if new_ref_3_end_offset < 0:
        new_ref_3_end_offset = ref_3_end_offset
    if new_ref_5_end_offset + new_ref_3_end_offset + MIN_READ_LEN >= isoform_len:
        new_ref_5_end_offset = ref_5_end_offset
        new_ref_3_end_offset = ref_3_end_offset
    return new_ref_5_end_offset,new_ref_3_end_offset
def get_sequence(start,end,isoform_name,reference_transcriptome,error_type,error_prob,original_isoform_name):
    sequence = Seq(reference_transcriptome.fetch(original_isoform_name,start - 1,end)).upper()
    mutated_sequence = mutate_read(str(sequence),error_type,error_prob)
    return mutated_sequence
def simulated_reads_for_isoform(kde,reference_transcriptome,error_type,error_prob,error_rate_5_end,error_rate_3_end,args):
    num_reads,ref_trx,ref_trx_len,strand,original_isoform_name = args
    output_str = ''
    
    for i in range(num_reads):
        sampled_3d = get_start_end_kde(kde, 100000)
        ref_start_offset,ref_end_offset = select_nearest_kde3d(sampled_3d, ref_trx_len)
        if error_rate_5_end !=0 or error_rate_3_end !=0:
            ref_start_offset,ref_end_offset = randomize_offset(ref_start_offset,ref_end_offset,error_rate_5_end,error_rate_3_end,ref_trx_len)
        # if strand == '+':
           
        # else:
        #     ref_end_offset,ref_start_offset = select_nearest_kde3d(sampled_3d, ref_trx_len)
        #     if error_rate_5_end !=0 or error_rate_3_end !=0:
        #         ref_end_offset,ref_start_offset = randomize_offset(ref_end_offset,ref_start_offset,error_rate_5_end,error_rate_3_end,ref_trx_len)
        # 1 based coord
        ref_start_pos,ref_end_pos = ref_start_offset + 1,ref_trx_len - ref_end_offset

        if ref_start_pos < 0 or ref_end_pos < ref_start_pos:
            continue
        output_ref_trx = ref_trx.split('.')[0]
        output_ref_start_pos = ref_start_pos - 1
        read_len = ref_end_pos - ref_start_pos + 1
        # if strand == '+':
        #     direction = 'F'
        # else:
        #     direction = 'R'
        direction = 'F'
        read_name = f'>{output_ref_trx}_{output_ref_start_pos}_aligned_{i}_{direction}_{ref_end_pos}_{read_len}_{i}'
        try:
            read_sequence = get_sequence(ref_start_pos,ref_end_pos,ref_trx,reference_transcriptome,error_type,error_prob,original_isoform_name)
        except Exception as e:
            print(e)
            break
        fasta_lines = read_name + '\n' + read_sequence + '\n'
        output_str += fasta_lines
    return output_str
def submit_worker(output_dir,kde_model_path,cDNA_path,error_type,error_prob,error_rate_5_end,error_rate_3_end,list_of_args,worker_id):
    buffer_size = 1024 * 1024 * 1024
    output_str = ''
    num_reads = 0
    temp_dir = f'{output_dir}/temp/'
    Path(temp_dir).mkdir(exist_ok=True,parents=True)
    reference_transcriptome = pysam.FastaFile(cDNA_path)
    with open(kde_model_path,'rb') as f:
        kde = pickle.load(f)
#     num_reads = isoform_num_reads_dict[ref_trx]
    with open(f'{output_dir}/temp/simulated_reads_{worker_id}.fasta','w') as f:
        for args in list_of_args:
            num_reads += args[0]
            isoform_output_str = simulated_reads_for_isoform(kde,reference_transcriptome,error_type,error_prob,error_rate_5_end,error_rate_3_end,args)
            output_str += isoform_output_str
            if len(output_str.encode('utf-8')) > buffer_size:
                f.write(output_str)
                print(f'Worker {worker_id} simulated {num_reads} reads',flush=True)
                num_reads = 0
                output_str = ''
        if output_str != '':
            f.write(output_str)
def submit_worker_helper(args):
    output_dir,kde_model_path,cDNA_path,error_type,error_prob,error_rate_5_end,error_rate_3_end,list_of_args,worker_id = args
    submit_worker(output_dir,kde_model_path,cDNA_path,error_type,error_prob,error_rate_5_end,error_rate_3_end,list_of_args,worker_id)
    return worker_id
def sample_reads(sampling_rates,total_num_reads,variance=False,no_sampling=False):
    if no_sampling:
        read_counts = np.around(sampling_rates* total_num_reads).astype(int)
        read_counts[read_counts<0] = 0
        return read_counts
    if variance:
        variance_sampling_rates = sampling_rates + 0.25 * sampling_rates * np.random.normal(0, 1.0, sampling_rates.shape[0])
        variance_sampling_rates[variance_sampling_rates<0] = 0
        variance_sampling_rates = variance_sampling_rates/variance_sampling_rates.sum()
        return np.random.multinomial(total_num_reads, variance_sampling_rates)
    else:
        return np.random.multinomial(total_num_reads, sampling_rates)
def main():
    st = time.time()
    args = parser.parse_args()
    error_rate_5_end = args.coord_error_in_5_end
    error_rate_3_end = args.coord_error_in_3_end
    insertion_rate = args.insertion_rate
    deletion_rate = args.deletion_rate
    substitution_rate = args.substitution_rate
    num_reads = args.num_reads
    expression_profile = args.expression_profile
    cDNA_path = args.reference_transcriptome
    kde_model_path = args.kde_model_path
    output_dir = args.output
    if args.threads is None:
        threads = len(os.sched_getaffinity(0))
    else:
        threads = args.threads
    error_type,error_prob = extract_error_rate(substitution_rate,insertion_rate,deletion_rate)
    
    # strand_dict = {}
    isoform_length_dict = {}
    # with open(ref_file_path,'r') as f:
    #     for line in f:
    #         if line.lstrip()[0] == "#":
    #             continue
    #         fields = line.split('\t')
    #         if (fields[2] != 'exon'):
    #             continue
    #         strand = fields[6]
    #         gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
    #         chr_name = fields[0]
    #         isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
    #         strand_dict[isoform_name] = strand
    #         start_pos = int(fields[3])
    #         end_pos = int(fields[4])
    #         if isoform_name not in isoform_length_dict:
    #             isoform_length_dict[isoform_name] = 0
    #         isoform_length_dict[isoform_name] += end_pos - start_pos + 1
    isoform_name_mapping = {}
    reference_transcriptome = pysam.FastaFile(cDNA_path)
    for isoform_name,isoform_length in zip(reference_transcriptome.references,reference_transcriptome.lengths):
        isoform_name_mapping[isoform_name.split('|')[0]] = isoform_name
        isoform_length_dict[isoform_name.split('|')[0]] = isoform_length
    # isoform_set = set(isoform_length_dict.keys())
    isoform_set = set(isoform_name_mapping.keys())
    df = pd.read_csv(expression_profile,sep='\t',skiprows=1,header=None)
    df.columns = ['target_id','est_counts','tpm']
    df = df.set_index('target_id')
    isoform_expr_index = df.index.intersection(isoform_set)
    if len(isoform_expr_index) < len(df.index):
        error_set = set(df.index)-set(isoform_expr_index)
        print(f'Isoforms {error_set} in the expression profile does not exists in the reference annotation!')
    df = df.loc[isoform_expr_index]
    df['sampling_rates'] = df['tpm']/df['tpm'].sum()
    df['Read_count'] = sample_reads(df['sampling_rates'].values,num_reads)
    df.index.name = 'target_id'
    df = df.sort_values('Read_count',ascending=False).reset_index()
    df['Worker_id'] = df.index.map(lambda x:x%threads)
    df = df.sort_values(by='Worker_id').set_index('target_id')
    isoform_num_reads_dict = df['Read_count'].to_dict()
    list_of_args = []
    for isoform,num_reads in isoform_num_reads_dict.items():
        if num_reads > 0:
            strand = '+'
            original_isoform_name = isoform_name_mapping[isoform]
            isoform_len = isoform_length_dict[isoform]
            args = num_reads,isoform,isoform_len,strand,original_isoform_name
            list_of_args.append(args)
    all_worker_id = []
    if threads > 1:
        chunksize, extra = divmod(len(list_of_args), threads)
        if extra:
            chunksize += 1
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = []
            for i in range(threads - 1):
                submit_args = output_dir,kde_model_path,cDNA_path,error_type,error_prob,error_rate_5_end,error_rate_3_end,list_of_args[i*chunksize:(i+1)*chunksize],i+1
                futures.append(executor.submit(submit_worker_helper,submit_args))
            submit_args = output_dir,kde_model_path,cDNA_path,error_type,error_prob,error_rate_5_end,error_rate_3_end,list_of_args[(threads-1)*chunksize:],threads
            futures.append(executor.submit(submit_worker_helper,submit_args))
            for fut in concurrent.futures.as_completed(futures):
                worker_id = fut.result()
                all_worker_id.append(worker_id)
    else:
        submit_worker(output_dir,kde_model_path,cDNA_path,error_type,error_prob,error_rate_5_end,error_rate_3_end,list_of_args,1)
        all_worker_id.append(1)

    with open(f'{output_dir}/simulated_reads.fasta','wb') as wfd:
        for worker_id in all_worker_id:
            try:
                with open(f'{output_dir}/temp/simulated_reads_{worker_id}.fasta','rb') as fd:
                    shutil.copyfileobj(fd, wfd)
            except Exception as e:
                print(e)
                continue
    try:
        shutil.rmtree(f'{output_dir}/temp/')
    except OSError as e:
        print(e)
    elapsed = time.time() - st
    print(f"DONE in {elapsed} s")
main()
    
        
