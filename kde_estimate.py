import argparse
from pathlib import Path
parser = argparse.ArgumentParser(description="Estimate weights",add_help=True)
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('--kde_model_path', type=str, help="Path of KDE model",required=True)
requiredNamed.add_argument('--rname', type=str, help="Which chormosome",required=True)
requiredNamed.add_argument('--annotation', type=str, help="GTF annotation",required=True)
requiredNamed.add_argument('-o','--output', type=str, help="The path of output weight",required=True)
args = parser.parse_args()
output_dir = args.output
kde_model_path = args.kde_model_path
selected_rname = args.rname
# temp_dir = f"{output_dir}_temp/"
Path(output_dir).mkdir(exist_ok=True,parents=True)
import dill as pickle
import numpy as np
import sys
from tqdm import tqdm
import os
sys.path.append('/fs/project/PCON0009/Au-scratch2/haoran/_projects/isoform_quantification/isoform_quantification')
from parse_annotation import parse_annotation
ref_file_path = args.annotation
READ_LEN=0
READ_JUNC_MIN_MAP_LEN= 1
short_read_alignment_file_path = None
threads = len(os.sched_getaffinity(0))
print(f'Using {threads} threads',flush=True)
[gene_exons_dict,gene_points_dict, gene_isoforms_dict,_,
            _, gene_regions_dict, gene_isoforms_length_dict,_,_,_] = parse_annotation(ref_file_path,threads)
import re
strand_dict = {}
with open(ref_file_path,'r') as f:
    for line in f:
        if line.lstrip()[0] == "#":
            continue
        fields = line.split('\t')
        if (fields[2] != 'exon'):
            continue
        strand = fields[6]
        gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
        strand_dict[gene_name] = strand
def get_start_end_kde(kde, num):
    return kde.sample(num)
def select_nearest_kde3d(sampled_3d, isoform_len):
    fc = sampled_3d[:, 0]
    dist = np.abs(fc - isoform_len)
    indices = np.where(dist == dist.min())[0]
    idx = np.random.choice(indices, 1)[0]
    return int(sampled_3d[idx][1]),int(sampled_3d[idx][2])
def simulated_reads_for_gene(kde,num_reads,gene_isoforms_length_dict,strand_dict,rname,gname):
    simulation_reads = {}
    for ref_trx, ref_trx_len in gene_isoforms_length_dict[rname][gname].items():
        simulation_reads[ref_trx] = []
        for i in range(num_reads):
            sampled_3d = get_start_end_kde(kde, 100000)
            if strand_dict[gname] == '+':
                ref_start_offset,ref_end_offset = select_nearest_kde3d(sampled_3d, ref_trx_len)
            else:
                ref_end_offset,ref_start_offset = select_nearest_kde3d(sampled_3d, ref_trx_len)
            ref_start_pos,ref_end_pos = ref_start_offset,ref_trx_len - ref_end_offset
            if ref_start_pos < 0 or ref_end_pos < ref_start_pos:
                continue
            simulation_reads[ref_trx].append([ref_start_pos,ref_end_pos])
    return simulation_reads
def assign_reads_to_region(simulation_reads,gene_isoforms_length_dict,gene_points_dict,gene_exons_dict,strand_dict,rname,gname):
    weight_dict = {}
    for isoform in simulation_reads:
        weight_dict[isoform] = {}
        isoform_length = gene_isoforms_length_dict[rname][gname][isoform]
        for read in simulation_reads[isoform]:
            [read_start_offset,read_end_offset] = read
            read_start_offset += 1
            read_length = read_end_offset - read_start_offset + 1
            exons = []
            for exon in gene_exons_dict[rname][gname]:
                [exon_start,exon_end,exon_len] = exon
                region = 'P'+str(gene_points_dict[rname][gname][exon_start])+':'+'P'+str(gene_points_dict[rname][gname][exon_end])
                if region in gene_regions_dict[rname][gname]:
                    if isoform in gene_regions_dict[rname][gname][region]:
                        exons.append(exon)
            exons = sorted(exons,key=lambda x:(x[0],x[1]))
            if strand_dict[gname] == '+':
                curr_read_start_offset = read_start_offset
            else:
                curr_read_start_offset = isoform_length - (read_start_offset + read_length) + 1 + 1
            read_start_exon_index = -1
            read_start_pos = 1
            for i in range(len(exons)):
                exon_len = exons[i][2]
                exon_start =  exons[i][0]
                exon_end =  exons[i][1]
                if curr_read_start_offset - exon_len > 0:
                    curr_read_start_offset -= exon_len
                else:
                    read_start_pos = exon_start + curr_read_start_offset - 1
                    read_start_exon_index = i
                    break
            all_read_exons = []
            remaining_read_length = read_length
            for i in range(read_start_exon_index,len(exons)):
                exon_start = exons[i][0]
                exon_end = exons[i][1]
                exon_len = exons[i][2]
                if read_start_pos > exon_start:
                    this_exon_start_pos = read_start_pos
                else:
                    this_exon_start_pos = exon_start
                if this_exon_start_pos + remaining_read_length - 1 > exon_end:
                    remaining_read_length -= exon_end - this_exon_start_pos + 1
                    all_read_exons.append((this_exon_start_pos,exon_end))

                else:
                    read_end_pos = this_exon_start_pos + remaining_read_length - 1
                    all_read_exons.append((this_exon_start_pos,read_end_pos))
                    break
            region = ''
            if len(all_read_exons) == 1:
                last_coord = 0
                last_point = 0
                for coord,point in gene_points_dict[rname][gname].items():
                    if coord > all_read_exons[0][0] and last_coord <= all_read_exons[0][0]:
                        region = 'P' + str(last_point) + ':P' + str(last_point+1)
                    last_coord = coord
                    last_point = point
            else:
                for i in range(len(all_read_exons)):
                    if i == 0:
                        first_point = 'P' + str(gene_points_dict[rname][gname][all_read_exons[0][1]] - 1)
                        end_point = 'P' + str(gene_points_dict[rname][gname][all_read_exons[0][1]])
                        region = first_point + ':' + end_point
                    elif i == len(all_read_exons) - 1:
                        first_point = 'P' + str(gene_points_dict[rname][gname][all_read_exons[-1][0]])
                        end_point = 'P' +  str(gene_points_dict[rname][gname][all_read_exons[-1][0]] + 1)
                        region += '-' + first_point + ':' + end_point
                    else:
                        first_point = 'P' + str(gene_points_dict[rname][gname][all_read_exons[i][0]])
                        end_point = 'P' +  str(gene_points_dict[rname][gname][all_read_exons[i][1]])
                        region += '-' + first_point + ':' + end_point
            if region not in weight_dict[isoform]:
                weight_dict[isoform][region] = 0
            weight_dict[isoform][region] += 1
    return weight_dict
import glob
def worker(rname,gname,kde_model,num_reads_each_isoform):
    simulation_reads = simulated_reads_for_gene(kde_model,num_reads_each_isoform,gene_isoforms_length_dict,strand_dict,rname,gname)
    gene_weight_dict = assign_reads_to_region(simulation_reads,gene_isoforms_length_dict,gene_points_dict,gene_exons_dict,strand_dict,rname,gname)
    # with open(f'{output_dir}/{gname}','wb') as f:
    #     pickle.dump(gene_weight_dict,f)
    return gene_weight_dict
def submit_worker(list_of_args, cstart, cstop,id):
    kde_model_path = list_of_args[0][2]
    with open(kde_model_path,'rb') as f:
        kde_model = pickle.load(f)
    all_gene_weight_dict = {}
    for args in list_of_args[cstart:cstop]:
        rname,gname,_,num_reads_each_isoform = args
        gene_weight_dict = worker(rname,gname,kde_model,num_reads_each_isoform)
        if rname not in all_gene_weight_dict:
            all_gene_weight_dict[rname] = {}
        if gname not in all_gene_weight_dict[rname]:
            all_gene_weight_dict[rname][gname] = gene_weight_dict
    with open(f'{output_dir}/{id}','wb') as f:
        pickle.dump(all_gene_weight_dict,f)
list_of_args = []
num_reads_each_isoform = 1000
for kde_model_path in [kde_model_path]:
    for rname in gene_isoforms_dict:
        if rname != selected_rname:
            continue
        for gname in gene_isoforms_dict[rname]:
            list_of_args.append((rname,gname,kde_model_path,num_reads_each_isoform))
print(selected_rname)
print(len(list_of_args))
chunksize, extra = divmod(len(list_of_args), threads)
if extra:
    chunksize += 1
import concurrent
# all_gene_weight_dict = {}
with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    with tqdm(total=len(list_of_args)) as pbar:
        futures = []
        # for args in list_of_args:
        #     futures.append(executor.submit(worker,args))
        for i in range(threads):
            cstart = chunksize * i
            cstop = chunksize * (i + 1) if i != threads - 1 else len(list_of_args)
            futures.append(executor.submit(submit_worker,list_of_args, cstart, cstop,i))
        for fut in concurrent.futures.as_completed(futures):
            _ = fut.result()
            # rname,gname,_,_ = args 
            # if rname not in all_gene_weight_dict:
            #     all_gene_weight_dict[rname] = {}
            # if gname not in all_gene_weight_dict[rname]:
            #     all_gene_weight_dict[rname][gname] = future
            # for rname in future:
            #     if rname not in all_gene_weight_dict:
            #         all_gene_weight_dict[rname] = {}
            #     for gname in future[rname]:
            #         all_gene_weight_dict[rname][gname] = future[rname][gname]
            pbar.update(chunksize)
    # for i in range(threads):
    #     cstart = chunksize * i
    #     cstop = chunksize * (i + 1) if i != threads - 1 else len(list_of_args)
    #     futures.append(executor.submit(submit_worker,list_of_args, cstart, cstop))




    # results = list(tqdm(executor.map(worker,list_of_args,chunksize=chunksize),total=len(list_of_args)))
    # for args, gene_weight_dict in zip(list_of_args,results):
    #     rname,gname,_,_ = args
    #     if rname not in all_gene_weight_dict:
    #         all_gene_weight_dict[rname] = {}
    #     all_gene_weight_dict[rname][gname] = gene_weight_dict
# with open(output_dir,'wb') as f:
#     pickle.dump(all_gene_weight_dict,f)