
import pysam
import time
from pathlib import Path
from parse_annotation_main import parse_reference_annotation, process_annotation_for_alignment
from parse_alignment_main import parse_alignment
import re
import numpy as np
import dill as pickle
import argparse
import config
import shutil
def infer_read_len(short_read_alignment_file_path):
    with pysam.AlignmentFile(short_read_alignment_file_path, "r") as samfile:
            for read in samfile:
                READ_LEN = read.infer_query_length()
                if READ_LEN is not None:
                    break
    return READ_LEN
def map_long_reads(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,filtering,output_path,multi_mapping_filtering,threads,raw_isoform_exons_dict):
    print('Mapping long read to regions...',flush=True)
    start_time = time.time()
    if multi_mapping_filtering == 'unique_only':
        pysam.view('-F','2820','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        long_read_alignment_file_path = f'{output_path}/temp_lr.sam'
    elif multi_mapping_filtering == 'best':
        pysam.view('-F','2820','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        long_read_alignment_file_path = f'{output_path}/temp_lr.sam'
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,gene_regions_read_mapping,total_long_read_length,num_LRs,filtered_gene_regions_read_length = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,filtering,threads)
    try:
        Path(f'{output_path}/temp_sr.sam').unlink()
    except:
        pass
    try:
        Path(f'{output_path}/temp_lr.sam').unlink()
    except:
        pass
    print('Mapped {} long reads'.format(num_LRs,flush=True))
    # print('Constructing matrix and calculating condition number...',flush=True)
    # long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict, LR_gene_regions_dict, long_read_gene_regions_read_count, long_read_gene_regions_read_length,
    #                                                                    LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, num_LRs, total_long_read_length, READ_JUNC_MIN_MAP_LEN, output_path, threads)
    end_time = time.time()
    print('Done in %.3f s'%(end_time-start_time),flush=True)
    return gene_regions_read_mapping
def parse(ref_file_path, READ_LEN, READ_JUNC_MIN_MAP_LEN, short_read_alignment_file_path, threads):
    print('Start parsing annotation...')
    start_time = time.time()
    if short_read_alignment_file_path is not None:
        READ_LEN = infer_read_len(short_read_alignment_file_path)
    gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict = parse_reference_annotation(
        ref_file_path, threads, READ_LEN, READ_JUNC_MIN_MAP_LEN, 'read_length')
    _, gene_range, gene_interval_tree_dict = process_annotation_for_alignment(
        gene_exons_dict, gene_points_dict)
    end_time = time.time()
    print('Done in %.3f s' % (end_time-start_time), flush=True)
    return gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict, gene_range, gene_interval_tree_dict
def get_singleton_gene_dist(gene_regions_read_mapping,singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict,gene_isoform_dict):
    singleton_genes_dist = []
    for rname in gene_regions_read_mapping:
        for gname in gene_regions_read_mapping[rname]:
            if gname in singleton_genes:
                for region in gene_regions_read_mapping[rname][gname]:
                    for read_mapping in gene_regions_read_mapping[rname][gname][region]:
                        read_start,read_end = read_mapping['read_pos']
                        isoform = next(iter(gene_isoform_dict[gname]))
                        isoform_len = isoform_len_dict[isoform]
                        start_offset = 0
                        for [exon_start,exon_end] in isoform_exon_dict[isoform]:
                            if exon_end > read_start:
                                start_offset += read_start - exon_start + 1
                                break
                            else:
                                start_offset += exon_end - exon_start + 1
                        end_offset = isoform_len_dict[isoform] - start_offset - read_mapping['read_length']
                        if start_offset < 0:
                            start_offset = 0
                        if end_offset < 0:
                            end_offset = 0
                        if strand_dict[gname]  == '+':
                            singleton_genes_dist.append([isoform_len,start_offset,end_offset])
                        else:
                            singleton_genes_dist.append([isoform_len,end_offset,start_offset])
    return singleton_genes_dist
        
def get_multi_isoform_dist(gene_regions_read_mapping,LR_gene_regions_dict,singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict):
    multi_isoform_genes_dist = []
    for rname in gene_regions_read_mapping:
        for gname in gene_regions_read_mapping[rname]:
            if gname not in singleton_genes:
                for region in gene_regions_read_mapping[rname][gname]:
                    if len(gene_regions_read_mapping[rname][gname][region]) == 0:
                        continue
                    if len(LR_gene_regions_dict[rname][gname][region]) == 1:
                        isoform = next(iter(LR_gene_regions_dict[rname][gname][region]))
                        if '-' in isoform:
                            continue
                        isoform_len = isoform_len_dict[isoform]
                        for read_mapping in gene_regions_read_mapping[rname][gname][region]:
                            read_start,read_end = read_mapping['read_pos']
                            start_offset = 0
                            for [exon_start,exon_end] in isoform_exon_dict[isoform]:
                                if exon_end > read_start:
                                    start_offset += read_start - exon_start + 1
                                    break
                                else:
                                    start_offset += exon_end - exon_start + 1
                            end_offset = isoform_len_dict[isoform] - start_offset - read_mapping['read_length']
                            if start_offset < 0:
                                start_offset = 0
                            if end_offset < 0:
                                end_offset = 0
                            if strand_dict[gname]  == '+':
                                multi_isoform_genes_dist.append([isoform_len,start_offset,end_offset])
                            else:
                                multi_isoform_genes_dist.append([isoform_len,end_offset,start_offset])
    return multi_isoform_genes_dist
def get_singleton_genes(annotation):
    gene_exon_dict = {}
    gene_isoform_dict = {}
    isoform_exon_dict = {}
    strand_dict = {}
    with open(annotation,'r') as f:
        for line in f:
            if line.lstrip()[0] == "#":
                continue
            fields = line.split('\t')
            if (fields[2] != 'exon'):
                continue
            strand = fields[6]
            chr_name = fields[0]
            gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
            isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            if gene_name not in gene_exon_dict:
                gene_exon_dict[gene_name] = []
                gene_isoform_dict[gene_name] = set()
            if isoform_name not in isoform_exon_dict:
                isoform_exon_dict[isoform_name] = []
            gene_exon_dict[gene_name].append([start_pos,end_pos])
            gene_isoform_dict[gene_name].add(isoform_name)
            isoform_exon_dict[isoform_name].append([start_pos,end_pos])
            strand_dict[gene_name] = strand
    for isoform in isoform_exon_dict:
        isoform_exon_dict[isoform] = sorted(isoform_exon_dict[isoform],key=lambda x:(x[0],x[1]))
    isoform_len_dict = {}
    for isoform in isoform_exon_dict:
        isoform_len_dict[isoform] = 0
        for exon in isoform_exon_dict[isoform]:
            isoform_len_dict[isoform] += exon[1] - exon[0] + 1
    singleton_genes = set()
    for gname in gene_isoform_dict:
        if len(gene_isoform_dict[gname]) == 1:
            singleton_genes.add(gname)
    return singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict,gene_isoform_dict
from sklearn.neighbors import KernelDensity
def kde3d(x, y,z):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    xyz = np.vstack([x, y,z])
    d = xyz.shape[0]
    n = xyz.shape[1]
    bw = (n * (d + 2) / 4.) ** (-1. / (d + 4))  # silverman
    kde_3d = KernelDensity(bandwidth=bw).fit(xyz.T)

    return kde_3d
'''
We use reads mapped to singleton gene(gene with single isoform) [get_singleton_gene_dist] 
and reads mapped to the unique region of isoform (the region not shared by any other isoform within the sam gene)[get_multi_isoform_dist]
'''
def main(ref_file_path,long_read_alignment_file_path,output_path,filtering='False',multi_mapping_filtering='best',threads=1,READ_LEN=0,READ_JUNC_MIN_MAP_LEN=15):
    Path(output_path).mkdir(exist_ok=True,parents=True)
    try:
        shutil.rmtree(f'{output_path}/temp/')
    except:
        pass
    _,gene_points_dict,gene_isoforms_dict,\
        _,_,LR_gene_regions_dict,LR_genes_regions_len_dict,\
            gene_isoforms_length_dict,raw_isoform_exons_dict,_,\
                _,_,gene_range,gene_interval_tree_dict = \
                    parse(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,None,threads)
    gene_regions_read_mapping = map_long_reads(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,filtering,output_path,multi_mapping_filtering,threads,raw_isoform_exons_dict)
    singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict,gene_isoform_dict = get_singleton_genes(ref_file_path)
    singleton_genes_dist = get_singleton_gene_dist(gene_regions_read_mapping,singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict,gene_isoform_dict)
    multi_isoform_genes_dist = get_multi_isoform_dist(gene_regions_read_mapping,LR_gene_regions_dict,singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict)
    dist_lst = multi_isoform_genes_dist+singleton_genes_dist
    dist_arr = np.array(dist_lst)
    # with open(f'{output_path}/dist','wb') as f:
    #     pickle.dump([dist_arr,gene_regions_read_mapping,singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict,gene_isoform_dict,LR_gene_regions_dict],f)
    num_reads = dist_arr.shape[0]
    print(f'Number of reads to build KDE model: {num_reads}')
    # x: isoform length y: 5' end truncation distance z: 3' end truncation distance
    x,y,z = dist_arr[:,0],dist_arr[:,1],dist_arr[:,2]
    kde = kde3d(x,y,z)
    with open(f'{output_path}/kde_model','wb') as f:
        pickle.dump(kde,f)
    try:
        shutil.rmtree(f'{output_path}/temp/')
    except:
        pass
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Build kde model on real data",add_help=True)
    
    parser.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    parser.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    parser.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    parser.add_argument('-t','--threads', type=int, default=1,help="Threads",required=True)
    args = parser.parse_args()
    config.output_path = args.output_path
    main(args.gtf_annotation_path,args.long_read_sam_path,args.output_path,threads=args.threads)
if __name__ == "__main__":
    parse_arguments()