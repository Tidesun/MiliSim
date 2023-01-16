import argparse
import config
from pathlib import Path
import numpy as np
import dill as pickle
import argparse
import config
import shutil
from build_kde_model import parse,map_long_reads,get_singleton_genes
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
                            singleton_genes_dist.append([isoform,isoform_len,start_offset,end_offset])
                        else:
                            singleton_genes_dist.append([isoform,isoform_len,end_offset,start_offset])
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
                                multi_isoform_genes_dist.append([isoform,isoform_len,start_offset,end_offset])
                            else:
                                multi_isoform_genes_dist.append([isoform,isoform_len,end_offset,start_offset])
    return multi_isoform_genes_dist
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
    # dist_arr [[isoform name, isoform length, 5' end truncation distance, 3' end truncation distance]]
    dist_arr = np.array(dist_lst)
    
    # with open(f'{output_path}/dist','wb') as f:
    #     pickle.dump([dist_arr,gene_regions_read_mapping,singleton_genes,isoform_len_dict,isoform_exon_dict,strand_dict,gene_isoform_dict,LR_gene_regions_dict],f)
    num_reads = dist_arr.shape[0]
    print(f'Number of reads to calculate the distribution: {num_reads}')
    with open(f'{output_path}/truncation_and_read_length_singleton_gene','wb') as f:
        pickle.dump(np.array(singleton_genes_dist),f)
    with open(f'{output_path}/truncation_and_read_length_multi_isoforms_gene_unique_region','wb') as f:
        pickle.dump(np.array(multi_isoform_genes_dist),f)
    with open(f'{output_path}/truncation_and_read_length','wb') as f:
        pickle.dump(dist_arr,f)
    try:
        shutil.rmtree(f'{output_path}/temp/')
    except:
        pass
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Calculate truncation and read length distribution of the data",add_help=True)
    
    parser.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    parser.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    parser.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    parser.add_argument('-t','--threads', type=int, default=1,help="Threads",required=True)
    args = parser.parse_args()
    config.output_path = args.output_path
    main(args.gtf_annotation_path,args.long_read_sam_path,args.output_path,threads=args.threads)
if __name__ == "__main__":
    parse_arguments()