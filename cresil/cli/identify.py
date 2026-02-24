#!/usr/bin/env python
import os, sys, signal
import math
import datetime
import pathlib
import subprocess
import operator
import pandas as pd
import numpy as np
import pysam, shutil
import pybedtools as bt
import graphviz as gv
import networkx as nx
import matplotlib.pyplot as plt
from multiprocessing import cpu_count, Pool
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Platform-specific configurations
PLATFORM_CONFIG = {
    "ont": {
        "use_medaka": True,
    },
    "hifi": {
        "use_medaka": False,  # HiFi has high accuracy, skip sequence correction
    }
}

#####################################################################
## utility fuctions #################################################

def lenLoci(loci):
    length = 0
    for region in loci.split(','):
        chrom, start, end, strand = region.rsplit('_', 3)
        length += int(end)-int(start)
    return length

def format_merge_region(merge_region):
    list_formatted = []
    for region in merge_region.split(','):
        chrom, start, end, strand = region.rsplit('_', 3)
        list_formatted.append("{}:{}-{}_{}".format(chrom, int(start) + 1, end, strand))
    return ','.join(list_formatted)

def draw_graph(graph, filename="test"):
    import graphviz as gv
    d = gv.Digraph(name=filename)
    for pair, weight in graph.items():
        d.edge(pair[0], pair[1], label=str(weight))
    d.render()
    return filename

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

def majority_strand(df_merged):
    strand = '+'
    list_str_strand = df_merged['strand'].astype(str).tolist()
    if list_str_strand.count('-1') > list_str_strand.count('1'):
        strand = '-'
    return strand

def chk_ovl_nodes(df_final_read_mergeid, list_nodes):
    list_uniq_nodes = list(set(list_nodes))
    l1 = df_final_read_mergeid['mergeid'].isin(list_uniq_nodes)
    df_x = df_final_read_mergeid[l1].groupby(by='readid')[['mergeid']].nunique().sort_values(by='mergeid', ascending=False)
    list_readid = list(set(df_x[df_x['mergeid'] == len(list_uniq_nodes)].index))

    chk_return = False
    if len(list_readid) > 0:
        df_y = df_final_read_mergeid[df_final_read_mergeid['readid'].isin(list_readid)].groupby(by='mergeid').agg({'ovl_5end':['sum'], 'ovl_3end':['sum']})
        list_boo = []
        for node in list_uniq_nodes:
            list_boo.append(df_y.loc[node,'ovl_3end']['sum'] > 0)
            list_boo.append(df_y.loc[node,'ovl_5end']['sum'] > 0)
        if all(list_boo):
            chk_return = True

    return chk_return

def chk_circular_subgraph(G, graph, dict_pair_strand):
    list_nodes = list(graph.nodes)
    num_nodes = len(list_nodes)
    print_nodes = ",".join(list_nodes)
    sl = True if len(list(nx.selfloop_edges(graph))) > 0 else False
    solved = True
    cyclic = True
    if num_nodes == 2:
        rm_sl_subgraph = graph.copy()
        rm_sl_subgraph.remove_edges_from(nx.selfloop_edges(rm_sl_subgraph))
        solved = all(x <= 2 for x in [rm_sl_subgraph.degree[node] for node in list_nodes])
        key_forward = tuple(list_nodes)
        key_reverse = tuple(list_nodes[::-1])
        cyclic = False
        if key_forward in dict_pair_strand and key_reverse in dict_pair_strand:
            for fstrands in dict_pair_strand[key_forward]:
                fl, fr = fstrands.split('_')
                for rstrands in dict_pair_strand[key_reverse]:
                    rl, rr = rstrands.split('_')
                    if fr == rl and fl == rr:
                        cyclic = True
                        break
    elif num_nodes > 2:
        rm_sl_subgraph = graph.copy()
        rm_sl_subgraph.remove_edges_from(nx.selfloop_edges(rm_sl_subgraph))
        solved = all(x <= 2 for x in [rm_sl_subgraph.degree[node] for node in list_nodes])
        cyclic = True if len(nx.cycle_basis(nx.DiGraph(rm_sl_subgraph).to_undirected())) > 0 else False
    else:
        cyclic = True if len(nx.cycle_basis(nx.DiGraph(graph).to_undirected())) > 0 else False

    if not solved:
        cyclic = False

    test_graph = graph.copy()
    test_graph.remove_edges_from(nx.selfloop_edges(test_graph))
    list_traversal = nx.cycle_basis(nx.DiGraph(test_graph).to_undirected())
    if len(list_traversal) > 0:
        list_traversal_ini = list_traversal[0]
        if len(list_traversal_ini) == len(list_nodes):
            print_nodes = ",".join(list_traversal_ini)

    return (print_nodes, num_nodes, solved, sl, cyclic)

def reverse_strand(strand):
    return ''.join(['+' if x == '-' else '-' if x == '+' else x for x in list(strand)])

def get_longest_node(chk_subgraph):
    list_nodes = list(chk_subgraph.nodes)
    longest_node = list_nodes[0]
    for node in list_nodes[1:]:
        if lenLoci(node) > lenLoci(longest_node):
            longest_node = node
    return longest_node

def get_dict_weight_subgraph(directed_graph, chk_subgraph):
    dict_return = {}
    for node_name in list(chk_subgraph.nodes):
        for k, v in directed_graph.adj[node_name].items():
            key = "{},{}".format(node_name, k)
            dict_return[key] = v[0]['weight']
    return dict_return

def nodes_to_merge_regions(list_nodes):
    list_merge_regions = []
    for node in list_nodes:
        chr_, start_, end_ = node.rsplit('_', 2)
        region_ = "{}:{}-{}".format(chr_, int(start_) + 1, end_)
        list_merge_regions.append(region_)
    return ','.join(list_merge_regions)

def get_fasta_length(fasta_path):
    total_len = 0
    for record in SeqIO.parse(fasta_path, 'fasta'):
        total_len += len(str(record.seq))
    return total_len

def get_consensus_status(value, dict_medaka_seq_len):
    status = '0|Unknown'
    id_ = value['id']
    merge_len = value['merge_len']
    consensus_len = dict_medaka_seq_len.get(id_, None)
    if consensus_len:
        if merge_len > consensus_len:
            status = "{}|Deletion".format(consensus_len)
        elif merge_len < consensus_len:
            status = "{}|Insertion".format(consensus_len)
        else:
            status = "{}|No_indel".format(consensus_len)
    return status

def any_overlapping_range(start1, end1, start2, end2):
    if start1 <= end2 and start2 <= end1:
        return True
    else:
        return False

def check_breakpoint_direction(df_check):

    list_pairs = []

    df_check.reset_index(drop=True, inplace=True)

    list_order_check = [(i , i+1) for i in df_check.index.tolist() if i + 1 < len(df_check)]

    for tup_check in list_order_check:

        l1 = df_check.index.isin(tup_check)

        q_starts = df_check[l1]['q_start'].tolist()
        q_ends = df_check[l1]['q_end'].tolist()

        if any_overlapping_range(q_starts[0], q_ends[0] + 50, q_starts[1], q_ends[1]):

            list_mergeid = df_check[l1]['mergeid'].tolist()
            list_check = df_check[l1]['strand'].tolist()

            if list_check[0] == 1 and list_check[1] == 1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_3[0] == 1 and list_5[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], '+_+' , True)
                    list_pairs.append(pair)

            if list_check[0] == -1 and list_check[1] == -1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_3[1] == 1 and list_5[0] == 1:
                    pair = (list_mergeid[1], list_mergeid[0], '-_-', True)
                    list_pairs.append(pair)


            if list_check[0] == -1 and list_check[1] == 1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_5[0] == 1 and list_5[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], '-_+', True)
                    list_pairs.append(pair)


            if list_check[0] == 1 and list_check[1] == -1:
                list_5 = df_check[l1]['ovl_5end'].tolist()
                list_3 = df_check[l1]['ovl_3end'].tolist()
                if list_3[0] == 1 and list_3[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], '+_-', True)
                    list_pairs.append(pair)
        else:
            break

    return list_pairs

def check_abs_ovl(value):

    cut_off=1

    r_start = value['r_start']
    r_end = value['r_end']
    mergeid = value['mergeid']
    merge_5end = int(mergeid.split('_')[1])
    merge_3end = int(mergeid.split('_')[2])

    ovl_5end_value = abs(r_start - merge_5end)
    ovl_3end_value = abs(r_end - merge_3end)

    ovl_5end = 1 if ovl_5end_value <= cut_off else 0
    ovl_3end = 1 if ovl_3end_value <= cut_off else 0
    sum_ends = ovl_5end + ovl_3end

    return ovl_5end, ovl_3end, sum_ends

def cal_skip_variant(gname, merge_region, num_region, is_cyclic, read_merged_ins_df):

    lengthLoci = lenLoci(merge_region)
    ctc = 'True' if is_cyclic == True else 'False'

    ## collect region lengths of overlapping reads using vectorized ops
    nodes = [node[:-2] for node in merge_region.split(',')]
    mask = read_merged_ins_df['mergeid'].isin(nodes)
    filtered = read_merged_ins_df.loc[mask]
    total_base = int((filtered['q_end'] - filtered['q_start']).sum())
    num_reads = filtered['readid'].astype(str).nunique()

    expectCov = "{:.2f}".format((float(total_base) / lengthLoci))

    return (gname, merge_region, lengthLoci, num_region, ctc, num_reads, total_base, expectCov)

def prepare_identify_seq(gname, merge_region, num_region, is_cyclic, assemGraph, fastaRef, fastaName, read_merged_ins_df):

    fa_ref = pysam.Fastafile(fastaRef)
    fa_reads = pysam.FastaFile(fastaName)

    lengthLoci = lenLoci(merge_region)
    ctc = 'True' if is_cyclic == True else 'False'

    ## create a subgraph folder
    assemFol = "{}/{}".format(assemGraph, gname)
    pathlib.Path(assemFol).mkdir(parents=True, exist_ok=True)

    ## write a reference sequence from regions
    ref_fasta_path = "{}/reference_regions.fa".format(assemFol)
    seq_ = ""
    with open(ref_fasta_path, 'w') as w_f:
        w_f.write(">{}\n".format(gname))
        for node in merge_region.split(','):
            chr_, start_, end_, strand_ = node.rsplit('_', 3)
            chr_ = str(chr_)
            if strand_ == '-':
                seq_ += str(Seq(fa_ref.fetch(chr_, int(start_), int(end_))).reverse_complement()).upper()
            else:
                seq_ += str(fa_ref.fetch(chr_, int(start_), int(end_))).upper()
        w_f.write(seq_ + "\n")

    ## write all regions in reads associating eccDNA regions
    total_base = 0
    list_readid = []

    nodes = [node[:-2] for node in merge_region.split(',')]
    mask = read_merged_ins_df['mergeid'].isin(nodes)
    filtered = read_merged_ins_df.loc[mask]

    out_filename_final = "{}/final_reads.fa".format(assemFol)
    with open(out_filename_final, mode='w') as foutF:
        for row in filtered.itertuples():
            readO_readid = str(row.readid)
            readO_start = row.q_start
            readO_end = row.q_end
            readO_name = "{}_{}_{}_{}".format(gname, readO_readid, readO_start, readO_end)
            readO_seq = str(fa_reads.fetch(readO_readid, readO_start, readO_end)).upper()
            foutF.write(">{}\n{}\n".format(readO_name, readO_seq))

            total_base += len(readO_seq)
            list_readid.append(readO_readid)

    expectCov = "{:.2f}".format((float(total_base) / lengthLoci))

    return (gname, merge_region, lengthLoci, num_region, ctc, len(set(list_readid)), total_base, expectCov)

def get_a_sequence(fa_path):
    seq = ''
    for record in SeqIO.parse(fa_path, "fasta"):
        seq = str(record.seq)
    return seq

def assemToGFA(tup_value):
    fa_path, tup_ = tup_value
    ec_id, assembly_len, coverage, ctc = tup_
    seq = get_a_sequence(fa_path)

    gfaS = []
    gfaL = []
    if ctc == "True":
        gfaS.append(f"S\t{ec_id}\t{seq}\tLN:i:{assembly_len}\tdp:f:{coverage}")
        gfaL.append(f"L\t{ec_id}\t+\t{ec_id}\t+\t0M")
    else:
        gfaS.append(f"S\t{ec_id}\t{seq}\tLN:i:{assembly_len}\tdp:f:{coverage}")

    write_path = "{}/{}.gfa".format(str(pathlib.Path(fa_path).parent), str(pathlib.Path(fa_path).stem))
    with open(write_path, 'w') as outf:
        for s in gfaS:
            outf.write(f"{s}\n")
        for l in gfaL:
            outf.write(f"{l}\n")

def dist_work(cmd):
    process = subprocess.call("{}".format(cmd), shell=True)

def argparser():

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Identify and verify eccDNA from enriched data')
    general = parser.add_argument_group(title='General options')
    general.add_argument('-platform', "--platform", dest='platform',
                            help="sequencing platform: ont or hifi [ont]",
                            type=str, choices=['ont', 'hifi'], default='ont')
    general.add_argument('-t', "--threads",
                            help="Number of threads [all CPU cores]",
                            type=int, default=0)
    general.add_argument('-s', "--skip-variant", dest='skipvariant',
                            help="skip sequence correction and variant calling steps [False]",
                            action='store_true')
    general.add_argument('-sg', "--skip-gfa", dest='skipgfa',
                            help="skip creating GFA files for each eccDNA (effective without -s) [False]",
                            action='store_false')
    general.add_argument('-fa', "--fa-ref", dest='faref',
                            help="reference sequence .fasta",
                            type=str, default=None)
    general.add_argument('-fai', "--fa-index", dest='fai',
                            help="genome size file e.g. hg19.chrom.sizes or file.fasta.fai",
                            type=str, default=None)
    general.add_argument('-fq', "--fq-input",  dest='fqinput',
                            help="input fasta/fastq",
                            type=str, default=None)
    general.add_argument('-trim', "--trim-input", dest='triminput',
                            help="CReSIL trim table",
                            type=str, default=None)
    general.add_argument('-minrsize', "--minimum-region-size", dest='minrsize',
                            help="minimum size of region of eccDNA [200]",
                            type=int, default=200)
    general.add_argument('-ovl', "--ovl-size", dest='ovlsize',
                            help="size of 5' and 3' regions on reads for a breakpiont overlapping check [50]",
                            type=int, default=50)
    general.add_argument('-break', "--brakpoint-depth", dest='breakpointdepth',
                            help="lowest number of supported breakpoints of eccDNA [3]",
                            type=int, default=3)
    general.add_argument('-depth', "--average-depth", dest='averagedepth',
                            help="lowest average depth of eccDNA [5]",
                            type=int, default=5)
    general.add_argument('-cm', "--consensus-model", dest='consensus_model',
                            help="Medaka_consensus model (only for ONT) [r1041_e82_400bps_sup_v4.3.0]",
                            type=str, default='r1041_e82_400bps_sup_v4.3.0')
    return parser

def main(args):

    # Get platform configuration
    platform = args.platform
    config = PLATFORM_CONFIG[platform]

    threads = cpu_count() if args.threads == 0 else args.threads

    # For HiFi platform, automatically skip Medaka (sequence correction and variant calling)
    if platform == "hifi":
        skipvariant = True  # HiFi has high accuracy, no need for Medaka
    else:
        skipvariant = args.skipvariant

    skipgfa = args.skipgfa
    fastaRef = args.faref
    chromSizes = args.fai
    fastaName = args.fqinput
    minrsize = args.minrsize
    check_ovl_size = args.ovlsize
    breakpointdepth = args.breakpointdepth
    regiondepth = args.averagedepth
    consensus_model = args.consensus_model

    if not args.triminput:
        sys.stderr.write("[ABORT] a table from trimming step (trim.txt) is needed\n")
        exit(1)

    if args.threads == 0:
        threads = os.cpu_count()

    bname = str(pathlib.Path(args.triminput).parent)

    readTrim = pd.read_csv(args.triminput, sep="\t", header=0)

    outDir="{}/cresil_run".format(bname)
    tmpDir="{}/tmp".format(outDir)
    assemGraph = "{}/assemGraph".format(outDir)

    ## input ############################################################

    #####################################################################
    ## output ###########################################################

    pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(tmpDir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(assemGraph).mkdir(parents=True, exist_ok=True)

    ## output ###########################################################

    ct = datetime.datetime.now()
    print("\n######### CReSIL : Start identifying process (platform: {}, thread: {})\n[{}] total trimmed region : {}".format(platform.upper(), threads, ct, len(readTrim)), flush=True)
    if platform == "hifi":
        print("[{}] HiFi mode: skipping Medaka sequence correction (high accuracy data)".format(ct), flush=True)
    if len(readTrim) == 0:
        sys.exit("[ABORT] zero trimmed region\n")

    #####################################################################
    ## Identify potential eccDNA regions ################################

    ## Calculate read coverage and filter out restuion with depth <= 5x
    ord_header = ['ref', 'r_start', 'r_end', 'readid', 'q_start', 'q_end', 'match',
                  'mapBlock', 'mapq', 'strand', 'qlenTrimmed', 'freqCov', 'order']
    readTrim = readTrim.loc[:,ord_header]

    ## prepare aligned reads
    aln_reads = bt.BedTool.from_dataframe(readTrim).sort()

    ## calculate reads coverage by bedtools genomecov
    genome_cov = aln_reads.genome_coverage(bg=True, g=chromSizes)

    ## extract merged regions based on aligned reads
    aln_reads_merge = aln_reads.merge()

    ## join merged regions with read coverage and perform filtering
    merge_genomecov = aln_reads_merge.intersect(genome_cov, output='{}/merge_genomecov.txt'.format(tmpDir), wo=True)
    merge_genomecov_df = pd.read_csv("{}/merge_genomecov.txt".format(tmpDir),sep="\t",header=None,dtype=str)
    merge_genomecov_df.columns = ['m_chrom','m_start','m_end','bg_chrom','bg_start','bg_end','depth','d_ovl']
    merge_genomecov_df[['depth','bg_start','d_ovl']] = merge_genomecov_df[['depth','bg_start','d_ovl']].apply(pd.to_numeric, errors='coerce')

    ## generate mergeID
    merge_genomecov_df['mergeid'] = merge_genomecov_df['m_chrom'] + "_" + merge_genomecov_df['m_start'] + "_" + merge_genomecov_df['m_end']

    ## trim low-coverage regions
    merge_genomecov_df_filt = merge_genomecov_df[merge_genomecov_df['depth'] >= regiondepth]
    merge_genomecov_df = merge_genomecov_df.sort_values(['mergeid','bg_start'])
    ## group reads coverage of each mergeid and calculate mean coverage
    merge_bg_filt = merge_genomecov_df_filt.groupby(['mergeid']).agg(
        {'bg_chrom': 'max', 'bg_start': 'min', 'bg_end': 'max', 'depth': 'mean'}).reset_index()
    ## create final potential eccDNA regions and
    ## re-create mergeid based on trimming low-coverage regions
    merge_bg_filt[['bg_start','bg_end']] = merge_bg_filt[['bg_start','bg_end']].apply(pd.to_numeric, errors='coerce')
    merge_bg_filt['length'] = merge_bg_filt.bg_end - merge_bg_filt.bg_start
    merge_bg_filt["bg_start"]= merge_bg_filt["bg_start"].astype(str)
    merge_bg_filt["bg_end"]= merge_bg_filt["bg_end"].astype(str)
    merge_bg_filt['mergeid'] = merge_bg_filt.bg_chrom+"_"+merge_bg_filt.bg_start+"_"+merge_bg_filt.bg_end
    merge_bg_filt = merge_bg_filt.loc[:,['bg_chrom','bg_start','bg_end','depth','length','mergeid']]

    # filter regions with the length >= 200 bp
    merge_bg_filt = merge_bg_filt[merge_bg_filt['length'] >= minrsize]

    ## Identify potential eccDNA regions ################################

    #####################################################################
    ## Assign readid to each potential eccDNA regions (mergeid) #########

    ct = datetime.datetime.now()
    print("[{}] calculating breakpoints and merging regions".format(ct), flush=True)

    ## create merge region
    aln_reads_merge = bt.BedTool.from_dataframe(merge_bg_filt)
    read_merged_intersect = aln_reads_merge.intersect(aln_reads, output='{}/reads_merge_intersect.bed'.format(tmpDir), wo=True)
    ## read merge region and create mergeid

    if os.stat("{}/reads_merge_intersect.bed".format(tmpDir)).st_size == 0:
        sys.exit("[ABORT] no merge region was detected (lower -breakpoint, -depth; or region too short)\n")


    read_merged_ins_df_org = pd.read_csv("{}/reads_merge_intersect.bed".format(tmpDir),sep="\t",header=None,dtype=str)
    header = list(merge_bg_filt.columns)+list(readTrim.columns)+["ovl"]
    read_merged_ins_df_org.columns = header
    header_select = ['ref','r_start','r_end','mergeid','depth','length','readid',
                     'q_start','q_end','match','mapBlock','mapq','strand','qlenTrimmed',
                     'freqCov','order','ovl']
    read_merged_ins_df_org = read_merged_ins_df_org[header_select]
    ## Assign readid to each potential eccDNA regions (mergeid) #########

    #####################################################################
    ## Annotate 5'end and 3'end region of each readid ###################

    ## collect 200bp from 5'end and 3'end of each potential eccDNA
    end_size = check_ovl_size if minrsize >= 200 else int(round(minrsize * 0.3, 0))
    end_size = max(end_size, 15)
    with open("{}/end5_merge_region.bed".format(tmpDir), 'w') as end5, \
         open("{}/end3_merge_region.bed".format(tmpDir), 'w') as end3:
        for row in merge_bg_filt.itertuples():
            chrom = row.bg_chrom
            start = row.bg_start
            end = row.bg_end
            end5.write("{}\t{}\t{}\n".format(chrom, start, int(start) + end_size))
            end3.write("{}\t{}\t{}\n".format(chrom, int(end) - end_size, end))
    bt_5end = bt.BedTool("{}/end5_merge_region.bed".format(tmpDir))
    bt_3end = bt.BedTool("{}/end3_merge_region.bed".format(tmpDir))
    ## Annotate each read that has 5'end and 3'end overlap
    ## add column ovl_5end and ovl_3end (0=no overlap, 1=presence of overlap)
    read_merged_ins_df_org_bt = bt.BedTool.from_dataframe(read_merged_ins_df_org)
    read_merged_final_bt = read_merged_ins_df_org_bt.annotate(files=["{}/end5_merge_region.bed".format(tmpDir),"{}/end3_merge_region.bed".format(tmpDir)], counts=True)

    read_merged_ins_df = read_merged_final_bt.to_dataframe(names=header_select + ['ovl_5end','ovl_3end'])
    ## sum 5'end and 3'end overlap (0=no overlap, 1=either 5'end or 5'end overlap, 2=both ends overlap)
    read_merged_ins_df['sum_ends'] = read_merged_ins_df.ovl_5end + read_merged_ins_df.ovl_3end

    ## create dictionary of nodes and their temporary strands
    dict_majority_strand = {}
    for group_name, group in read_merged_ins_df.groupby(by='mergeid'):
        dict_majority_strand[group_name] = majority_strand(group)

    ## create dictionary of node pairs and all graphs
    ct = datetime.datetime.now()
    print("[{}] analyzing graphs".format(ct), flush=True)

    dict_pair_strand = {}
    graph = {}

    for readid, group in read_merged_ins_df.groupby(by='readid'):

        df_check = group.sort_values(by='order')
        df_check.reset_index(drop=True, inplace=True)

        if len(df_check) > 1:
            for tup_result in check_breakpoint_direction(df_check):
                if tup_result[3] == True:

                    pair = (tup_result[0], tup_result[1])

                    if pair in dict_pair_strand:
                        dict_pair_strand[pair].add(tup_result[2])
                    else:
                        dict_pair_strand[pair] = {tup_result[2]}

                    graph[pair] = graph.get(pair, 0) + 1

    ## filter graphs with low breakpoints
    graphFilt = {k: v for k, v in graph.items() if v >= breakpointdepth}

    ## create all graph objects
    G = nx.MultiDiGraph()
    for pair, weight in graphFilt.items():
        G.add_edge(pair[0], pair[1], weight=weight)

    ## create list of subgraphs
    subgraphs = list(connected_component_subgraphs(G.to_undirected()))

    ct = datetime.datetime.now()
    print("[{}] initial subgraphs : {}".format(ct, len(subgraphs)), flush=True)

    ## correct eccDNA strands and create a graph summary file
    list_graph_summary = []

    for idx, graph in enumerate(subgraphs, start=0):
        nodes = list(graph.nodes())
        gname = "ec{}".format(idx + 1)

        regions, num_nodes, can_be_solved, contain_selfloop, is_cyclic = chk_circular_subgraph(G, graph, dict_pair_strand)

        if len(nodes) == 1:
            pair = (nodes[0], nodes[0])
            regions = "{}_{}".format(pair[0], list(dict_pair_strand[pair])[0][0])

        elif len(nodes) == 2:
            pair = tuple(nodes)
            list_region = list(pair)

            if pair not in dict_pair_strand:
                pair = tuple(nodes[::-1])

            list_strand = list(dict_pair_strand[pair])[0].split('_')
            regions = ','.join("{}_{}".format(x[0], x[1]) for x in zip(list_region, list_strand))

        else:
            test_graph = graph.copy()
            test_graph.remove_edges_from(nx.selfloop_edges(test_graph))

            list_traversal = nx.cycle_basis(nx.DiGraph(test_graph).to_undirected())
            if len(list_traversal) == 0:
                regions = ','.join(["{}_{}".format(node, dict_majority_strand[node]) for node in nodes])

            else:
                list_traversal = list_traversal[0]

                list_order_check = [(i, i+1) for i in range(len(list_traversal)) if i + 1 < len(list_traversal)]

                list_temp_all = []
                for tup_check in list_order_check:
                    l_region = list_traversal[tup_check[0]]
                    r_region = list_traversal[tup_check[1]]
                    pair = (l_region, r_region)

                    list_this_level = []

                    if pair in dict_pair_strand:
                        for str_strand in dict_pair_strand[pair]:
                            list_str_strand = ['_'.join(x) for x in zip([l_region, r_region], str_strand.split('_'))]
                            list_this_level.append(list_str_strand)
                    else:
                        rev_pair = (r_region, l_region)
                        for str_strand in [reverse_strand(x) for x in dict_pair_strand[rev_pair]]:
                            list_str_strand = ['_'.join(x) for x in zip([l_region, r_region], str_strand.split('_'))]
                            list_this_level.append(list_str_strand)

                    list_temp_all.append(list_this_level)

                list_traverse = list_temp_all[0]
                for list_level in list_temp_all[1:]:
                    tail_to_idx = {tl[-1]: i for i, tl in enumerate(list_traverse)}
                    for list_ in list_level:
                        idx = tail_to_idx.get(list_[0])
                        if idx is not None:
                            list_traverse[idx].append(list_[1])

                regions = ','.join(max(list_traverse, key=len))

        tup_graph_summary = (gname, regions, num_nodes, can_be_solved, contain_selfloop, is_cyclic)

        list_graph_summary.append(tup_graph_summary)

    if len(list_graph_summary) > 0:
        header = ['id', 'regions', 'num_nodes', 'can_be_solved', 'contain_selfloop', 'is_cyclic']
        df_graph_summary = pd.DataFrame(list_graph_summary, columns=header)
    else:
        sys.exit("[ABORT] no eccDNA was detected (lower -breakpoint, -depth; or region too short)\n")

    ## if skip sequence correction and variant calling steps
    if skipvariant:

        ct = datetime.datetime.now()
        print("[{}] skipped sequence correction and variant calling steps\n[{}] collecting data from subgraphs : {}".format(ct, ct, len(subgraphs)), flush=True)

        list_identify_results = []
        total_subgraphs = len(df_graph_summary)

        for counter, row in enumerate(df_graph_summary.itertuples(), start=1):

            if counter % 100 == 0 or counter == total_subgraphs:
                ct = datetime.datetime.now()
                print("[{}] {}/{}".format(ct, counter, total_subgraphs), flush=True)

            result = cal_skip_variant(
                row.id, row.regions, row.num_nodes, row.is_cyclic,
                read_merged_ins_df)
            list_identify_results.append(result)

        ## write a subGraph summary file
        header = ['id', 'merge_region', 'merge_len', 'num_region', 'ctc', 'numreads', 'totalbase', 'coverage']
        selected_header = ['id', 'merge_region', 'merge_len', 'num_region', 'can_be_solved', 'contain_selfloop', 'ctc', 'numreads', 'totalbase', 'coverage']
        df_temp_graph_summary = pd.DataFrame(list_identify_results, columns=header)
        df_final_graph_summary = pd.merge(left=df_graph_summary, right=df_temp_graph_summary, left_on='id', right_on='id', how='inner')
        df_final_graph_summary = df_final_graph_summary[selected_header]
        df_final_graph_summary['merge_region'] = df_final_graph_summary['merge_region'].apply(lambda x: format_merge_region(x))

        write_path = "{}/subGraphs.summary.txt".format(outDir)
        df_final_graph_summary.to_csv(write_path, sep='\t', index=None)

        ## write identified eccDNA summary
        header = ['id', 'merge_region', 'merge_len', 'num_region', 'ctc', 'numreads', 'totalbase', 'coverage']
        df_identify_summary = df_final_graph_summary[header].copy()
        eccdna_final_path = "{}/eccDNA_final.txt".format(bname)
        df_identify_summary.to_csv(eccdna_final_path, sep='\t', index=None)

        ct = datetime.datetime.now()
        print("[{}] finished identifying process\n".format(ct), flush=True)
        sys.exit()

    ## Run consensus sequence and variants
    ct = datetime.datetime.now()
    print("[{}] preparing data for correcting sequences from subgraphs : {}".format(ct, len(subgraphs)), flush=True)

    list_identify_results = []
    total_subgraphs = len(df_graph_summary)

    for counter, row in enumerate(df_graph_summary.itertuples(), start=1):

        if counter % 100 == 0 or counter == total_subgraphs:
            ct = datetime.datetime.now()
            print("[{}] {}/{}".format(ct, counter, total_subgraphs), flush=True)

        result = prepare_identify_seq(
            row.id, row.regions, row.num_nodes, row.is_cyclic,
            assemGraph, fastaRef, fastaName, read_merged_ins_df)
        list_identify_results.append(result)

    ## add parameters for Medaka
    list_params = []
    list_chk_consensus = []
    list_chk_variant = []

    for row in df_graph_summary.itertuples():
        gname = row.id
        assemFol = "{}/{}".format(assemGraph, gname)
        out_filename_final = "{}/final_reads.fa".format(assemFol)
        ref_fasta_path = "{}/reference_regions.fa".format(assemFol)
        consensusFol = "{}/medaka_consensus".format(assemFol)
        consensus_hdf_path = "{}/consensus_probs.hdf".format(consensusFol)
        variant_path = "{}/variant.vcf".format(consensusFol)
        tup_params = (out_filename_final, ref_fasta_path, consensusFol, consensus_model, consensus_hdf_path, variant_path)
        list_params.append(tup_params)

        consensus_path = "{}/consensus.fasta".format(consensusFol)
        dest_path = "{}/{}_consensus.fa".format(assemFol, gname)
        list_chk_consensus.append((consensus_path, dest_path))

        dest_path = "{}/{}_variant.vcf".format(assemFol, gname)
        list_chk_variant.append((variant_path, dest_path))

    ## write a subGraph summary file
    header = ['id', 'merge_region', 'merge_len', 'num_region', 'ctc', 'numreads', 'totalbase', 'coverage']
    selected_header = ['id', 'merge_region', 'merge_len', 'num_region', 'can_be_solved', 'contain_selfloop', 'ctc', 'numreads', 'totalbase', 'coverage']
    df_temp_graph_summary = pd.DataFrame(list_identify_results, columns=header)
    df_final_graph_summary = pd.merge(left=df_graph_summary, right=df_temp_graph_summary, left_on='id', right_on='id', how='inner')
    df_final_graph_summary = df_final_graph_summary[selected_header]
    df_final_graph_summary['merge_region'] = df_final_graph_summary['merge_region'].apply(lambda x: format_merge_region(x))

    write_path = "{}/subGraphs.summary.txt".format(outDir)
    df_final_graph_summary.to_csv(write_path, sep='\t', index=None)

    ct = datetime.datetime.now()
    print("[{}] finished preparing data".format(ct), flush=True)

    ## write all parameters for correcting sequences
    params_path = "{}/params.txt".format(tmpDir)
    with open(params_path, 'w') as w_f:
        for tup in list_params:
            write_line = "{}\n".format("\t".join(tup))
            w_f.write(write_line)

    ## run sequence correction
    ct = datetime.datetime.now()
    total_medaka = len(list_params)
    print("[{}] running sequence correction ({} eccDNA candidates, {} processes)".format(ct, total_medaka, threads), flush=True)

    list_medaka_cmds = []
    for tup in list_params:
        cmd = "medaka_consensus -i {} -d {} -o {} -m {} >/dev/null 2>&1".format(tup[0], tup[1], tup[2], tup[3])
        list_medaka_cmds.append(cmd)

    pool = Pool(processes=threads)
    for i, _ in enumerate(pool.imap_unordered(dist_work, list_medaka_cmds), 1):
        if i % 100 == 0 or i == total_medaka:
            ct = datetime.datetime.now()
            print("[{}] sequence correction: {}/{}".format(ct, i, total_medaka), flush=True)
    pool.close()
    pool.join()

    ## copy corrected sequences
    dict_medaka_seq_len = {}
    for tup_consensus in list_chk_consensus:
        consensus_path, dest_path = tup_consensus
        if os.path.isfile(consensus_path):
            shutil.copy(consensus_path, dest_path)

            seq_len = get_fasta_length(dest_path)
            key_seq = dest_path.split('/')[-1].split('_consensus')[0]
            dict_medaka_seq_len[key_seq] = seq_len

    ## run variant calling
    ct = datetime.datetime.now()
    total_variant = len(list_params)
    print("[{}] running variant calling ({} eccDNA candidates, {} processes)".format(ct, total_variant, threads), flush=True)

    list_variant_cmds = []
    for tup in list_params:
        cmd = "medaka variant --quiet {} {} {}".format(tup[1], tup[4], tup[5])
        list_variant_cmds.append(cmd)

    pool = Pool(processes=threads)
    for i, _ in enumerate(pool.imap_unordered(dist_work, list_variant_cmds), 1):
        if i % 100 == 0 or i == total_variant:
            ct = datetime.datetime.now()
            print("[{}] variant calling: {}/{}".format(ct, i, total_variant), flush=True)
    pool.close()
    pool.join()

    ## copy variant results
    for tup_variant in list_chk_variant:
        variant_path, dest_path = tup_variant
        if os.path.isfile(variant_path):
            shutil.copy(variant_path, dest_path)

    ## write identified eccDNA summary
    header = ['id', 'merge_region', 'merge_len', 'num_region', 'ctc', 'numreads', 'totalbase', 'coverage']
    df_identify_summary = pd.DataFrame(list_identify_results, columns=header)
    df_identify_summary['merge_region'] = df_identify_summary['merge_region'].apply(lambda x: format_merge_region(x))
    df_identify_summary['combined_consensus_status'] = df_identify_summary.apply(lambda x: get_consensus_status(x, dict_medaka_seq_len), axis = 1)
    df_identify_summary['consensus_len'] = df_identify_summary['combined_consensus_status'].apply(lambda x: int(x.split('|')[0]))
    df_identify_summary['consensus_status'] = df_identify_summary['combined_consensus_status'].apply(lambda x: x.split('|')[1])
    df_identify_summary.drop(columns=['combined_consensus_status'], inplace=True)

    eccdna_final_path = "{}/eccDNA_final.txt".format(bname)
    df_identify_summary.to_csv(eccdna_final_path, sep='\t', index=None)

    if skipgfa:
        ## create GFA files
        ct = datetime.datetime.now()
        print("[{}] creating GFA files".format(ct), flush=True)

        for row in df_identify_summary.itertuples():
            tup_ = (row.id, row.consensus_len, row.coverage, row.ctc)
            fa_path = "{}/{}/{}_consensus.fa".format(assemGraph, row.id, row.id)
            assemToGFA((fa_path, tup_))

    ct = datetime.datetime.now()
    print("[{}] finished identifying process\n".format(ct), flush=True)
