import pysam
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import csv
import sys
from collections import defaultdict

def parse_gtf(gtf_file):
    exons = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'exon':
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                exon = (start, end)
                exons.setdefault(chrom, []).append(exon)
    for chrom in exons:
        exons[chrom].sort(key=lambda x: x[0])
    return exons

def find_nearest_exon(exons, chrom, position):
    nearest_exon = None
    min_dist = float('inf')
    for exon in exons.get(chrom, []):
        exon_start, exon_end = exon
        if exon_end < position:
            dist = position - exon_end
        elif exon_start > position:
            dist = exon_start - position
        else:
            return exon
        if dist < min_dist:
            min_dist = dist
            nearest_exon = exon
    return nearest_exon

def create_exon_chain(exons, chrom, blocks):
    chain = []
    for block_start, block_end in blocks:
        start_exon = find_nearest_exon(exons, chrom, block_start)
        end_exon = find_nearest_exon(exons, chrom, block_end)
        # Choose the closer exon for the start of the block
        if start_exon and end_exon:
            if abs(start_exon[0] - block_start) < abs(end_exon[1] - block_end):
                chosen_exon = start_exon
            else:
                chosen_exon = end_exon
        elif start_exon:
            chosen_exon = start_exon
        elif end_exon:
            chosen_exon = end_exon
        else:
            continue
        # Append if not already in the chain
        if not chain or chain[-1] != chosen_exon:
            chain.append(chosen_exon)
    return tuple(chain)

def is_alternative_splicing(blocks):
    # Check if there are multiple blocks separated by at least the minimum intron size
    # and less than the maximum intron size.
    for i in range(len(blocks) - 1):
        intron_size = blocks[i+1][0] - blocks[i][1]
        if intron_size > min_intron_size and intron_size < max_intron_size:
            return True
    return False

def deduplicate_intervals(intervals_dict, error_margin=1):
    deduplicated_dict = {}
    for chromosome, intervals in intervals_dict.items():
        deduplicated_intervals = []
        for interval in sorted(intervals):  # Sort intervals to compare them
            if not deduplicated_intervals:
                # If list is empty, add the first interval
                deduplicated_intervals.append(interval)
            else:
                # Get the last interval in the deduplicated list
                last_interval = deduplicated_intervals[-1]
                # Check if the current interval is within the error_margin
                if (abs(interval[0] - last_interval[0]) <= error_margin and
                        abs(interval[1] - last_interval[1]) <= error_margin):
                    # If it is within the error_margin, skip it (it's a duplicate)
                    continue
                else:
                    # Otherwise, add the current interval to the list
                    deduplicated_intervals.append(interval)
        # Add the deduplicated list of intervals to the dictionary
        deduplicated_dict[chromosome] = deduplicated_intervals
    return deduplicated_dict

def parse_chain_from_bam(bam, exons):
    exon_chains = {}
    all_rec_cnt=0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        chrom = bam.get_reference_name(read.reference_id)
        blocks = read.get_blocks()
        if not is_alternative_splicing(blocks):
            continue
        all_rec_cnt+=1
        # if len(blocks) >= 2:
        #     print(blocks, read.query_name)
        exon_chain = create_exon_chain(exons, chrom, blocks)
        
        if exon_chain:
            exon_chains[exon_chain] = exon_chains.get(exon_chain, 0) + 1
    print(f"Total records: {all_rec_cnt}")
    
    return exon_chains

def parse_inv(inv_bed):
    # file format :
    # chr22:16390136-16390220	test:16387694-16390887,test	chr22:16390136-16390220-15.33
    # parse the first col, get the inv region, and store in a dict, key is chrom, value is a list of (start, end)
    # some duplicated inv regions, only keep one
    inv_rgs = {}
    with open(inv_bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            inv_rg_str = parts[0]
            chrom = inv_rg_str.split(':')[0]
            start = int(inv_rg_str.split(':')[1].split('-')[0])
            end = int(inv_rg_str.split(':')[1].split('-')[1])
            # use set to store the inv regions
            inv_rgs.setdefault(chrom, set()).add((start, end))

    return inv_rgs

# def calculate_region_depth(bam, chrom, start, end):
#     # 打开BAM文件
#     # 初始化总深度为0
#     total_depth = 0
#     # 初始化覆盖的碱基数量为0
#     covered_bases = 0
    
#     # 遍历指定区域的pileup列
#     for pileupcolumn in bam.pileup(chrom, start, end, nofilter=True):
#         # 确保pileup列在指定的区域内
#         if pileupcolumn.reference_pos >= start and pileupcolumn.reference_pos < end:
#             # 累加该位置的读段数量到总深度
#             total_depth += pileupcolumn.nsegments
#             # 增加覆盖的碱基数量
#             covered_bases += 1
    
#     # 如果有覆盖的碱基，则计算平均深度
#     if covered_bases > 0:
#         average_depth = total_depth / covered_bases
#     else:
#         average_depth = 0

#     return average_depth
def calculate_region_depth(bam, chrom, start, end):
    # 使用count_coverage方法获取指定区域的覆盖度信息
    coverage = bam.count_coverage(chrom, start, end, quality_threshold=0)
    total_depth = 0
    covered_bases = 0

    # coverage包含四个元素的元组，每个元素对应一个碱基（A, C, G, T）的覆盖度列表
    for base_idx in range(4):  # 遍历每个碱基
        for pos_depth in coverage[base_idx]:  # 遍历该碱基对应的覆盖度列表
            total_depth += pos_depth  # 累加到总深度
            if pos_depth > 0:
                covered_bases += 1  # 如果该位置有覆盖，则增加覆盖的碱基数量

    # 如果有覆盖的碱基，则计算平均深度
    if covered_bases > 0:
        average_depth = total_depth / covered_bases
    else:
        average_depth = 0

    return average_depth

def gen_seg(exons):
    # seg format: SEG 16387697-16387840 34.24305555555556 1
    # finaldepth equals recover depth + remap depth
    segs_strs= []
    for chrom in exons:
        for exon in exons[chrom]:
            exon_start, exon_end = exon
            recover_depth = calculate_region_depth(recover_bam, chrom, exon_start, exon_end)
            remap_depth = calculate_region_depth(remap_bam, chrom, exon_start, exon_end)
            final_depth = recover_depth*4 + remap_depth
            print(exon, remap_depth, recover_depth, final_depth)
            segs_strs.append(f"SEG {exon_start}-{exon_end} {final_depth} 1\n")
    return segs_strs

def gen_junc(recover_exon_chains, remap_exon_chains, inv_exons):
    # JUNC 16387697-16387840 + 16390137-16390137 - 24
    # process recover exon chains
    juncs_strs = []
    inv_exons_lst = inv_exons["chr"+str(chrom)]
    juncs_dict = {}
    # for recovered bam
    for chain, count in recover_exon_chains.items():
        if len(chain) == 1:
            continue
        chain_len = len(chain)
        for i in range(len(chain) - 1):
            exon1_start, exon1_end = chain[i]
            exon2_start, exon2_end = chain[i+1]
            if chain[i] in inv_exons_lst and chain[i+1] not in inv_exons_lst:
                junc_str = f"JUNC {exon1_start}-{exon1_end} - {exon2_start}-{exon2_end} +"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
            elif chain[i] not in inv_exons_lst and chain[i+1] in inv_exons_lst:
                junc_str = f"JUNC {exon1_start}-{exon1_end} + {exon2_start}-{exon2_end} -"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
            elif chain[i] in inv_exons_lst and chain[i+1] in inv_exons_lst:
                junc_str = f"JUNC {exon1_start}-{exon1_end} - {exon2_start}-{exon2_end} -"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
            else:
                junc_str = f"JUNC {exon1_start}-{exon1_end} + {exon2_start}-{exon2_end} +"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
    # for remapped bam
    for chain, count in remap_exon_chains.items():
        if len(chain) == 1:
            continue
        for i in range(len(chain) - 1):
            exon1_start, exon1_end = chain[i]
            exon2_start, exon2_end = chain[i+1]
            junc_str = f"JUNC {exon1_start}-{exon1_end} + {exon2_start}-{exon2_end} +"
            if junc_str in juncs_dict:
                juncs_dict[junc_str] += count
            else:
                juncs_dict[junc_str] = count
    for junc_str, count in juncs_dict.items():
        juncs_strs.append(f"{junc_str} {count}\n")

    return juncs_strs

def wite_graph(segs, juncs, output_file):
    with open(output_file, 'w') as f:
        for seg in segs:
            f.write(seg)
        for junc in juncs:
            f.write(junc)           

def parse_graph_file(filepath):
    with open(filepath, 'r') as file:
        graphs = []
        current_graph = None
        
        for line in file:
            # Remove any leading/trailing whitespace
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Check if the line is the start of a new graph
            if line.startswith('Graph:'):
                # If we are already building a graph, add it to the graphs list
                if current_graph is not None:
                    graphs.append(current_graph)
                # Start a new graph
                current_graph = []
            else:
                start, end = line.split('\t')[0].split('-')
                cov = float(line.split('\t')[1])
                current_graph.append((int(start), int(end), cov))
        
        # Don't forget to add the last graph to the list if it exists
        if current_graph is not None:
            graphs.append(current_graph)
    
    return graphs

def parse_transfrg_file(filepath):
    with open(filepath, 'r') as file:
        graphs = []
        current_graph = None
        
        for line in file:
            # Remove any leading/trailing whitespace
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Check if the line is the start of a new graph
            if line.startswith('Graph:'):
                # If we are already building a graph, add it to the graphs list
                if current_graph is not None:
                    graphs.append(current_graph)
                # Start a new graph
                current_graph = []
            else:
                strand, path, aboundance = line.split('\t')
                strand = int(strand)
                aboundance = float(aboundance)
                path = [int(x) for x in path.split(',') if x!='']
                
                current_graph.append((strand, path, aboundance))
        
        # Don't forget to add the last graph to the list if it exists
        if current_graph is not None:
            graphs.append(current_graph)
    
    return graphs

def assign_in_out_through_nodes(gnodes, gtransfrag):
    # gnodes: [(start, end, cov), ...]
    # gtransfrag: [(strand, path, aboundance), ...]
    # nodes_ab_status: {node: (in_aboundance, out_aboundance)}
    node_status= {}
    for gi, graph in enumerate(gnodes):
        node_status[gi] = {}
        for ni, node in enumerate(graph):
            if ni == 0 or ni == len(graph)-1:
                continue
            node_status[gi][ni]= {
                "in_aboundance" : 0.0,
                "out_aboundance" : 0.0,
                "through_aboundance" : 0.0
                }
    for gi, graph in enumerate(gtransfrag):
        for strand, path, aboundance in graph:
            for idx, node in enumerate(path):
                if node in [0, len(gnodes[gi])-1]:
                    continue
                if node == path[0]:
                    node_status[gi][node]["out_aboundance"] += aboundance
                elif node==path[-1]:
                    node_status[gi][node]["in_aboundance"] += aboundance
                else:
                    node_status[gi][node]["through_aboundance"] += aboundance

    return node_status

def assign_in_out_through_nodes_inv(gnodes, inv_transfrg, nodes_ab_status):
    gi=0
    for strand, path, aboundance in inv_transfrg:
        for idx, node in enumerate(path):
            if node in [0, len(gnodes[gi])-1]:
                continue
            if node == path[0]:
                nodes_ab_status[gi][node]["out_aboundance"] += aboundance
            elif node==path[-1]:
                nodes_ab_status[gi][node]["in_aboundance"] += aboundance
            else:
                nodes_ab_status[gi][node]["through_aboundance"] += aboundance
    return nodes_ab_status

def generate_normal_edges(gnodes, gtransfrag):
    g_edges=[]
    for gi, graph in enumerate(gtransfrag):
        sub_g_edges = {}
        for strand, paths, aboundance in graph:
            # genrate edges for neibour nodes
            for pi, p in enumerate(paths):
                if pi == len(paths) - 1:
                    break
                if p in [0, len(gnodes[gi])-1 ] or paths[pi+1] in [0, len(gnodes[gi])-1 ]:
                    continue
                exon1_start, exon1_end, cov1 = gnodes[gi][p]
                exon2_start, exon2_end, cov2 = gnodes[gi][paths[pi+1]]
                exon1_str = f"{exon1_start}-{exon1_end}"
                exon2_str = f"{exon2_start}-{exon2_end}"
                edge_str = f"JUNC {exon1_str} + {exon2_str} +"

                if edge_str in sub_g_edges:
                    sub_g_edges[edge_str] += aboundance
                else:
                    sub_g_edges[edge_str] = aboundance
        g_edges.append(sub_g_edges)
    return g_edges

def generate_normal_graph(gnodes, gtransfrag):
    # generate network flow graph fron gtransfrag
    for gi, graph in enumerate(gtransfrag):
        sub_g_edges = {}
        for strand, paths, aboundance in graph:
            # genrate edges for neibour nodes
            for pi, p in enumerate(paths):
                if pi == len(paths) - 1:
                    break
                
def generate_normal_graph(gnodes, gtransfrag):
    # Create a new directed graph
    G = nx.DiGraph()
    
    # Add nodes to the graph
    G.add_nodes_from(range(len(gnodes[0])))
    
    # Generate network flow graph from gtransfrag
    for gi, graph in enumerate(gtransfrag):
        for strand, paths, abundance in graph:
            # Generate edges for neighboring nodes
            for pi in range(len(paths) - 1):
                # Add an edge from the current path node to the next path node
                # with 'abundance' as the capacity of the edge
                G.add_edge(paths[pi], paths[pi + 1], capacity=abundance)
    # Draw the graph
    pos = nx.spring_layout(G)  # positions for all nodes
    nx.draw(G, pos, with_labels=True)
    edge_labels = nx.get_edge_attributes(G, 'capacity')
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    # Show the plot
    # plt.show()

    sink_node = max(G.nodes)
    # 获取直接连接到sink的所有节点（前驱节点）
    predecessors = list(G.predecessors(sink_node))

    print(f"The sink node is: {sink_node}")
    print("Nodes connected to the sink node:", predecessors)
    
    return G, predecessors


def modify_last_exons(nodes_ab_status, last_exon_ids):
    for gi, graph in nodes_ab_status.items():
        for ni, node in graph.items():
            if ni in last_exon_ids:
                nodes_ab_status[gi][ni]['cap'] = node['in_aboundance']+node['through_aboundance']
            else:
                nodes_ab_status[gi][ni]['cap'] = node['out_aboundance']+node['through_aboundance']

    return nodes_ab_status


            
def write_graph_file(nodes, normal_edges, output_file):
    with open(output_file, 'w') as f:
        for node in nodes:
            f.write(f"{node}\n")
        for edge, val in normal_edges.items():
            f.write(f"{edge} {val}\n")

def main():

    # process for each candidate inv region
    # so the chromosome is fixed
    exons = parse_gtf(args.gtf_file)
    exons = deduplicate_intervals(exons, error_margin=5)
    print(exons)
    recover_exon_chains = parse_chain_from_bam(recover_bam, exons)
    inv_rg = parse_inv(args.inv_bed)
    # get exons in inv regions
    inv_exons = {}
    for chrom in inv_rg:
        chrom_nochr = chrom.replace('chr', '')
        inv_exons[chrom] = []
        for inv_start, inv_end in inv_rg[chrom]:
            for exon in exons[chrom_nochr]:
                exon_start, exon_end = exon
                if exon_start >= inv_start and exon_end <= inv_end:
                    inv_exons[chrom].append(exon)
    # print(inv_exons)
    # Print the exon chains and their counts
    # for chain, count in recover_exon_chains.items():
    #     print(f"Exon chain {chain} occurs {count} times")
    # Print the exon chains and their counts

    segs = gen_seg(exons)

    # link_to_graph(recover_bam, remap_bam, recover_exon_chains, remap_exon_chains, inv_exons)

def gen_inv_edges(recover_exon_chains):
    juncs_strs = []
    inv_exons_lst = inv_exons[str(chrom)]
    juncs_dict = {}
    # for recovered bam
    for chain, count in recover_exon_chains.items():
        if len(chain) == 1:
            continue
        chain_len = len(chain)
        for i in range(len(chain) - 1):
            exon1_start, exon1_end = chain[i]
            exon2_start, exon2_end = chain[i+1]
            if chain[i] in inv_exons_lst and chain[i+1] not in inv_exons_lst:
                junc_str = f"JUNC {exon1_start}-{exon1_end} - {exon2_start}-{exon2_end} +"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
            elif chain[i] not in inv_exons_lst and chain[i+1] in inv_exons_lst:
                junc_str = f"JUNC {exon1_start}-{exon1_end} + {exon2_start}-{exon2_end} -"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
            elif chain[i] in inv_exons_lst and chain[i+1] in inv_exons_lst:
                junc_str = f"JUNC {exon1_start}-{exon1_end} - {exon2_start}-{exon2_end} -"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
            else:
                junc_str = f"JUNC {exon1_start}-{exon1_end} + {exon2_start}-{exon2_end} +"
                if junc_str in juncs_dict:
                    juncs_dict[junc_str] += count/theta[chain_len]
                else:
                    juncs_dict[junc_str] = count/theta[chain_len]
    return juncs_dict

def merge_edges(normal_edges, inversion_edges):
    merged_dict = normal_edges.copy()  # Start with the first dictionary

    for key, value in inversion_edges.items():
        if key in merged_dict:
            merged_dict[key] += value  # Sum values if key is found in both dictionaries
        else:
            merged_dict[key] = value  # Add the new key-value pair from the second dictionary
    return merged_dict

def gen_segs(gnodes, recover_bam):
    for gi, graph in enumerate(gnodes):
        for ni, node in enumerate(graph):
            seg_str = ""

def chain_to_frg(gnodes, recover_exon_chains, seg_to_idx):
    res = []
    for chain, count in recover_exon_chains.items():
        path = []
        for seg in chain:
            path.append(seg_to_idx[seg])
        res.append((1, path, count/theta[len(chain)]))

    
    return res

def gen_final_node_str(gnodes, nodes_ab_status, recover_bam):
    seg_strs = []
    for gi, graph in nodes_ab_status.items():
        for ni, node in graph.items():
            exon_start, exon_end, exon_cov = gnodes[gi][ni]
            recover_depth = calculate_region_depth(recover_bam, chrom_nochr, exon_start, exon_end)
            node_str = f"SEG {exon_start}-{exon_end} {node['cap']} 1 {exon_cov+recover_depth*2}"
            seg_strs.append(node_str)
    return seg_strs

def modify_last_exons(nodes_ab_status, last_exon_ids):
    for gi, graph in nodes_ab_status.items():
        for ni, node in graph.items():
            if ni in last_exon_ids:
                nodes_ab_status[gi][ni]['cap'] = node['in_aboundance']+node['through_aboundance']
            else:
                nodes_ab_status[gi][ni]['cap'] = node['out_aboundance']+node['through_aboundance']

    return nodes_ab_status
      
def merge_data(input_file_path):
    # 使用defaultdict来存储合并的数据
    all_nodes= []
    transfrgs = defaultdict(float)
    strand = 0
    # 读取输入文件
    with open(input_file_path, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue  # 跳过格式不正确的行
            paths_str = row[1].strip(',')
            paths = paths_str.split(',')
            strand = row[0].strip()
            value = float(row[2].strip())  # 第三列转换为浮点数
            transfrgs[paths_str] += value
            for p in paths:
                if int(p) not in all_nodes:
                    all_nodes.append(int(p))

    # sort all_nodes
    all_nodes.sort()

    # trans segid to matrix id
    new_paths = {}
    for k, v in transfrgs.items():
        paths = k.split(',')
        new_path_str = ""
        for pi, p in enumerate(paths):
            if pi == len(paths)-1:
                break
            n1 = int(p)
            n2 = int(paths[pi+1])
            if n1 == 0:
                n1_me = (len(all_nodes)-2)*4+1
                n2_ms = n2*4-4
                n2_me = n2*4-3
                new_path_str += f"{n1_me},{n2_ms},{n2_me},"


            elif n2 == all_nodes[-1]:
                n2_ms = (len(all_nodes)-2)*4+2
                n1_ms = n1*4-4
                n1_me = n1*4-3
                if pi == 1:
                    new_path_str += f"{n2_ms},"
                else:
                    new_path_str += f"{n1_ms},{n1_me},{n2_ms},"
            else:
                n1_ms = n1*4-4
                n1_me = n1*4-3
                n2_ms = n2*4-4
                n2_me = n2*4-3
                if pi ==0:
                    new_path_str += f"{n1_ms},{n1_me},{n2_ms},{n2_me},"
                else:
                    new_path_str += f"{n2_ms},{n2_me},"
        new_paths[new_path_str]=v
    return new_paths



def transf_to_file(inv_transfrg, output_trans):
    normal_paths = merge_data(args.transfrg_file)
    inv_paths={}
    for invf in inv_transfrg:
        new_path_str = ""
        if len(inv_exons_ids[chrom]) ==1:
            for pi, p in enumerate(invf[1]):
                if pi == len(invf[1])-1:
                    break
                n1 = int(p)
                n2 = int(invf[1][pi+1])
                if n1 in inv_exons_ids[chrom] and n2 not in inv_exons_ids[chrom]:
                    n1_ms = n1*4-2
                    n1_me = n1*4-1
                    n2_ms = n2*4-4
                    n2_me = n2*4-3
                    if pi==0:
                        new_path_str += f"{n1_ms},{n1_me},{n2_ms},{n2_me},"
                    else:
                        new_path_str += f"{n2_ms},{n2_me},"
                elif n1 not in inv_exons_ids[chrom] and n2 in inv_exons_ids[chrom]:
                    n1_ms = n1*4-4
                    n1_me = n1*4-3
                    n2_ms = n2*4-2
                    n2_me = n2*4-1
                    if pi==0:
                        new_path_str += f"{n1_ms},{n1_me},{n2_ms},{n2_me},"
                    else:
                        new_path_str += f"{n2_ms},{n2_me},"
                elif n1 in inv_exons_ids[chrom] and n2 in inv_exons_ids[chrom]:
                    n1_ms = n1*4-2
                    n1_me = n1*4-1
                    n2_ms = n2*4-2
                    n2_me = n2*4-1
                    if pi==0:
                        new_path_str += f"{n2_ms},{n2_me},{n1_ms},{n1_me},"
                    else:
                        new_path_str += f"{n1_ms},{n1_me},"
                
                
            inv_paths[new_path_str]=float(invf[2])
    # print(inv_transfrg)
    print(inv_paths)
    with open(output_trans, 'w') as f:
        for k, v in normal_paths.items():
            f.write(f"{1} {k} {v}\n")
        for k, v in inv_paths.items():
            f.write(f"{1} {k} {v}\n")

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Track exon chains in sequencing reads')
    parser.add_argument('recover_bam_file', help='Path to the BAM file')
    parser.add_argument('gtf_file', help='Path to the output file')
    parser.add_argument('seg_file', help='Path to the GTF file')
    parser.add_argument('inv_bed', help='Path to the BAM file')
    parser.add_argument('transfrg_file', help='Path to the output file')
    parser.add_argument('output_file', help='Path to the output file')
    parser.add_argument('output_trans', help='Path to the output file')


    min_intron_size=50
    max_intron_size=100000
    theta={2:0.3,3:0.3}
    
    # Parse the arguments
    args = parser.parse_args()
    gnodes = parse_graph_file(args.seg_file)
    # for normal part
    gtransfrag = parse_transfrg_file(args.transfrg_file)
    nodes_ab_status = assign_in_out_through_nodes(gnodes, gtransfrag)
    seg_to_idx = {}
    for gi, graph in enumerate(gnodes):
        for ni, node in enumerate(graph):
            seg_to_idx[(node[0], node[1])]=ni

            # print(ni, node)
    # for gi, graph in enumerate(gtransfrag):
    #     for strand, paths, ab in graph:
    #         print(paths, ab)
    normal_edges = generate_normal_edges(gnodes, gtransfrag)
    normal_g, last_exon_ids = generate_normal_graph(gnodes, gtransfrag)

    # for inv part 
    recover_bam=pysam.AlignmentFile(args.recover_bam_file, "rb")
    chrom = recover_bam.references[0]
    exons = parse_gtf(args.gtf_file)
    exons = deduplicate_intervals(exons, error_margin=5)
    print("exons:", exons)
    recover_exon_chains = parse_chain_from_bam(recover_bam, exons)
    inv_rg = parse_inv(args.inv_bed)
    inv_exons = {}
    inv_exons_ids = {}
    chrom_nochr = ""
    for chrom in inv_rg:
        chrom_nochr = chrom
        inv_exons[chrom] = []
        inv_exons_ids[chrom] = []
        for inv_start, inv_end in inv_rg[chrom]:
            for exon_id, exon in enumerate(exons[chrom_nochr]):
                exon_start, exon_end = exon
                if exon_start >= inv_start and exon_end <= inv_end:
                    inv_exons[chrom].append(exon)
                    inv_exons_ids[chrom].append(exon_id+1)
    # print(inv_exons)
    # Print the exon chains and their counts
    print("recover_exon_chains", recover_exon_chains)
    print(seg_to_idx)
    inv_transfrg = chain_to_frg(gnodes, recover_exon_chains, seg_to_idx)
    # print(nodes_ab_status)
    nodes_ab_status_full = assign_in_out_through_nodes_inv(gnodes, inv_transfrg, nodes_ab_status)
    # print(nodes_ab_status_full)


    # for chain, count in recover_exon_chains.items():
    #     print(f"Exon chain {chain} occurs {count} times")
    # # generate edges for inversion
    inversion_edges = gen_inv_edges(recover_exon_chains)
    full_edges = merge_edges(normal_edges[0], inversion_edges)
    # print(inversion_edges)
    # print(normal_edges)
    # print(full_edges)
    # print(gtransfrag)
    # print(inv_transfrg)
    # full_segs=gen_segs(gnodes, recover_bam)
    nodes_ab_status_full = modify_last_exons(nodes_ab_status_full, last_exon_ids)
    nodes = gen_final_node_str(gnodes, nodes_ab_status_full, recover_bam)
    # print(nodes)
    write_graph_file(nodes, full_edges, args.output_file)
    transf_to_file(inv_transfrg, args.output_trans)
    




    # for gi, graph in enumerate(normal_edges):
    #     for edge, ab in graph.items():
    #         print(edge, ab)


    # for gi, graph in nodes_ab_status.items():
    #     for ni, node in graph.items():
    #         print(f"node {ni}: {node}, rate: {node['out_aboundance']/node['in_aboundance']}, cap: {node['out_aboundance']+node['through_aboundance']}, cap2: {node['in_aboundance']+node['through_aboundance']}")
    


    # recover_bam=pysam.AlignmentFile(args.recover_bam_file, "rb")
    # chrom = recover_bam.references[0]
    
    # Run the main function with the provided arguments
    # main()