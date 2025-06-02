import csv
import networkx as nx
import argparse
import matplotlib.pyplot as plt
import pysam
import matplotlib.patches as mpatches
from pybedtools import BedTool

def parse_graphs(file_name, bed_inv):
    graphs = []
    G = nx.DiGraph()
    with open(file_name, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)  # Skip the first header
        for row in reader:
            if row[0] == 'Node_id':  # Start a new graph
                if G.number_of_nodes() > 0:  # If there are nodes in the graph, add it to the list
                    source_node = min(G.nodes, key=lambda x: G.nodes[x]['id'])
                    sink_node = max(G.nodes, key=lambda x: G.nodes[x]['id'])
                    G.nodes[source_node]['source'] = True
                    G.nodes[sink_node]['sink'] = True
                    graphs.append(G)
                G = nx.DiGraph()  # Create a new graph
            else:
                node_id = int(row[0])
                start = int(row[1])
                end = int(row[2])
                coverage = float(row[3])  # Assuming the coverage value is in the third column
                children_ids = [int(child_id) for child_id in row[7].split(',') if child_id]  # split and remove empty strings
                # G.add_node(node_id, id=node_id, start=start, end=end, reads=coverage)  # add the node to the graph with attributes
                G.add_node(node_id, id=node_id, start=start, end=end, reads=coverage, reads_count=0)  # add the node to the graph with attributes
                for child_id in children_ids:  # add edges to the children
                    G.add_edge(node_id, child_id)
        if G.number_of_nodes() > 0:  # Add the last graph if it's not empty
            source_node = min(G.nodes, key=lambda x: G.nodes[x]['id'])
            sink_node = max(G.nodes, key=lambda x: G.nodes[x]['id'])
            G.nodes[source_node]['source'] = True
            G.nodes[sink_node]['sink'] = True
            graphs.append(G)
    return graphs

def parse_graphs2(file_name, bed_inv):
    graphs = []
    G = nx.DiGraph()
    node_id = 0
    source_node = None
    sink_node = None
    with open(file_name, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip headers
        for row in reader:
            node_id = int(row[0])
            start = int(row[1])
            end = int(row[2])
            coverage = float(row[5]) 
            children_ids = [int(child_id) for child_id in row[7].split(',') if child_id]
            length = end - start +1 

            is_inv = False
            if start != 0 and end != 0:  # ignore source and sink nodes for inversion check
                node_bed = BedTool(f"{args.chrom}\t{start}\t{end}", from_string=True)
                intersection = bed_inv.intersect(node_bed)
                if len(intersection) > 0:
                    is_inv = True

            G.add_node(node_id, id=node_id, start=start, end=end, reads=coverage, reads_count=0, is_inv=is_inv, length=length)
            for child_id in children_ids:  
                G.add_edge(node_id, child_id)

        if G.number_of_nodes() > 0: 
            source_node = 0  # based on provided data
            sink_node = node_id
            G.nodes[source_node]['source'] = True
            G.nodes[sink_node]['sink'] = True
            graphs.append(G)
    return graphs

def count_conreads_through_edges(graphs, chrom):
    for G in graphs:
        for (node1_id, node2_id) in G.edges:
            node1_data = G.nodes[node1_id]
            node2_data = G.nodes[node2_id]
            if not node1_data.get('source', False) and not node2_data.get('sink', False):
                print("Processing edge:", node1_data, node2_data)

def parse_juncs(junc_file):
    juncs = []
    with open(junc_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            juncs.append((int(row[0]), int(row[1])))
    return juncs

def add_inv_info(graphs, bed_inv):
    start_pos_to_node = {}
    end_pos_to_node = {}
    # for loop each node in graph
    for G in graphs:
        for node_id in G.nodes:
            node_data = G.nodes[node_id]
            if node_data.get('source', False):
                print("debug source node:", node_data)
                node_data['is_inv'] = False
                continue
            if node_data.get('sink', False):
                node_data['is_inv'] = False
                continue
            # Create a BedTool object for the node
            node_bed = BedTool(f"{args.chrom}\t{node_data['start']}\t{node_data['end']}", from_string=True)
            print("debug node bed:", node_bed, node_id)
            # Check if the node intersects with any region in bed_inv
            intersection = bed_inv.intersect(node_bed)
            if len(intersection) > 0:
                node_data['is_inv'] = True
            else:
                node_data['is_inv'] = False
            start_pos_to_node[node_data['start']] = node_data
            end_pos_to_node[node_data['end']] = node_data
    return graphs, start_pos_to_node, end_pos_to_node

def add_juncs_to_graphs(graphs, juncs, start_pos_to_node, end_pos_to_node):
    # juncs = []
    juncs_dict = {}
    for junc in juncs:
        # JUNC 0-0 + 1-1 - 3
        if junc[0] not in end_pos_to_node or junc[1] not in start_pos_to_node:
            continue
        node1 =  end_pos_to_node[junc[0]]
        node2 = start_pos_to_node[junc[1]]
        seg1= f"{node1['start']}-{node1['end']}"
        seg2= f"{node2['start']}-{node2['end']}"
        sign1 = ""
        sign2 = ""
        sign1 = "-" if node1['is_inv'] else "+"
        sign2 = "-" if node2['is_inv'] else "+"
        junc_str = f"JUNC {seg1} {sign1} {seg2} {sign2}"
        if junc_str not in juncs_dict:
            juncs_dict[junc_str] = 1
        else:
            juncs_dict[junc_str] += 1
    return juncs_dict

def write_graphs(graphs, juncs, file_name):
    of = open(file_name, 'w')
    for G in graphs:
        for node_id in G.nodes:
            node_data = G.nodes[node_id]
            if node_data.get('source', False) or node_data.get('sink', False):
                continue
            # print(node_data)
            node_str = f"SEG {node_data['start']}-{node_data['end']} {node_data['coverage']} 1"
            of.write(node_str + '\n')
        for u, v in G.edges:
            node_u = G.nodes[u]
            node_v = G.nodes[v]
            if node_u.get('source', False) or node_v.get('sink', False):
                continue
            node_u_str = f"{node_u['start']}-{node_u['end']}"
            node_v_str = f"{node_v['start']}-{node_v['end']}"
            sign_u = "-" if node_u['is_inv'] else "+"
            sign_v = "-" if node_v['is_inv'] else "+"
            edge_str = f"JUNC {node_u_str} {sign_u} {node_v_str} {sign_v} {G[u][v]['support']}"
            of.write(edge_str + '\n')



        # for junc in juncs:
        #     of.write(f"{junc} {juncs[junc]}\n")
    of.close()    

def visualize_graph(G):
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G)  # positions for all nodes

    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=700)

    # edges
    nx.draw_networkx_edges(G, pos, width=6)

    # labels
    nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

    plt.axis('off')
    plt.show()


def update_coverage_and_edge_support(graphs, bam_file_name, start_pos_to_node, end_pos_to_node):
    bam_file = pysam.AlignmentFile(bam_file_name, 'rb')

    # Initialize coverage and edge support dictionaries
    coverage = {node_id: [0]*G.nodes[node_id]['length'] for G in graphs for node_id in G.nodes}
    edge_support = {(u, v): 0 for G in graphs for u, v in G.edges}
    all_junc_cnt = 0
    used_junc_cnt = 0
    # Iterate over each read in the BAM file
    for read in bam_file.fetch():
        blocks = read.get_blocks()
        # print(blocks)
        new_blocks = []
        for block in blocks:
            new_blocks.append((block[0]+1, block[1]))

        # Update node coverage
        for block_start, block_end in new_blocks:
            for G in graphs:
                for node_id, node_data in G.nodes(data=True):
                    node_start = node_data['start']
                    node_end = node_data['end']
                    
                    if node_start >= block_start and node_end >= block_end:
                        overlap_start = 0
                        overlap_end = block_end - node_start + 1
                    elif node_start <= block_start and node_end >= block_end:
                        overlap_start = block_start - node_start
                        overlap_end = block_end - node_start + 1
                    elif node_start <= block_start and node_end <= block_end:
                        overlap_start = block_start - node_start
                        overlap_end = node_end - node_start +1
                    elif node_start >= block_start and node_end <= block_end:
                        
                        overlap_start = 0
                        overlap_end = node_end - node_start+1
                        # if node_id == 2:
                        #     print("debug node_start node_end:", node_start, node_end, node_id, block_start, block_end, overlap_start, overlap_end)
                    else:
                        continue
                    for i in range(overlap_start, overlap_end):
                        coverage[node_id][i] += 1

                    # if node_start <= block_start and node_end >= block_end:
                    #     for i in range(block_start - node_start, block_end - node_start):
                    #         coverage[node_id][i] += 1

        # Update edge support for splice junctions
        for i in range(len(new_blocks) - 1):
            start_block_end = new_blocks[i][1]
            end_block_start = new_blocks[i + 1][0]

            all_junc_cnt += 1
            # print("debug start_block_end end_block_start:", start_block_end, end_block_start, new_blocks)
            if start_block_end in end_pos_to_node and end_block_start in start_pos_to_node:
                used_junc_cnt += 1
                start_node = end_pos_to_node[start_block_end]
                end_node = start_pos_to_node[end_block_start]
                edge_support[(start_node['id'], end_node['id'])] += 1
            else:
                print("not used junc:", new_blocks[i], new_blocks[i+1])
        for block in new_blocks:
            for u, v in graphs[0].edges:
                if (block[0] <= graphs[0].nodes[u]['end'] and block[1] >= graphs[0].nodes[u]['end']) and \
                    (block[0] <= graphs[0].nodes[v]['start'] and block[1] >= graphs[0].nodes[v]['start']):
                    edge_support[(u, v)] += 1
                    #
        # update edge support for cut site junctions
    bam_file.close()
    print("debug coverage:", coverage)
    print("debug edge_support:", edge_support)
    # Update graphs with coverage and edge support
    for G in graphs:
        for node_id in G.nodes:
            # Calculate average coverage for each node
            G.nodes[node_id]['coverage'] = sum(coverage[node_id]) / len(coverage[node_id])
        for u, v in G.edges:
            G[u][v]['support'] = edge_support[(u, v)]
    print("debug all_junc_cnt used_junc_cnt:", all_junc_cnt, used_junc_cnt)

    return graphs

def main():
    
    bed_inv = BedTool(args.inv_bed)
    print(bed_inv[0].start, type(bed_inv))
    # Use the function
    graphs = parse_graphs2(args.file_name, bed_inv)
    graphs, start_pos_to_node, end_pos_to_node = add_inv_info(graphs, bed_inv)
    juncs = parse_juncs(args.junc_file)
    # juncs = add_juncs_to_graphs(graphs, juncs, start_pos_to_node, end_pos_to_node)
    

    update_coverage_and_edge_support(graphs, args.bam_file, start_pos_to_node, end_pos_to_node)
    write_graphs(graphs, juncs, args.output)
    # for g in graphs:
    #     for edge in g.edges:
    #         print(g.nodes[edge[0]], g.nodes[edge[1]])
    #     print(g[4][5])
    # for i, G in enumerate(graphs):
    #     visualize_graph(G)
    # print(f'Parsed {len(graphs)} graphs from {args.file_name}.')
    # count_conreads_through_edges(graphs, args.chrom)

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse a graph file.')
    parser.add_argument('file_name', type=str, help='The name of the file to parse.')
    parser.add_argument('inv_bed', type=str, help='inv bed.')
    parser.add_argument('junc_file', type=str, help='junc file.')
    parser.add_argument('chrom', type=str, help='Chromosome name.')
    parser.add_argument('bam_file', type=str, help='BAM file.')
    parser.add_argument('-o', '--output', help='Output file', required=True)
    args = parser.parse_args()
    main()