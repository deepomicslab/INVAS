import csv
import networkx as nx
import argparse
import matplotlib.pyplot as plt
import pysam
import matplotlib.patches as mpatches


def parse_graphs(file_name):
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
                coverage = float(row[5])  # Assuming the coverage value is in the third column
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


def count_and_print_reads_in_interval(samfile, chrom, pos1, pos2):
    # samfile = pysam.AlignmentFile(bam_file, "rb")
    count = 0
    total_reads = 0
    for read in samfile.fetch(chrom, pos1-1, pos2):
        if read.is_unmapped:
            continue
        # Ensure the read covers the whole interval
        # print(read.query_name, read.reference_start, read.reference_end, read.cigarstring, pos1, pos2)
        if read.reference_start <= pos1 and read.reference_end >= pos2:
            total_reads += 1
            # Check the positions are not clipped
            # print(read.query_name, read.reference_start, read.preference_end, read.cigarstring)

            cigar_tuples = read.cigartuples
            read_start = read.reference_start+1
            pos1_matched = False
            pos2_matched = False
            for (op, length) in cigar_tuples:
                # print(op, length)
                if op not in [0, 7]:  # 0 and 7 correspond to match operations
                    read_start += length
                    continue
                # print("parse cigar:",read_start, pos1, read_start + length)
                if read_start <= pos1 < read_start + length:
                    pos1_matched = True
                if read_start <= pos2 < read_start + length:
                    pos2_matched = True
                read_start += length
            if pos1_matched and pos2_matched:
                # print(read.query_name, read.cigarstring)
                count += 1
            # else:
                # print("Read {} does not cover the interval [{}, {}]".format(read.query_name, pos1, pos2))
    print(pos1, pos2, count, total_reads)
    return count, total_reads

def mapped_node(read, node):
    if (read.reference_start <= node['end'] and read.reference_start >= node['start']) or \
       (read.reference_end <= node['end'] and read.reference_end >= node['start']) or \
        (read.reference_start <= node['start'] and read.reference_end >= node['end']):
        return True
    return False

def get_supplementary_alignments(read):
    try:
        sa_tag = read.get_tag('SA')
    except KeyError:
        return None

    sa_fields = sa_tag.split(';')
    supplementary_alignments = []

    for field in sa_fields:
        if field:
            chrom, pos, strand, cigar, mapq, nm = field.split(',')
            supplementary_alignments.append({
                'chrom': chrom,
                'pos': int(pos),
                'strand': strand,
                'cigar': cigar,
                'mapq': int(mapq),
                'nm': int(nm)
            })

    return supplementary_alignments
def check_connection(samfile, chrom, node1, node2):
    # check directed conected edges and inverted conected edges
    # p_p_count : 5+'-- 3+' |----->| 5+'-- 3+' or  3-'---5-' |<-----| 3-'---5-'
    # p_n_count : 5+'-- 3+' |----->| 5-'---3-' or  3-'---5-' |<-----| 3+'-- 5+'
    # p_p_count : one read map two edges
    # p_n_count : one read map two edges with different direction
    # Initialize the counts
    p_p_count = 0
    p_n_count = 0
    read_orientation = "+"
    for read in samfile.fetch(chrom, node1['start'], node2['end']):
        # if read.is_supplementary:
        #     print(read, 'supp')
        if read.is_unmapped:
            continue
        

    return p_p_count, p_n_count


def count_conreads_through_edges(graphs, bam_file, chrom):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for G in graphs:
        for (node1_id, node2_id) in G.edges:
            node1_data = G.nodes[node1_id]
            node2_data = G.nodes[node2_id]
            if not node1_data.get('source', False) and not node2_data.get('sink', False):
                print("Processing edge:", node1_data, node2_data)
                p_p_count, p_n_count = check_connection(samfile, chrom, node1_data, node2_data)
                print(p_p_count, p_n_count)


def main():
    parser = argparse.ArgumentParser(description='Parse a graph file.')
    parser.add_argument('file_name', type=str, help='The name of the file to parse.')
    parser.add_argument('bam_file', type=str, help='Path to BAM file.')
    parser.add_argument('chrom', type=str, help='Chromosome name.')
    args = parser.parse_args()

    # Use the function
    graphs = parse_graphs(args.file_name)
    print(f'Parsed {len(graphs)} graphs from {args.file_name}.')
    # count_reads_through_edges(graphs, args.bam_file, args.chrom)
    count_conreads_through_edges(graphs, args.bam_file, args.chrom)
    # draw_graphs2(graphs)
    # print(graphs[0].nodes, len(graphs[0].nodes))
    # print(graphs[0].edges, len(graphs[0].edges))
    # for node in graphs[0].nodes:
    #     print(graphs[0].nodes[node])
    # for edge in graphs[0].edges:
    #     print(f"Edge from {edge[0]} to {edge[1]} with attributes {graphs[0].edges[edge]}")

if __name__ == '__main__':
    main()