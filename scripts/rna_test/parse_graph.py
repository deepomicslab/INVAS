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
                G.add_node(node_id, id=node_id, start=start, end=end, reads=coverage)  # add the node to the graph with attributes
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

# def check_conection(samfile, chrom, node1, node2):
#     # check directed conected edges and inverted conected edges
#     # p_p_count : 5+'-- 3+' |----->| 5+'-- 3+' or  3-'---5-' |<-----| 3-'---5-'
#     # p_n_count : 5+'-- 3+' |----->| 5-'---3-' or  3-'---5-' |<-----| 3+'-- 5+'
#     # p_p_count : one read map two edges
#     # p_n_count : one read map two edges with different direction




def count_reads_through_edges(graphs, bam_file, chrom):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for G in graphs:
        for (node1_id, node2_id) in G.edges:
            node1_data = G.nodes[node1_id]
            node2_data = G.nodes[node2_id]
            if not node1_data.get('source', False) and not node2_data.get('sink', False):
                print("Processing edge", node1_data, node2_data)
                map_cnt, all_cnt = count_and_print_reads_in_interval(samfile, chrom, node1_data['end'], node2_data['start'])
                G.edges[node1_id, node2_id]['reads'] = map_cnt
                G.edges[node1_id, node2_id]['capacity'] = map_cnt

                # Add the number of reads to the node's 'reads' attribute
                G.nodes[node2_id]['reads'] += map_cnt

        # Find the source and sink nodes
        source_node = [node for node, data in G.nodes(data=True) if data.get('source', False)][0]
        sink_node = [node for node, data in G.nodes(data=True) if data.get('sink', False)][0]

        # Handle case where there is only one node between source and sink
        if G.number_of_nodes() == 3:  
            single_node_id = [node for node in G.nodes if node != source_node and node != sink_node][0]
            G.edges[source_node, single_node_id]['capacity'] = G.nodes[single_node_id]['reads']
            G.edges[single_node_id, sink_node]['capacity'] = G.nodes[single_node_id]['reads']
        else:
            # Calculate the capacity for edges from source node
            for _, target_node in G.edges(source_node):
                # Sum the capacity of all outgoing edges from the target_node
                G.edges[source_node, target_node]['capacity'] = sum([G.edges[target_node, node].get('capacity', 0) for node in G.successors(target_node)])

            # Calculate the capacity for edges to sink node
            for source_node, _ in G.in_edges(sink_node):
                # Sum the capacity of all outgoing edges from the source_node
                G.edges[source_node, sink_node]['capacity'] = sum([G.edges[node, source_node].get('capacity', 0) for node in G.predecessors(source_node)])
    samfile.close()

def calculate_max_flows(graphs):
    max_flows = []
    max_flow_graphs = []
    for G in graphs:
        # Find the source and sink nodes
        source_node = [node for node, data in G.nodes(data=True) if data.get('source', False)][0]
        sink_node = [node for node, data in G.nodes(data=True) if data.get('sink', False)][0]

        # Calculate the maximum flow from the source to the sink
        flow_value, flow_dict = nx.maximum_flow(G, source_node, sink_node, capacity='capacity')

        # Store the maximum flow value and dict for later use
        max_flows.append((flow_value, flow_dict))

        # Create a new graph from the flow dict
        max_flow_graph = nx.DiGraph()
        for node, edges in flow_dict.items():
            for edge, flow in edges.items():
                if flow > 0:
                    max_flow_graph.add_edge(node, edge, flow=flow)
        max_flow_graphs.append(max_flow_graph)
    return max_flows, max_flow_graphs


def draw_graphs(graphs, prefix='graph'):
    for i, G in enumerate(graphs):
        print("Processing graph", i)
        for node in graphs[i].nodes:
            print(graphs[i].nodes[node])
        for edge in graphs[i].edges:
            print(f"Edge from {edge[0]} to {edge[1]} with attributes {graphs[i].edges[edge]}")

        plt.figure(figsize=(8, 6))

        # Find the maximum capacity to scale the widths and colors
        capacities = nx.get_edge_attributes(G, 'capacity').values()
        max_capacity = max(capacities)

        # Calculate the layout once
        pos = nx.shell_layout(G)

        # Draw edges with width proportional to capacity, color representing capacity and dashed style
        for (u, v, capacity) in G.edges.data('capacity'):
            # Make sure the color is not too light
            color = plt.cm.Blues(max(0.5, capacity/max_capacity))

            # Create a FancyArrowPatch, which allows us to control the shape of the arrow
            arrow = mpatches.FancyArrowPatch(pos[u], pos[v], arrowstyle='->', mutation_scale=30,
                                             connectionstyle='arc3', color=color, linewidth=2*capacity/max_capacity)
            plt.gca().add_patch(arrow)

        # Draw nodes
        nx.draw_networkx_nodes(G, pos=pos)

        # Draw node labels
        nx.draw_networkx_labels(G, pos=pos)

        plt.savefig(f'{prefix}_{i}.png')
        plt.close()

import matplotlib.pyplot as plt

# Function to find max flow path
def find_max_flow_path(flow_graph, source, sink):
    max_flow = 0
    max_flow_path = None
    for path in nx.all_simple_paths(flow_graph, source, sink):
        path_flow = min(flow_graph[u][v]['capacity'] for u, v in zip(path[:-1], path[1:]))
        if path_flow > max_flow:
            max_flow = path_flow
            max_flow_path = path
    return max_flow_path, max_flow

def draw_graphs2(graphs, prefix='graph'):
    for i, graph in enumerate(graphs):

        # print graph info
        

        # Source and target nodes
        source = [node for node, data in graph.nodes(data=True) if data.get('source', False)][0]
        sink = [node for node, data in graph.nodes(data=True) if data.get('sink', False)][0]

        iteration = 0
        while True:
            # Create a new figure for each iteration
            print("Processing iter'", iteration)
            print("##################")
            for node in graphs[i].nodes:
                print(graphs[i].nodes[node])
            for edge in graphs[i].edges:
                print(f"Edge from {edge[0]} to {edge[1]} with attributes {graphs[i].edges[edge]}")
            plt.figure(figsize=(8, 12))

            # Compute the maximum flow network
            flow_value, flow_dict = nx.maximum_flow(graph, source, sink, capacity='capacity')

            # If no more augmenting paths are found, break the loop
            if flow_value == 0:
                break

            # Draw the original graph
            plt.subplot(3, 1, 1)
            pos = nx.spring_layout(graph)  # positions for all nodes
            nx.draw(graph, pos, with_labels=True, node_color='skyblue', arrows=True)
            edge_labels = nx.get_edge_attributes(graph, 'capacity')
            nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
            plt.title(f'Graph {i+1}, Iteration {iteration+1} O riginal Graph')

            # Draw the graph with maximum flow
            plt.subplot(3, 1, 2)
            max_flow_graph = nx.DiGraph(flow_dict)
            edge_labels_flow = nx.get_edge_attributes(max_flow_graph, 'capacity')
            nx.draw(max_flow_graph, pos, with_labels=True, node_color='skyblue', arrows=True)
            nx.draw_networkx_edge_labels(max_flow_graph, pos, edge_labels=edge_labels)

            plt.title(f'Graph {i+1}, Iteration {iteration+1} with Maximum Flow (Flow={flow_value})')

            # Draw the graph with the path of maximum flow at each iteration
            plt.subplot(3, 1, 3)
            max_flow_path, max_flow = find_max_flow_path(graph, source, sink)
            edge_labels_path = {(u, v): d['capacity'] for u, v, d in graph.edges(data=True) if (u,v) in zip(max_flow_path, max_flow_path[1:])}
            J = graph.edge_subgraph(edge_labels_path.keys())
            pos_line = {node: (idx, 1) for idx, node in enumerate(max_flow_path)}  # positions in a line
            nx.draw(J, pos_line, with_labels=True, node_color='skyblue', edge_color="red", arrows=True)
            nx.draw_networkx_edge_labels(J, pos_line, edge_labels=edge_labels_path)
            plt.title(f'Graph {i+1}, Iteration {iteration+1} with Maximum Flow Path (Flow={max_flow})')

            # Update the 'capacity' attribute of the edges in the max-flow path
            for u, v in zip(max_flow_path, max_flow_path[1:]):
                if graph.has_edge(u, v):
                    graph[u][v]['capacity'] -= max_flow
            
            print("Left flow:", flow_value)

            iteration += 1

            # Save the figure
            plt.savefig(f'{prefix}_{i}_iteration_{iteration}.png')

            # Display the figure
            # plt.show()
        break
        # plt.show()

def main():
    parser = argparse.ArgumentParser(description='Parse a graph file.')
    parser.add_argument('file_name', type=str, help='The name of the file to parse.')
    parser.add_argument('bam_file', type=str, help='Path to BAM file.')
    parser.add_argument('chrom', type=str, help='Chromosome name.')
    args = parser.parse_args()

    # Use the function
    graphs = parse_graphs(args.file_name)
    print(f'Parsed {len(graphs)} graphs from {args.file_name}.')
    count_reads_through_edges(graphs, args.bam_file, args.chrom)
    # draw_graphs2(graphs)
    # print(graphs[0].nodes, len(graphs[0].nodes))
    # print(graphs[0].edges, len(graphs[0].edges))
    # for node in graphs[0].nodes:
    #     print(graphs[0].nodes[node])
    # for edge in graphs[0].edges:
    #     print(f"Edge from {edge[0]} to {edge[1]} with attributes {graphs[0].edges[edge]}")

if __name__ == '__main__':
    main()