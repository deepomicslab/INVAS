import argparse
from collections import defaultdict
def main():
    all_nodes=[]
    trf_recs=[]
    with open(args.trf_in, 'r') as f:
        for line in f:
            line = line.strip()
            strand, paths, cov = line.split(' ')
            trf_recs.append((strand, paths, cov))
            paths=paths.split(',')[:-1]
            for p in paths:
                all_nodes.append(int(p))
    all_nodes=list(set(all_nodes))
    all_nodes.sort()
    print(all_nodes)
    source_node = all_nodes[-2]
    sink_node = all_nodes[-1]
    first_node_dict={}
    last_node_dict={}
    for trf in trf_recs:
        if trf[1].strip(",").startswith(str(source_node)):
            paths=trf[1].strip(",").split(',')
            first_node=paths[1]
            first_node_dict[first_node]=True
        if trf[1].strip(",").endswith(str(sink_node)):
            paths=trf[1].strip(",").split(',')
            last_node=paths[-2]
            last_node_dict[last_node]=True
    remain_first_node_dict=defaultdict(float)
    remain_last_node_dict=defaultdict(float)
    print(first_node_dict)
    print(last_node_dict)
    for trf in trf_recs:
        if not trf[1].strip(",").startswith(str(source_node)):
            paths=trf[1].strip(",").split(',')
            first_node=paths[0]
            if first_node in first_node_dict:
                remain_first_node_dict[first_node]+=float(trf[2])
        if not trf[1].strip(",").endswith(str(sink_node)):
            print("here not add end", trf[1].strip(","))
            paths=trf[1].strip(",").split(',')
            last_node=paths[-1]
            print(last_node )
            if last_node in last_node_dict:
                print("find last node", last_node)
                remain_last_node_dict[last_node]+=float(trf[2])
    print(remain_first_node_dict)
    print(remain_last_node_dict)

    
    
    with open(args.trf_out, 'w') as f:
        for trf in trf_recs:
            f.write(f"{trf[0]} {trf[1]} {trf[2]}\n")
        for k, v in remain_first_node_dict.items():
            f.write(f"{1} {source_node},{k}, {v}\n")
        for k, v in remain_last_node_dict.items():
            f.write(f"{1} {k},{sink_node}, {v}\n")
            
    
    
    
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse transcript exons and BAM status files.")
    parser.add_argument("trf_in", help="File containing transcript exons information.")
    parser.add_argument("trf_out", help="File containing BAM status information.")


    args = parser.parse_args()
    main()