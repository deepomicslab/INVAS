import csv
import sys
from collections import defaultdict

def merge_data(input_file_path, output_file_path):
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
    
    with open(output_file_path, 'w') as f:
        for k, v in new_paths.items():
            f.write(f"{strand} {k} {v}\n")

            





        

    print(f"Merged data has been written to {output_file_path}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    merge_data(input_file, output_file)