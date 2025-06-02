import csv
import sys
from collections import defaultdict

def merge_data(input_file_path, output_file_path):
    # 使用defaultdict来存储合并的数据
    merged_data = defaultdict(float)
    all_nodes= []
    # 读取输入文件
    with open(input_file_path, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue  # 跳过格式不正确的行
            key = row[0].strip() + " " +row[1].strip()  # 第二列作为键，去除可能的空白字符
            paths = row[1].strip(',').split(',')
            for p in paths:
                if p not in all_nodes:
                    all_nodes.append(p)
            print(key)
            value = float(row[2].strip())  # 第三列转换为浮点数
            merged_data[key] += value  # 合并第三列的值

    print(merged_data)
    # sort all_nodes
    all_nodes.sort()
    print(all_nodes)

    # 将合并后的数据写入输出文件
    with open(output_file_path, 'w',) as outfile:
        for k,v in merged_data.items():
            outfile.write(f"{k} {v}\n")
        

    print(f"Merged data has been written to {output_file_path}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    merge_data(input_file, output_file)