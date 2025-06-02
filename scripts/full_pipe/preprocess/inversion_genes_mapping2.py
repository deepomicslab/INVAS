# inversion_genes_mapping.py
import argparse

def map_inversions_to_genes(input_file, output_file):
    inversion_gene_map = {}

    with open(input_file, 'r') as infile:
        for line in infile:
            fields = line.strip().split('\t')
            inversion_key = f"{fields[0]}:{fields[1]}-{fields[2]}"
            gene_name = fields[6]
            gene_start = int(fields[4])
            gene_end = int(fields[5])
            gene_length = gene_end - gene_start

            # 将基因信息存储为元组并添加到字典中
            gene_info = (gene_name, gene_start, gene_end, gene_length)

            if inversion_key in inversion_gene_map:
                inversion_gene_map[inversion_key].append(gene_info)
            else:
                inversion_gene_map[inversion_key] = [gene_info]

    with open(output_file, 'w') as outfile:
        for inversion, genes in inversion_gene_map.items():
            # 根据基因的长度降序排序
            genes_sorted = sorted(genes, key=lambda x: x[3], reverse=True)
            
            # 获取长度最大的基因及其位置，格式为：基因名:起始位置-终止位置
            largest_gene = genes_sorted[0]
            largest_gene_str = f"{largest_gene[0]}:{largest_gene[1]}-{largest_gene[2]}"
            
            # 将剩余基因名称合并为一个字符串
            other_genes_str = ','.join([gene[0] for gene in genes_sorted[1:]])
            
            # 如果还有其他基因，则将它们加入到字符串中
            genes_str = largest_gene_str if not other_genes_str else f"{largest_gene_str},{other_genes_str}"
            
            outfile.write(f"{inversion}\t{genes_str}\n")

def main():
    parser = argparse.ArgumentParser(description="Map inversions to genes.")
    parser.add_argument("input_file", help="Input file path containing inversions")
    parser.add_argument("output_file", help="Output file path to write the mappings")
    
    args = parser.parse_args()

    map_inversions_to_genes(args.input_file, args.output_file)

if __name__ == "__main__":
    main()