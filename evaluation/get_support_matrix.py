import os
import pandas as pd
import sys

def parse_support_file(file_path):
    """解析单个support.txt文件，返回基因-支持数的字典"""
    gene_support = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                # 获取基因名称 (格式如 NOL10:10710891-10830101 中提取 NOL10)
                gene = parts[1].split(':')[0]
                support_count = int(parts[-1])
                gene_support[gene] = support_count
    return gene_support

def main():
    res_dir = "../res_plus/"
    sample_dir = sys.argv[1]
    all_support = {}
    all_genes = set()
    
    # 获取所有样本目录
    # samples = [d for d in os.listdir(res_dir) 
    #           if os.path.isdir(os.path.join(res_dir, d))]
    # print(f"共有 {len(samples)} 个样本")
    # print(f"样本列表: {samples}")
    samples=[sample_dir]
    
    # 处理每个样本的support.txt
    for sample in samples:
        support_file = os.path.join(res_dir, sample, "support.txt")
        if os.path.exists(support_file):
            gene_support = parse_support_file(support_file)
            all_support[sample] = gene_support
            all_genes.update(gene_support.keys())
            print(f"样本 {sample} 的支持数已加载")
            print(f"共有 {len(gene_support)} 个基因")
    
    # 创建数据框
    df = pd.DataFrame(0, index=sorted(list(all_genes)), 
                     columns=sorted(all_support.keys()))
    
    # 填充数据
    for sample in all_support:
        for gene, support in all_support[sample].items():
            df.loc[gene, sample] = support
    
    # 保存到文件
    output_file = os.path.join(res_dir, sample_dir,"gene_support_matrix.txt")
    df.to_csv(output_file, sep='\t')
    print(f"结果已保存到: {output_file}")

if __name__ == "__main__":
    main()
