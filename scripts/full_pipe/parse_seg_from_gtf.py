import csv

def parse_stringtie_gtf(gtf_file, output_file):
    """
    解析 StringTie 生成的 GTF 文件，提取 chr, exon_start, exon_end, gene_id, strand 信息。
    """
    with open(gtf_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        # 写入输出文件的标题行
        writer.writerow(['chr', 'exon_start', 'exon_end', 'gene_id', 'strand'])
        
        for line in infile:
            # 跳过注释行（以 # 开头）
            if line.startswith('#'):
                continue
            
            # 将 GTF 行分割为字段
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # 确保行包含至少 9 个字段
            
            # 提取字段
            chrom = fields[0]          # 染色体编号
            feature_type = fields[2]   # 功能类型（例如 exon, transcript 等）
            start = fields[3]          # 起始位置
            end = fields[4]            # 结束位置
            strand = fields[6]         # 链（+ 或 -）
            attributes = fields[8]     # 属性字段
            
            # 过滤只提取 "exon" 类型的行
            if feature_type != 'exon':
                continue
            
            # 从属性字段中提取 gene_id
            gene_id = None
            for attribute in attributes.split(';'):
                attribute = attribute.strip()
                if attribute.startswith('gene_id'):
                    # 提取 gene_id 的值（去掉引号）
                    gene_id = attribute.split(' ')[1].replace('"', '')
                    break
            
            # 确保 gene_id 存在
            if gene_id is None:
                continue
            
            # 写入提取的数据到输出文件
            writer.writerow([chrom, start, end, gene_id, strand])

    print(f"提取完成，结果已写入：{output_file}")


# 示例调用
# 输入 GTF 文件路径
gtf_file = 'stringtie/stringtie.app.a2.gtf'
# 输出文件路径
output_file = 'exons_output.tsv'

# 调用函数解析 GTF 文件
parse_stringtie_gtf(gtf_file, output_file)