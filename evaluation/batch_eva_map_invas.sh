#!/bin/bash
# eva_map_invas.sh
#
# Usage:
#   ./eva_map_invas.sh <sample_list.txt>
#
# Description:
#   对于样本列表中的每个 sample，
#     1. 从 ../<sample_name>/truth/inversions.tsv 中去除 header 并提取第一列得到 unique 的 gene_id 列表。
#     2. 利用 all_truth_inversions.tsv（第一列：gene_id，第二列：gene_name）将 gene_id 转换为 gene_name。
#     3. 检查 haps 组装转录本文件:
#            ../<sample_name>/haps/<gene_name>/haps2/transcripts_haps.fa
#        是否存在。
#     4. 如果存在，则使用对应的 paired-end reads 文件:
#            ../<sample_name>/rna/genes/<gene>_1.fq 与 ../<sample_name>/rna/genes/<gene>_2.fq
#        调用 zsh eva_map_single3.sh  进行 mapping 评价，参数顺序为：
#           <fq1> <fq2> <assembly.fasta> <output_folder>
#        输出文件夹固定为:
#           ../<sample_name>/eva/<gene_id>/analysis/invas/seq_map
#
# Note:
#   请确保 all_truth_inversions.tsv 文件与 zsh eva_map_single3.sh  脚本在运行目录下可正常访问。
#   如果在 all_truth_inversions.tsv 中无法通过 gene_id 转换到 gene_name，则跳过该 gene 的评估。

set -euo pipefail

# 检查输入参数
if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    exit 1
fi

filename="$1"

while read -r sample_name wgs_sample1 rna_sample; do
    echo "Processing sample '$sample_name' with WGS sample '$wgs_sample1' and RNA sample '$rna_sample'..."
    
    gene_file="../${sample_name}/truth/inversions.tsv"
    if [ ! -f "$gene_file" ]; then
        echo "  Gene file '$gene_file' not found for sample '$sample_name'. Skipping sample."
        continue
    fi

    # 去除 header 并提取第一列生成 unique 的 gene_id 列表
    awk 'NR>1 {print $1}' "$gene_file" | sort | uniq > "${gene_file}.unique"

    # 遍历每个 gene_id
    while read -r gene; do
        echo "  Processing gene_id '$gene'..."
        
        # 利用 all_truth_inversions.tsv 将 gene_id 转换为 gene_name
        gene_name=$(awk -F'\t' -v g="$gene" '$1 == g {print $2; exit}' all_truth_inversions.tsv)
        if [ -z "$gene_name" ]; then
            echo "    Cannot convert gene_id '$gene' to gene_name using all_truth_inversions.tsv. Skipping."
            continue
        fi
        echo "    gene_id '$gene' converts to gene_name '$gene_name'."
        
        # 定义 paired-end reads 文件（以 gene_id 为依据）
        fq1="../${sample_name}/rna/genes/${gene}/${gene}_1.fq"
        fq2="../${sample_name}/rna/genes/${gene}/${gene}_2.fq"
        if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
            echo "    Reads file(s) not found for gene '$gene': '$fq1' or '$fq2'. Skipping mapping evaluation."
            continue
        fi

        # 定义 haps 组装转录本文件路径，使用 gene_name 构造路径
        haps_asm="../${sample_name}/haps/${gene_name}/haps2/transcripts_haps.fa"
        if [ -f "$haps_asm" ]; then
            echo "    Found haps assembly for gene_name '$gene_name' at '$haps_asm'."
            # 输出文件夹固定为 "../<sample_name>/eva/<gene_id>/analysis/invas/seq_map"
            out_dir="../${sample_name}/eva/${gene}/analysis/invas/seq_map"
            # if out_dir/mapping_metrics.tsv exists, skip

            mkdir -p "$out_dir"
            if [ -f "${out_dir}/mapping_metrics.tsv" ]; then
                echo "    Mapping evaluation already exists for gene_id '$gene'. Skipping."
                continue
            fi
            # 调用 zsh eva_map_single3.sh ，参数顺序：<fq1> <fq2> <assembly.fasta> <output_folder>
            zsh eva_map_single3.sh  "$fq1" "$fq2" "$haps_asm" "$out_dir"
        else
            echo "    Haps assembly file '$haps_asm' not found. Skipping evaluation for gene_id '$gene'."
        fi

    done < "${gene_file}.unique"

done < "$filename"