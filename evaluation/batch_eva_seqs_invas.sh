#!/bin/bash
# eva_invas_haps.sh
#
# Usage:
#   ./eva_invas_haps.sh <sample_list.txt>
#
# Description:
#   对于样本列表中的每个 sample（每行包含 sample_name、wgs_sample1、rna_sample）：
#     1. 从 ../${sample_name}/truth/inversions.tsv 中去掉 header 并提取第一列，得到唯一的 gene_id 列表 (inversions.tsv.unique)。
#     2. 对于每个 gene_id，通过 all_truth_inversions.tsv 将 gene_id 转换为 gene_name。
#     3. 构造 haps 组装转录本文件路径：
#            ../${sample_name}/haps/${gene_name}/haps2/transcripts_haps.fa
#        如果该文件存在，则调用 eva_seqs_single2.py 进行评估，参数如下：
#
#             --ref ../${sample_name}/truth_seq/${gene}/transcripts.fa
#             --asm ../${sample_name}/haps/${gene_name}/haps2/transcripts_haps.fa
#             --out ../${sample_name}/eva/${gene}/analysis/invas/seq_eva2
#             --threads 30
#
# Note:
#   - 其中 gene 变量代表 gene_id（用于 reads、truth 文件以及输出目录的命名），
#     而组装文件路径则需要转换为 gene_name 后再构造。
#   - 请确保 all_truth_inversions.tsv 与 eva_seqs_single2.py 脚本可在运行目录下正确访问。

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    exit 1
fi

sample_list="$1"

while read -r sample_name wgs_sample1 rna_sample; do
    echo "Processing sample '$sample_name' with WGS sample '$wgs_sample1' and RNA sample '$rna_sample'..."
    
    gene_file="../${sample_name}/truth/inversions.tsv"
    if [ ! -f "$gene_file" ]; then
        echo "  Gene file '$gene_file' not found for sample '$sample_name'. Skipping sample."
        continue
    fi

    # 去除 header 并提取第一列生成唯一的 gene_id 列表（保存为 inversions.tsv.unique）
    awk 'NR>1 {print $1}' "$gene_file" | sort -u > "${gene_file}.unique"

    # 对每个 gene_id 进行处理
    while read -r gene; do
        echo "  Processing gene_id '$gene'..."
        
        # 通过 all_truth_inversions.tsv 将 gene_id 转换为 gene_name
        gene_name=$(awk -F'\t' -v g="$gene" '$1 == g {print $2; exit}' all_truth_inversions.tsv)
        if [ -z "$gene_name" ]; then
            echo "    Unable to convert gene_id '$gene' to gene_name using all_truth_inversions.tsv. Skipping."
            continue
        fi
        echo "    gene_id '$gene' converts to gene_name '$gene_name'."
        
        # 定义 haps 组装转录本文件路径（使用转换后的 gene_name）
        asm_file="../${sample_name}/haps/${gene_name}/haps2/transcripts_haps.fa"
        if [ ! -f "$asm_file" ]; then
            echo "    Assembly file '$asm_file' not found. Skipping gene '$gene'."
            continue
        fi
        echo "    Found assembly file '$asm_file'. Running evaluation..."
        
        # 创建输出文件夹，固定为 ../${sample_name}/eva/${gene}/analysis/invas/seq_eva2 （使用 gene_id）
        out_dir="../${sample_name}/eva/${gene}/analysis/invas/seq_eva2"
        mkdir -p "$out_dir"
        
        # 调用 eva_seqs_single2.py 进行评估
        # if "$out_dir"/mapping_metrics.tsv exists, skip
        if [ -f "$out_dir/mapping_metrics.tsv" ]; then
            echo "    Evaluation results already exist. Skipping."
            continue
        fi
        python eva_seqs_single2.py \
            --ref "../${sample_name}/truth_seq/${gene}/transcripts.fa" \
            --asm "$asm_file" \
            --out "$out_dir" \
            --threads 30

    done < "${gene_file}.unique"

done < "$sample_list"