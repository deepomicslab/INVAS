#!/bin/bash
# mapping_evaluation.sh
#
# Usage:
#   ./mapping_evaluation.sh <sample_list.txt>
#
# Description:
#   对于样本列表文件中的每个样本：
#     1. 从 ../<sample_name>/truth/inversions.tsv 文件中去除 header 和重复值，得到 unique gene list；
#     2. 对于每个 gene，获取对应的 paired-end reads：
#            ../<sample_name>/rna/genes/<gene>_1.fq 与 ../<sample_name>/rna/genes/<gene>_2.fq
#     3. 分别检测三个组装软件的结果：
#            - stringtie: assembly 文件位于 ../<sample_name>/rna/genes/<gene>/analysis/stringtie/<gene>.stringtie.fasta
#            - cufflinks: assembly 文件位于 ../<sample_name>/rna/genes/<gene>/analysis/cufflinks/<gene>.cufflinks.fasta（以 transcripts.gtf 判断结果存在）
#            - trinity:   assembly 文件位于 ../<sample_name>/rna/genes/<gene>/analysis/trinity/Trinity.fasta
#     4. 如果对应结果存在，则调用 zsh eva_map_single3.sh 进行 mapping 评价，
#        调用命令格式为：
#            zsh eva_map_single3.sh <fq1> <fq2> <assembly.fasta> <output_folder>
#        其中 <output_folder> 为保存 mapping 结果的输出目录（自动创建）。
#
# Note:
#   请确保 zsh eva_map_single3.sh 脚本已赋予可执行权限，并且其参数顺序与上面描述一致。
#
# 设置严格模式：一旦命令失败则退出，并禁止使用未定义的变量
set -euo pipefail

# 检查样本列表文件是否提供
if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample_list.txt>"
    exit 1
fi

# 输入的样本列表文件，每行格式为: sample_name wgs_sample1 rna_sample
filename="$1"

while read -r sample_name wgs_sample1 rna_sample; do
    echo "Processing sample '$sample_name' with WGS sample '$wgs_sample1' and RNA sample '$rna_sample'..."
    
    # 定位 gene 信息文件，并判断是否存在
    gene_file="../${sample_name}/truth/inversions.tsv"
    # if [ ! -f "$gene_file" ]; then
    #     echo "  Gene file '$gene_file' not found for sample '$sample_name'. Skipping this sample."
    #     continue
    # fi

    # # 去除 header 并提取第一列生成 unique gene list
    # awk 'NR>1 {print $1}' "$gene_file" | sort | uniq > "${gene_file}.unique"
    
    # 遍历 gene 列表中的每个 gene
    while read -r gene; do
        echo "  Processing gene '$gene'..."
        
        # 定义 gene 对应的 paired-end reads 文件
        fq1="../${sample_name}/rna/genes/${gene}/${gene}_1.fq"
        fq2="../${sample_name}/rna/genes/${gene}/${gene}_2.fq"

        if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
            echo "    Reads file not found for gene '$gene': ($fq1 or $fq2). Skipping mapping evaluation for this gene."
            continue
        fi

        # ----------------- stringtie 评估 -----------------
        st_dir="../${sample_name}/rna/genes/${gene}/analysis/stringtie"
        st_asm="${st_dir}/${gene}.stringtie.fasta"
        if [ -d "$st_dir" ] && [ -f "$st_asm" ]; then
            echo "    Running mapping evaluation for gene '$gene' using stringtie assembly..."
            out_dir="../${sample_name}/eva/${gene}/analysis/stringtie/seq_map"
            mkdir -p "$out_dir"
            # if out_dir/mapping_metrics.tsv exists, skip
            if [ -f "$out_dir/mapping_metrics.tsv" ]; then
                echo "    Mapping evaluation already done for gene '$gene' using stringtie assembly. Skipping."
                continue
            fi
            # 调用 zsh eva_map_single3.sh，参数顺序：<fq1> <fq2> <assembly.fasta> <output_folder>
            zsh eva_map_single3.sh "$fq1" "$fq2" "$st_asm" "$out_dir"
        else
            echo "    Skipping gene '$gene' in stringtie: assembly data not found."
        fi

        # ----------------- cufflinks 评估 -----------------
        cl_dir="../${sample_name}/rna/genes/${gene}/analysis/cufflinks"
        cl_gtf="${cl_dir}/transcripts.gtf"
        cl_asm="${cl_dir}/${gene}.cufflinks.fasta"
        if [ -d "$cl_dir" ] && [ -f "$cl_gtf" ] && [ -f "$cl_asm" ]; then
            echo "    Running mapping evaluation for gene '$gene' using cufflinks assembly..."
            out_dir="../${sample_name}/eva/${gene}/analysis/cufflinks/seq_map"
            mkdir -p "$out_dir"
            if [ -f "$out_dir/mapping_metrics.tsv" ]; then
                echo "    Mapping evaluation already done for gene '$gene' using cufflinks assembly. Skipping."
                continue
            fi
            zsh eva_map_single3.sh "$fq1" "$fq2" "$cl_asm" "$out_dir"
        else
            echo "    Skipping gene '$gene' in cufflinks: assembly data not found."
        fi

        # ----------------- trinity 评估 -----------------
        tr_dir="../${sample_name}/rna/genes/${gene}/analysis/trinity"
        tr_asm="${tr_dir}/Trinity.fasta"
        if [ -d "$tr_dir" ] && [ -f "$tr_asm" ]; then
            echo "    Running mapping evaluation for gene '$gene' using trinity assembly..."
            out_dir="../${sample_name}/eva/${gene}/analysis/trinity/seq_map"
            mkdir -p "$out_dir"
            if [ -f "$out_dir/mapping_metrics.tsv" ]; then
                echo "    Mapping evaluation already done for gene '$gene' using trinity assembly. Skipping."
                continue
            fi
            zsh eva_map_single3.sh "$fq1" "$fq2" "$tr_asm" "$out_dir"
        else
            echo "    Skipping gene '$gene' in trinity: assembly data not found."
        fi

    done < "${gene_file}.unique"

done < "$filename"