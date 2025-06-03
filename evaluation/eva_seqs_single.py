#!/usr/bin/env python3
"""
Script: compare_rna_alignment_with_metrics.py
Description:
    使用 minimap2 将组装的 RNA 转录本比对到真实 RNA 转录本上，
    解析 PAF 输出，并统计：
      1. 真实转录本回收率 (Reference Transcript Recovery Rate, RTRR)
      2. 高质量组装转录本比例 (High-Quality Assembly Ratio, HQAR)
         （高质量定义为比对时查询覆盖率 (qcov) ≥ 80%，目标覆盖率 (tcov) ≥ 80%，且百分比相似性 (pid) ≥ 90%）
    
Usage:
    python compare_rna_alignment_with_metrics.py --ref reference.fasta --asm assembled.fasta --outdir results

Requirements:
    - Python 3
    - minimap2 (已加入系统 PATH)
    
Author: ChatGPT
Date: 2025-03-17
"""

import argparse
import subprocess
import os
import sys

def run_minimap2(query_fasta, target_fasta, out_file, preset="asm5", threads=4):
    """
    调用 minimap2 进行比对。
    
    Parameters:
        query_fasta: FASTA 文件，作为查询（组装转录本）
        target_fasta: FASTA 文件，作为目标（真实转录本）
        out_file: PAF 格式输出文件
        preset: minimap2 预设参数（例如 asm5 或 asm10，根据序列相似性选择）
        threads: 使用的线程数
    Returns:
        成功返回 True，否则返回 False。
    """
    cmd = [
        "minimap2",
        "-cx", preset,
        "-t", str(threads),
        target_fasta,
        query_fasta
    ]
    print("Running command:", " ".join(cmd))
    try:
        with open(out_file, "w") as fout:
            subprocess.run(cmd, stdout=fout, check=True)
    except subprocess.CalledProcessError as e:
        print("Error running minimap2:", e)
        return False
    return True

def parse_paf(paf_file):
    """
    解析 PAF 文件，对每个查询（组装转录本）保留比对区域最长的记录。
    
    每条记录中提取的信息包括：
      - qname: 组装转录本 ID
      - qlen: 组装转录本全长
      - qstart, qend: 比对在组装转录本中的起始与终止位置
      - tname: 真实转录本 ID
      - tlen, tstart, tend: 真实转录本长度及比对区间位置
      - matches: 比对匹配碱基数
      - aln_len: 比对区域长度
      - mapq: 映射质量
      - qcov: 查询覆盖率 = (qend - qstart) / qlen
      - tcov: 目标覆盖率 = (tend - tstart) / tlen
      - pid: 百分比相似性 approximated as matches / aln_len * 100
     
    Returns:
        dict {query_id: record}
    """
    best_alignments = {}
    with open(paf_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 12:
                continue
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tlen = int(fields[6])
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])
            qcov = (qend - qstart) / qlen if qlen > 0 else 0.0
            tcov = (tend - tstart) / tlen if tlen > 0 else 0.0
            pid = (matches / aln_len * 100) if aln_len > 0 else 0.0

            record = {
                "qname": qname,
                "qlen": qlen,
                "qstart": qstart,
                "qend": qend,
                "strand": strand,
                "tname": tname,
                "tlen": tlen,
                "tstart": tstart,
                "tend": tend,
                "matches": matches,
                "aln_len": aln_len,
                "mapq": mapq,
                "qcov": qcov,
                "tcov": tcov,
                "pid": pid
            }
            # 针对同一查询，保留比对区域最长的记录
            if qname in best_alignments:
                if aln_len > best_alignments[qname]["aln_len"]:
                    best_alignments[qname] = record
            else:
                best_alignments[qname] = record
    return best_alignments

def read_fasta_ids(fasta_file):
    """
    解析 FASTA 文件，返回所有序列的 ID 集合。
    """
    seq_ids = set()
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                seq_ids.add(header)
    return seq_ids

def main():
    parser = argparse.ArgumentParser(
        description=("使用 minimap2 比对组装 RNA 至真实 RNA，并计算两个关键指标：\n"
                     "1. 真实转录本回收率 (RTRR)。\n"
                     "2. 高质量组装转录本比例 (HQAR)，其中高质量定义为 qcov>=80%，tcov>=80%，pid>=90%。")
    )
    parser.add_argument("--ref", required=True, help="真实 RNA 转录本 FASTA 文件")
    parser.add_argument("--asm", required=True, help="组装 RNA 转录本 FASTA 文件")
    parser.add_argument("--outdir", default="results", help="结果输出目录")
    parser.add_argument("--threads", type=int, default=4, help="minimap2 使用的线程数")
    parser.add_argument("--preset", default="asm5", help="minimap2 预设参数（例如 asm5 或 asm10）")
    # 定义高质量标准阈值
    parser.add_argument("--min_qcov", type=float, default=0.7, help="高质量比对的最小查询覆盖率（0-1）")
    parser.add_argument("--min_tcov", type=float, default=0.7, help="高质量比对的最小目标覆盖率（0-1）")
    parser.add_argument("--min_pid", type=float, default=80.0, help="高质量比对的最小百分比相似性（%）")
    
    args = parser.parse_args()

    # 创建输出目录
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # 统计输入 FASTA 文件中的转录本数量
    print("Count input sequences...")
    asm_ids = read_fasta_ids(args.asm)
    ref_ids = read_fasta_ids(args.ref)
    total_asm = len(asm_ids)
    total_ref = len(ref_ids)
    print(f"Total assembled transcripts: {total_asm}")
    print(f"Total reference transcripts: {total_ref}")

    # 用 minimap2 将组装 RNA 比对到真实 RNA
    paf_file = os.path.join(args.outdir, "asm2ref.paf")
    print("\nRunning minimap2 alignment of assembled RNA to reference RNA...")
    if not run_minimap2(args.asm, args.ref, paf_file, preset=args.preset, threads=args.threads):
        sys.exit(1)

    # 解析 PAF 文件
    print("\nParsing minimap2 PAF output...")
    alignments = parse_paf(paf_file)
    aligned_asm_ids = set(alignments.keys())
    num_aligned = len(aligned_asm_ids)
    missing_asm = total_asm - num_aligned

    # 取出比对记录中所有被 hit 的真实转录本 ID
    hit_ref_ids = set(rec["tname"] for rec in alignments.values())
    missing_ref = total_ref - len(hit_ref_ids)
    # 计算真实转录本回收率 (RTRR)
    rtrr = (len(hit_ref_ids) / total_ref * 100) if total_ref > 0 else 0.0

    # 统计符合高质量标准的组装转录本数
    high_quality_count = 0
    for rec in alignments.values():
        if rec["qcov"] >= args.min_qcov and rec["tcov"] >= args.min_tcov and rec["pid"] >= args.min_pid:
            high_quality_count += 1
    # 可选：可以分别计算相对于所有组装序列和仅相对于比对上的 assembled 的比例
    hq_ratio_total = (high_quality_count / total_asm * 100) if total_asm > 0 else 0.0
    hq_ratio_aligned = (high_quality_count / num_aligned * 100) if num_aligned > 0 else 0.0

    # 计算比对质量的平均指标（仅针对比对上的记录）
    total_qcov = sum(rec["qcov"] for rec in alignments.values())
    total_tcov = sum(rec["tcov"] for rec in alignments.values())
    total_pid  = sum(rec["pid"] for rec in alignments.values())
    avg_qcov = (total_qcov / num_aligned) * 100 if num_aligned > 0 else 0.0
    avg_tcov = (total_tcov / num_aligned) * 100 if num_aligned > 0 else 0.0
    avg_pid  = total_pid / num_aligned if num_aligned > 0 else 0.0

    # 输出详细比对报告（TSV格式）
    report_file = os.path.join(args.outdir, "alignment_report.tsv")
    with open(report_file, "w") as outfh:
        header = [
            "Assembled_ID", "Reference_ID", "Assembled_Length", "Query_Cov(%)",
            "Target_Cov(%)", "Percent_Identity", "Alignment_Length", "MapQ"
        ]
        outfh.write("\t".join(header) + "\n")
        for qid, rec in alignments.items():
            line = [
                rec["qname"],
                rec["tname"],
                str(rec["qlen"]),
                f"{rec['qcov']*100:.2f}",
                f"{rec['tcov']*100:.2f}",
                f"{rec['pid']:.2f}",
                str(rec["aln_len"]),
                str(rec["mapq"])
            ]
            outfh.write("\t".join(line) + "\n")

    # 输出汇总统计报告
    summary_file = os.path.join(args.outdir, "summary_report.tsv")
    with open(summary_file, "w") as sf:
        # in tsv format
        sf.write("Metric\tValue\n")
        sf.write(f"Total assembled transcripts\t{total_asm}\n")
        sf.write(f"Total reference transcripts\t{total_ref}\n")
        sf.write(f"Assembled transcripts with alignment\t{num_aligned}\n")
        sf.write(f"Assembled transcripts without alignment\t{missing_asm}\n")
        sf.write(f"Reference transcripts hit\t{len(hit_ref_ids)}\n")
        sf.write(f"Reference transcripts missing\t{missing_ref}\n")
        sf.write(f"Reference Transcript Recovery Rate (RTRR)\t{rtrr:.2f}%\n")
        sf.write(f"Average Query Coverage\t{avg_qcov:.2f}%\n")
        sf.write(f"Average Target Coverage\t{avg_tcov:.2f}%\n")
        sf.write(f"Average Percent Identity\t{avg_pid:.2f}%\n")
        sf.write(f"High quality assembled transcripts\t{high_quality_count}\n")
        sf.write(f"High Quality Ratio (aligned transcripts)\t{hq_ratio_aligned:.2f}%\n")
        sf.write(f"High Quality Ratio (total assembled transcripts)\t{hq_ratio_total:.2f}%\n")
        
    # with open(summary_file, "w") as sf:
    #     sf.write("=== Alignment Summary ===\n")
    #     sf.write(f"Total assembled transcripts: {total_asm}\n")
    #     sf.write(f"Total reference transcripts: {total_ref}\n\n")
    #     sf.write(f"Assembled transcripts with alignment: {num_aligned}\n")
    #     sf.write(f"Assembled transcripts without alignment: {missing_asm}\n\n")
    #     sf.write(f"Reference transcripts hit by at least one alignment: {len(hit_ref_ids)}\n")
    #     sf.write(f"Reference transcripts without alignment: {missing_ref}\n")
    #     sf.write(f"Reference Transcript Recovery Rate (RTRR): {rtrr:.2f}%\n\n")
    #     sf.write("--- Average Quality Metrics for Alignments (aligned transcripts only) ---\n")
    #     sf.write(f"Average Query Coverage: {avg_qcov:.2f}%\n")
    #     sf.write(f"Average Target Coverage: {avg_tcov:.2f}%\n")
    #     sf.write(f"Average Percent Identity: {avg_pid:.2f}%\n\n")
    #     sf.write("--- High-Quality Assembly Ratio ---\n")
    #     sf.write(f"High quality (qcov>={args.min_qcov*100:.0f}%, tcov>={args.min_tcov*100:.0f}%, pid>={args.min_pid:.0f}%) assembled transcripts:\n")
    #     sf.write(f"Number: {high_quality_count}\n")
    #     sf.write(f"High Quality Ratio among aligned transcripts: {hq_ratio_aligned:.2f}%\n")
    #     sf.write(f"High Quality Ratio among total assembled transcripts: {hq_ratio_total:.2f}%\n")
    
    # 打印汇总结果
    print("\n=== Detailed Results ===")
    print("Alignment report written to:", report_file)
    print("Summary report written to:", summary_file)
    print("\nSummary:")
    print(f"Total assembled transcripts: {total_asm}")
    print(f"Total reference transcripts: {total_ref}")
    print(f"Assembled transcripts with alignment: {num_aligned}")
    print(f"Assembled transcripts without alignment: {missing_asm}")
    print(f"Reference transcripts hit: {len(hit_ref_ids)}")
    print(f"Reference transcripts missing: {missing_ref}")
    print(f"Reference Transcript Recovery Rate (RTRR): {rtrr:.2f}%")
    print(f"Average Query Coverage: {avg_qcov:.2f}%")
    print(f"Average Target Coverage: {avg_tcov:.2f}%")
    print(f"Average Percent Identity: {avg_pid:.2f}%")
    print(f"High quality assembled transcripts: {high_quality_count}")
    print(f"High Quality Ratio (aligned transcripts): {hq_ratio_aligned:.2f}%")
    print(f"High Quality Ratio (total assembled transcripts): {hq_ratio_total:.2f}%")
    print("Evaluation complete.")

if __name__ == "__main__":
    main()