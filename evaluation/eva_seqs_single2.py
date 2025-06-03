#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys

def run_minimap2(query_fasta, target_fasta, out_file, preset="asm5", threads=4):
    """ 运行 minimap2 进行 RNA 转录本比对 """
    cmd = ["minimap2", "-cx", preset, "-t", str(threads), target_fasta, query_fasta]
    print("Running command:", " ".join(cmd))
    try:
        with open(out_file, "w") as fout:
            subprocess.run(cmd, stdout=fout, check=True)
    except subprocess.CalledProcessError as e:
        print("Error running minimap2:", e)
        return False
    return True

def parse_paf(paf_file):
    """ 解析 PAF 文件，仅保留 Primary Alignment """
    best_alignments = {}
    with open(paf_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 12:
                continue
            
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            matches, aln_len, mapq = int(fields[9]), int(fields[10]), int(fields[11])

            qcov = (qend - qstart) / qlen if qlen > 0 else 0.0
            tcov = (tend - tstart) / tlen if tlen > 0 else 0.0
            pid = (matches / aln_len * 100) if aln_len > 0 else 0.0

            optional_fields = {field.split(":")[0]: field for field in fields[12:]}
            if "tp" in optional_fields and optional_fields["tp"] != "tp:A:P":
                continue  # 只保留 Primary Alignment

            record = {
                "qname": qname, "qlen": qlen, "qcov": qcov,
                "tname": tname, "tlen": tlen, "tcov": tcov,
                "matches": matches, "aln_len": aln_len, "pid": pid, "mapq": mapq
            }

            if qname in best_alignments:
                if aln_len > best_alignments[qname]["aln_len"]:
                    best_alignments[qname] = record
            else:
                best_alignments[qname] = record

    return best_alignments

def read_fasta_ids(fasta_file):
    """ 解析 FASTA 文件，返回所有转录本 ID """
    seq_ids = set()
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_ids.add(line[1:].strip().split()[0])
    return seq_ids

def main():
    parser = argparse.ArgumentParser(description="RNA 组装评估：RTRR 和 HQAR 计算")
    parser.add_argument("--ref", required=True, help="真实转录本 FASTA 文件")
    parser.add_argument("--asm", required=True, help="组装转录本 FASTA 文件")
    parser.add_argument("--outdir", default="results", help="结果输出目录")
    parser.add_argument("--threads", type=int, default=4, help="minimap2 线程数")
    parser.add_argument("--preset", default="asm5", help="minimap2 预设参数")
    parser.add_argument("--min_qcov", type=float, default=0.8, help="最小查询覆盖率")
    parser.add_argument("--min_tcov", type=float, default=0.8, help="最小目标覆盖率")
    parser.add_argument("--min_pid", type=float, default=90.0, help="最小百分比相似性")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    asm_ids = read_fasta_ids(args.asm)
    ref_ids = read_fasta_ids(args.ref)
    total_asm, total_ref = len(asm_ids), len(ref_ids)

    paf_file = os.path.join(args.outdir, "asm2ref.paf")
    if not run_minimap2(args.asm, args.ref, paf_file, args.preset, args.threads):
        sys.exit(1)

    alignments = parse_paf(paf_file)
    num_aligned = len(alignments)

    hit_ref_ids = {rec["tname"] for rec in alignments.values()
                   if rec["qcov"] >= args.min_qcov and rec["tcov"] >= args.min_tcov and rec["pid"] >= args.min_pid}
    rtrr = (len(hit_ref_ids) / total_ref * 100) if total_ref > 0 else 0.0

    missing_ref_ids = total_ref - len(hit_ref_ids)

    avg_qcov = sum(rec["qcov"] for rec in alignments.values()) / num_aligned * 100 if num_aligned > 0 else 0.0
    avg_tcov = sum(rec["tcov"] for rec in alignments.values()) / num_aligned * 100 if num_aligned > 0 else 0.0
    avg_pid = sum(rec["pid"] for rec in alignments.values()) / num_aligned if num_aligned > 0 else 0.0

    high_quality_count = len(hit_ref_ids)
    hq_ratio_total = (high_quality_count / total_asm * 100) if total_asm > 0 else 0.0
    hq_ratio_aligned = (high_quality_count / num_aligned * 100) if num_aligned > 0 else 0.0

    summary_file = os.path.join(args.outdir, "summary_report.tsv")
    with open(summary_file, "w") as sf:
        sf.write("Metric\tValue\n")
        sf.write(f"Total assembled transcripts\t{total_asm}\n")
        sf.write(f"Total reference transcripts\t{total_ref}\n")
        sf.write(f"Assembled transcripts with alignment\t{num_aligned}\n")
        sf.write(f"Assembled transcripts without alignment\t{total_asm - num_aligned}\n")
        sf.write(f"Reference transcripts hit\t{len(hit_ref_ids)}\n")
        sf.write(f"Reference transcripts missing\t{missing_ref_ids}\n")
        sf.write(f"Reference Transcript Recovery Rate (RTRR)\t{rtrr:.2f}%\n")
        sf.write(f"Average Query Coverage\t{avg_qcov:.2f}%\n")
        sf.write(f"Average Target Coverage\t{avg_tcov:.2f}%\n")
        sf.write(f"Average Percent Identity\t{avg_pid:.2f}%\n")
        sf.write(f"High quality assembled transcripts\t{high_quality_count}\n")
        sf.write(f"High Quality Ratio (aligned transcripts)\t{hq_ratio_aligned:.2f}%\n")
        sf.write(f"High Quality Ratio (total assembled transcripts)\t{hq_ratio_total:.2f}%\n")

    print(f"Summary report saved in {summary_file}")

if __name__ == "__main__":
    main()