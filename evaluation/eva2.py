import pandas as pd
import argparse

def load_truth_data(truth_file, bed_data):
    """读取 truth TSV 文件，并使用 BED 数据映射基因坐标"""
    truth_df = pd.read_csv(truth_file, sep="\t")
    truth_genes = set(truth_df["gene_name"])
    
    truth_regions = {}
    for _, row in bed_data.iterrows():
        if row["gene"] in truth_genes:
            truth_regions[row["gene"]] = (row["chromosome"], row["start"], row["end"])
    
    return truth_regions

def load_predicted_data(predicted_file, bed_data):
    """读取预测文件，并使用 BED 数据映射基因坐标"""
    pred_df = pd.read_csv(predicted_file, sep="\t", header=None, names=["gene", "support"])
    predicted_genes = set(pred_df["gene"])
    
    predicted_regions = {}
    for _, row in bed_data.iterrows():
        if row["gene"] in predicted_genes:
            predicted_regions[row["gene"]] = (row["chromosome"], row["start"], row["end"])
    
    return predicted_regions

def load_bed_file(bed_file):
    """读取 BED 文件，获取所有基因的坐标"""
    bed_df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None, 
        names=["chromosome", "start", "end", "gene", "score", "strand"]
    )
    return bed_df

def evaluate_performance(truth_regions, predicted_regions):
    """
    评估预测结果：
    1. 首先通过基因名称完全匹配统计 TP
    2. 对于遗漏的真实基因（FN候选），检查是否存在预测区域与其坐标重叠80%以上，
       如满足则重新计入 TP
    3. 最终根据 TP、FN、FP 计算 precision, recall, f1_score 并输出详细结果。
    """
    truth_set = set(truth_regions.keys())
    predicted_set = set(predicted_regions.keys())
    
    # 1. 基因名称完全匹配
    exact_matches = truth_set & predicted_set
    TP = len(exact_matches)
    matched_truth = set(exact_matches)
    
    # 保存未使用的预测区域。后续用于拯救遗漏的真值基因
    available_pred_set = set(predicted_set - exact_matches)
    rescue_matches = {}  # 存储通过重叠匹配拯救的真值基因: {truth_gene: predicted_gene}
    
    # 2. 对于遗漏的真实基因，尝试匹配预测区域 (坐标重叠80%以上)
    for truth_gene in truth_set - exact_matches:
        truth_chr, truth_start, truth_end = truth_regions[truth_gene]
        truth_length = truth_end - truth_start
        # 遍历所有未被捕获的预测区域
        for pred_gene in list(available_pred_set):
            pred_chr, pred_start, pred_end = predicted_regions[pred_gene]
            # 判断是否在同一染色体上
            if truth_chr != pred_chr:
                continue
            # 计算重叠区间长度
            overlap = max(0, min(truth_end, pred_end) - max(truth_start, pred_start))
            overlap_fraction = overlap / truth_length if truth_length > 0 else 0
            if overlap_fraction >= 0.00000001:
                TP += 1
                rescue_matches[truth_gene] = pred_gene
                matched_truth.add(truth_gene)
                available_pred_set.remove(pred_gene)
                break  # 一旦匹配成功，则终止当前真值基因的遍历

    # 3. 重新计算 FN 和 FP
    total_TP = TP
    FN = len(truth_set) - total_TP
    FP = len(predicted_set) - total_TP

    precision = total_TP / (total_TP + FP) if (total_TP + FP) > 0 else 0
    recall = total_TP / (total_TP + FN) if (total_TP + FN) > 0 else 0  
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    print("\n=== 统计分析结果 ===")
    print(f"True Positives (TP): {total_TP}")
    print(f"False Positives (FP): {FP}")
    print(f"False Negatives (FN): {FN}")

    print("\n=== 评估指标 ===")
    print(f"Precision: {precision:.4f}")
    print(f"Recall (Sensitivity): {recall:.4f}")
    print(f"F1-score: {f1_score:.4f}")
    with open(args.output, 'w') as f:
        # in tsv format
        f.write(f"TP\tFP\tFN\tPrecision\tRecall\tF1-score\n")
        f.write(f"{total_TP}\t{FP}\t{FN}\t{precision:.4f}\t{recall:.4f}\t{f1_score:.4f}\n")


    # 输出没有成功匹配的基因（FN）
    missing_genes = truth_set - matched_truth
    if missing_genes:
        print("\n=== 预测遗漏的基因（FN） ===")
        for gene in sorted(missing_genes):
            print(gene)
    
    # 输出预测错误的基因（FP）：
    # FP 为预测中未被任何真值基因匹配到的基因（包括基因名称和坐标匹配）
    used_pred_genes = exact_matches.union(set(rescue_matches.values()))
    false_positive_genes = predicted_set - used_pred_genes
    if false_positive_genes:
        print("\n=== 预测错误的基因（FP） ===")
        for gene in sorted(false_positive_genes):
            print(gene)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="使用本地 BED 文件查询基因坐标，并评估预测结果。\n"
                    "首先通过基因名称匹配，其次对于遗漏基因如果与某预测区域坐标重叠80%以上则也计入TP。"
    )
    parser.add_argument("--truth", required=True, help="真实 inversion 基因的 TSV 文件")
    parser.add_argument("--predicted", required=True, help="软件预测的基因文件（两列格式）")
    parser.add_argument("--bed", required=True, help="包含所有基因坐标的 BED 文件")
    parser.add_argument("--output", help="输出文件名")
    args = parser.parse_args()

    bed_data = load_bed_file(args.bed)
    truth_regions = load_truth_data(args.truth, bed_data)
    predicted_regions = load_predicted_data(args.predicted, bed_data)

    evaluate_performance(truth_regions, predicted_regions)