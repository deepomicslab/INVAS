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
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, 
                         names=["chromosome", "start", "end", "gene", "score", "strand"])
    return bed_df

def has_overlap(region1, region2):
    """判断两个基因区间是否有重叠"""
    chr1, start1, end1 = region1
    chr2, start2, end2 = region2

    return chr1 == chr2 and max(start1, start2) <= min(end1, end2)

def evaluate_performance(truth_regions, predicted_regions):
    """计算 TP, FN, FP，并进行统计分析"""
    # 从真实基因角度判断是否找到预测
    matched_truth = {
        truth_gene 
        for truth_gene, truth_region in truth_regions.items()
        if any(has_overlap(truth_region, pred_region) for pred_region in predicted_regions.values())
    }
    # 从预测基因角度判断是否其与真实基因重叠
    matched_pred = {
        pred_gene 
        for pred_gene, pred_region in predicted_regions.items()
        if any(has_overlap(pred_region, truth_region) for truth_region in truth_regions.values())
    }

    TP = len(matched_truth)
    FN = len(truth_regions) - len(matched_truth)
    FP = len(predicted_regions) - len(matched_pred)

    # 计算评估指标
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0  
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    print("\n=== 统计分析结果 ===")
    print(f"True Positives (TP): {TP}")
    print(f"False Positives (FP): {FP}")
    print(f"False Negatives (FN): {FN}")

    print("\n=== 评估指标 ===")
    print(f"Precision: {precision:.4f}")
    print(f"Recall (Sensitivity): {recall:.4f}")
    print(f"F1-score: {f1_score:.4f}")
    with open(args.output, 'w') as f:
        # in tsv format
        f.write(f"TP\tFP\tFN\tPrecision\tRecall\tF1-score\n")
        f.write(f"{TP}\t{FP}\t{FN}\t{precision:.4f}\t{recall:.4f}\t{f1_score:.4f}\n")

    # 输出遗漏的基因（FN）
    missing_genes = set(truth_regions.keys()) - matched_truth
    if missing_genes:
        print("\n=== 预测遗漏的基因（FN） ===")
        for gene in sorted(missing_genes):
            print(gene)
    
    # 输出预测错误的基因（FP）
    false_positive_genes = set(predicted_regions.keys()) - matched_pred
    if false_positive_genes:
        print("\n=== 预测错误的基因（FP） ===")
        for gene in sorted(false_positive_genes):
            print(gene)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="使用本地 BED 文件查询基因坐标，并评估预测结果")
    parser.add_argument("--truth", required=True, help="真实 inversion 基因的 TSV 文件")
    parser.add_argument("--predicted", required=True, help="软件预测的基因文件（两列格式）")
    parser.add_argument("--bed", required=True, help="包含所有基因坐标的 BED 文件")
    parser.add_argument("--output", help="输出文件名")
    args = parser.parse_args()

    bed_data = load_bed_file(args.bed)
    truth_regions = load_truth_data(args.truth, bed_data)
    predicted_regions = load_predicted_data(args.predicted, bed_data)

    evaluate_performance(truth_regions, predicted_regions)