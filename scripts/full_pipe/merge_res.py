import argparse
from collections import defaultdict

def parse_gtf(file_path, skip_header_lines=0):
    """
    解析 GTF 文件，跳过头部行，返回每个 transcript 的完整信息（包括 transcript 和 exon 行）。
    """
    transcripts = defaultdict(lambda: {"transcript_line": None, "exon_lines": []})

    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#") or skip_header_lines > 0:  # 跳过 header
                skip_header_lines -= 1
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # 忽略无效行
            
            feature_type = fields[2]  # "transcript" 或 "exon"
            attributes = fields[8]
            
            # 提取 transcript_id
            transcript_id = extract_attribute(attributes, "transcript_id")
            print(transcript_id)
            if not transcript_id:
                continue  # 如果没有 transcript_id，跳过该行
            
            # 根据 feature 类型存储
            if feature_type == "transcript":
                transcripts[transcript_id]["transcript_line"] = line.strip().replace("StringTie", "Seqflow")
            elif feature_type == "exon":
                transcripts[transcript_id]["exon_lines"].append(line.strip().replace("StringTie", "Seqflow"))
    
    return transcripts

def extract_attribute(attributes, key):
    """
    从 GTF 的 attribute 字段中提取指定的 key（如 transcript_id）。
    """
    for attr in attributes.split(";"):
        attr = attr.strip()
        if attr.startswith(key):
            return attr.split('"')[1]
    return None

def merge_gtf(file1, file2, output_file, skip_header_file1=1, skip_header_file2=2):
    """
    合并两个 GTF 文件，确保每个 transcript 的完整性（包括 transcript 和 exon）。
    """
    # 解析两个 GTF 文件
    transcripts1 = parse_gtf(file1, skip_header_lines=skip_header_file1)
    transcripts2 = parse_gtf(file2, skip_header_lines=skip_header_file2)
    
    # 合并数据
    merged_transcripts = {}

    # 先加入第一个文件的内容
    merged_transcripts.update(transcripts1)

    
    # 再加入第二个文件的内容（覆盖第一个文件中重复的 transcript）
    merged_transcripts.update(transcripts2)
    print(merged_transcripts)

    # format exon line , with start < end and filter exon length < 2
    for transcript_id, transcript_data in merged_transcripts.items():
        for i in range(len(transcript_data["exon_lines"])):
            fields = transcript_data["exon_lines"][i].strip().split("\t")
            start = int(fields[3])
            end = int(fields[4])
            if start > end:
                fields[3], fields[4] = fields[4], fields[3]
            if end - start < 2:
                transcript_data["exon_lines"][i] = None
            else:
                transcript_data["exon_lines"][i] = "\t".join(fields)
        transcript_data["exon_lines"] = list(filter(None, transcript_data["exon_lines"]))
        
    
    # 写入合并结果
    with open(output_file, "w") as out_f:
        for transcript_id, transcript_data in merged_transcripts.items():
            # 写入 transcript 行
            if transcript_data["transcript_line"]:
                out_f.write(transcript_data["transcript_line"] + "\n")
            # 写入对应的 exon 行
            for exon_line in transcript_data["exon_lines"]:
                out_f.write(exon_line + "\n")
    
    print(f"Merged GTF file has been saved to: {output_file}")

# 主函数
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge two GTF files and retain complete transcripts (transcript + exon lines).")
    parser.add_argument("file1", type=str, help="Path to the first GTF file (1-line header).")
    parser.add_argument("file2", type=str, help="Path to the second GTF file (2-line header).")
    parser.add_argument("output_file", type=str, help="Path to save the merged GTF file.")
    args = parser.parse_args()
    
    merge_gtf(args.file1, args.file2, args.output_file)