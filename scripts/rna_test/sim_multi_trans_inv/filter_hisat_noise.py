import argparse
import pysam

def filter_bam(input_bam_path, output_bam_path):
    # 打开原始BAM文件进行读取
    with pysam.AlignmentFile(input_bam_path, "rb") as original_bam:
        # 创建一个新的BAM文件用于输出过滤后的记录
        with pysam.AlignmentFile(output_bam_path, "wb", template=original_bam) as filtered_bam:
            # 遍历原始BAM文件中的每条记录
            for read in original_bam.fetch():
                # 获取CIGAR字符串
                cigar_string = read.cigarstring
                
                # 检查CIGAR字符串是否不包含'N'并且包含'S'
                if cigar_string and 'N' not in cigar_string and 'S' in cigar_string:
                    continue  # 如果满足条件，跳过这条记录
                else:
                    # 如果不满足条件，将记录写入新的BAM文件
                    filtered_bam.write(read)

    print(f"过滤完成。输出文件: {output_bam_path}")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Filter BAM records based on CIGAR string.')
    parser.add_argument('input_bam', type=str, help='Input BAM file path')
    parser.add_argument('output_bam', type=str, help='Output BAM file path')

    # 解析命令行参数
    args = parser.parse_args()

    # 过滤BAM文件
    filter_bam(args.input_bam, args.output_bam)

if __name__ == '__main__':
    main()