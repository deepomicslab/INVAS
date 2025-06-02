import os
import subprocess
import sys

def parse_bed_file_to_regions(bed_file):
    """
    解析 BED 文件，将所有基因的染色体区域合并为一个列表。
    
    Args:
        bed_file (str): BED 文件路径，格式为 "基因名\t染色体区域(chr:start-end)"。
    
    Returns:
        list: 所有染色体区域组成的列表，例如 ["chr1:10000-20000", "chr2:30000-40000"]。
    """
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                print(f"BED 文件格式错误：{line}")
                continue
            region = parts[1]  # 第二列为染色体区域
            regions.append(region)
    return regions

def extract_bam_regions(input_dir, bed_file, sample_names_file, output_dir):
    """
    根据样本名文件提取指定文件夹中样本的 BAM 文件的指定区域（所有基因合并到一个 BAM 文件），输出到目标文件夹。
    
    Args:
        input_dir (str): 样本文件夹主目录。
        bed_file (str): 包含基因名和染色体区域的 BED 文件路径。
        sample_names_file (str): 包含样本名的文本文件路径，每行一个样本名。
        output_dir (str): 提取后 BAM 文件的输出目录。
    """
    # 检查输入路径和 BED 文件是否存在
    if not os.path.exists(input_dir):
        print(f"输入目录 {input_dir} 不存在！")
        sys.exit(1)
    if not os.path.exists(bed_file):
        print(f"BED 文件 {bed_file} 不存在！")
        sys.exit(1)
    if not os.path.exists(sample_names_file):
        print(f"样本名文件 {sample_names_file} 不存在！")
        sys.exit(1)
    
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    # 解析 BED 文件，获取所有区域
    regions = parse_bed_file_to_regions(bed_file)
    if not regions:
        print("BED 文件为空或格式错误，请检查！")
        sys.exit(1)

    # 将区域合并为一个字符串，供 samtools 使用
    regions_string = " ".join(regions)
    print(f"从 BED 文件中读取到 {len(regions)} 个区域：")
    print(regions_string)
    
    # 读取样本名列表
    with open(sample_names_file, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]
    
    if not sample_names:
        print("样本名文件为空，请检查！")
        sys.exit(1)
    
    print(f"样本名文件中读取到 {len(sample_names)} 个样本名：")
    print("\n".join(sample_names))
    
    # 遍历样本名，处理每个样本
    for sample_name in sample_names:
        sample_path = os.path.join(input_dir, sample_name)
        
        # 判断样本文件夹是否存在
        if not os.path.isdir(sample_path):
            print(f"样本文件夹 {sample_name} 不存在，跳过。")
            continue
        
        # # 寻找 sample_name.bam 文件
        # bam_file = os.path.join(sample_path,f"{sample_name}.bam")
        # if not os.path.exists(bam_file):
        #     print(f"样本文件夹 {sample_name} 中未找到 {sample_name}.bam，跳过。")
        #     continue
        
        # search .bam file in the sample folderas the input bam file
        bam_file = ""
        for root, dirs, files in os.walk(sample_path):
            for file in files:
                if file.endswith(".bam"):
                    bam_file = os.path.join(root, file)
                    break
            if bam_file:
                break
        # 创建该样本的输出文件夹
        sample_output_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # 输出的 BAM 文件路径
        output_bam = os.path.join(sample_output_dir, f"{sample_name}_focus.bam")
        
        # samtools view 命令
        cmd = f"samtools view -b {bam_file} {regions_string} -o {output_bam}"
        print(f"正在处理样本 {sample_name}，运行命令：{cmd}")
        
        # 调用子进程运行命令
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"样本 {sample_name} 提取完成，结果保存到 {output_bam}")
            # sort the output bam file
            cmd = f"samtools sort -@ {threads} -o {output_bam.replace('.bam', '_sorted.bam')} {output_bam}"
            subprocess.run(cmd, shell=True, check=True)
            # index the output bam file
            cmd = f"samtools index {output_bam.replace('.bam', '_sorted.bam')}"
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"样本 {sample_name} 提取失败，错误：{e}")

if __name__ == "__main__":
    # 示例：修改为实际路径
    input_directory = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/"           # 样本文件夹主目录
    bed_file_path = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/focus.bed"           # BED 文件路径
    sample_names_path = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/all_wgs.samples"  # 样本名文件路径
    output_directory = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/sub_wgs"         # 输出目录
    threads=16
    
    # 调用主函数
    extract_bam_regions(input_directory, bed_file_path, sample_names_path, output_directory)
