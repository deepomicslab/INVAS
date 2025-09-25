import os
import subprocess
import sys

def parse_bed_file_to_regions(bed_file):
    """
    parse a BED file and extract the chromosome regions from the second column.
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
            region = parts[1]
            regions.append(region)
    return regions

def extract_bam_regions(input_dir, bed_file, sample_names_file, output_dir):
    """
    extract specified regions from BAM files for each sample listed in sample_names_file.
    """
    # check if input paths exist
    if not os.path.exists(input_dir):
        print(f"Input directory {input_dir} does not exist!")
        sys.exit(1)
    if not os.path.exists(bed_file):
        print(f"BED file {bed_file} does not exist!")
        sys.exit(1)
    if not os.path.exists(sample_names_file):
        print(f"Sample names file {sample_names_file} does not exist!")
        sys.exit(1)

    # create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # parse BED file to get all regions
    regions = parse_bed_file_to_regions(bed_file)
    if not regions:
        print("BED file is empty or incorrectly formatted, please check!")
        sys.exit(1)

    # make a single string of regions for samtools
    regions_string = " ".join(regions)
    print(f"regions to extract: {len(regions)}")
    print(regions_string)
    
    # read sample names from the file
    with open(sample_names_file, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]
    
    if not sample_names:
        print("Sample names file is empty, please check!")
        sys.exit(1)
    
    print(f"get {len(sample_names)} samples from {sample_names_file}:")
    print("\n".join(sample_names))
    
    # iterate over each sample and extract regions
    for sample_name in sample_names:
        sample_path = os.path.join(input_dir, sample_name)
        
        # check if sample directory exists
        if not os.path.isdir(sample_path):
            print(f"Sample directory {sample_name} does not exist, skipping.")
            continue        
        # search .bam file in the sample folderas the input bam file
        bam_file = ""
        for root, dirs, files in os.walk(sample_path):
            for file in files:
                if file.endswith(".bam"):
                    bam_file = os.path.join(root, file)
                    break
            if bam_file:
                break
        # create output directory for the sample
        sample_output_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)

        # output BAM file path
        output_bam = os.path.join(sample_output_dir, f"{sample_name}_focus.bam")
        
        # samtools view command
        cmd = f"samtools view -b {bam_file} {regions_string} -o {output_bam}"
        print(f"Processing sample: {sample_name}")

        # call subprocess to run command
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"Sample {sample_name} extraction completed, results saved to {output_bam}")
            # sort the output bam file
            cmd = f"samtools sort -@ {threads} -o {output_bam.replace('.bam', '_sorted.bam')} {output_bam}"
            subprocess.run(cmd, shell=True, check=True)
            # index the output bam file
            cmd = f"samtools index {output_bam.replace('.bam', '_sorted.bam')}"
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Sample {sample_name} extraction failed, error: {e}")

if __name__ == "__main__":

    input_directory = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/"        
    bed_file_path = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/focus.bed"          
    sample_names_path = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/all_wgs.samples"  
    output_directory = "/scratch/project/cs_shuaicli/wxd/TCGA_LIVER_LIHC/part1/sub_wgs"        
    threads=16
    
    # 调用主函数
    extract_bam_regions(input_directory, bed_file_path, sample_names_path, output_directory)
