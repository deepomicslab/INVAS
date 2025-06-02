#!/bin/bash
# 单样本SV分析脚本 - scripts_dir默认为当前脚本目录

# 获取脚本所在目录的绝对路径
DEFAULT_SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 定义SV caller处理函数
process_sv_caller() {
    local caller="$1"
    local vcf_file="$2"
    local filtered_vcf="$3"
    local res_dir="$4"
    local wgs_sample="$5"
    
    # 检查VCF文件是否存在
    if [ ! -f "$vcf_file" ]; then
        echo "Warning: $caller VCF not found: $vcf_file"
        return 1
    fi
    
    # 对于需要检查文件大小的caller
    if [[ "$caller" == "lumpy" || "$caller" == "svaba" ]]; then
        if [ ! -s "$vcf_file" ]; then
            echo "Warning: $caller VCF is empty: $vcf_file"
            return 1
        fi
    fi
    
    # 创建结果目录
    if [ ! -d "$res_dir" ]; then
        mkdir -p "$res_dir"
    fi
    
    echo "Processing $caller for $wgs_sample"
    
    # 根据caller执行特定的预处理
    case "$caller" in
        "delly")
            python "$scripts_dir/filter_delly.py" --output_vcf "$filtered_vcf" --input_vcf "$vcf_file"
            if [ $? -eq 0 ]; then
                bash "$scripts_dir/preprocess.sh" "$filtered_vcf" "$ref_bed" "$bwa_bam" "$hisat2_bam" "$res_dir"
            else
                echo "Error: DELLY filtering failed"
                return 1
            fi
            ;;
        "lumpy")
            bcftools view -i 'GT!="0/0" && GT!="." && GT!="./."' "$vcf_file" -o "$filtered_vcf"
            if [ $? -eq 0 ]; then
                bash "$scripts_dir/preprocess.sh" "$filtered_vcf" "$ref_bed" "$bwa_bam" "$hisat2_bam" "$res_dir"
            else
                echo "Error: LUMPY filtering failed"
                return 1
            fi
            ;;
        "manta"|"svaba")
            bash "$scripts_dir/preprocess.sh" "$vcf_file" "$ref_bed" "$bwa_bam" "$hisat2_bam" "$res_dir"
            ;;
        *)
            echo "Error: Unknown caller: $caller"
            return 1
            ;;
    esac
    
    if [ $? -eq 0 ]; then
        echo "✓ $caller processing completed"
        return 0
    else
        echo "✗ $caller processing failed"
        return 1
    fi
}

# 显示使用帮助
show_usage() {
    cat << EOF
Usage: $0 <sample_name> <sv_dir> <hisat2_bam> <bwa_bam> <output_dir> <ref_bed> <callers> [scripts_dir]

Parameters:
  sample_name   : Sample name (e.g., SRR5663703)
  sv_dir        : SV results directory (containing caller subdirs)
  hisat2_bam    : Path to HISAT2 BAM file
  bwa_bam       : Path to BWA unmapped BAM file
  output_dir    : Output directory path for results
  ref_bed       : Reference gene BED file path
  callers       : Comma-separated list of SV callers to use
  scripts_dir   : [Optional] Scripts directory path (default: script location)

Default scripts_dir: $DEFAULT_SCRIPTS_DIR

Expected SV Directory Structure:
  sv_dir/
  ├── delly/SAMPLE.delly.vcf
  ├── manta/SAMPLE.manta.invfmt.vcf
  ├── lumpy/SAMPLE.lumpy.svtyper.vcf
  └── svaba/SAMPLE.svaba.sv.info.vcf

Supported Callers:
  delly, manta, lumpy, svaba

Examples:
  # Use default scripts directory
  $0 SRR5663703 \\
     /path/to/sv_results \\
     /path/to/sample.hisat2.bam \\
     /path/to/sample.unmap.bam \\
     /path/to/output \\
     /path/to/genes.bed \\
     delly,manta,lumpy,svaba
  
  # Use custom scripts directory
  $0 SRR5663703 \\
     /path/to/sv_results \\
     /path/to/sample.hisat2.bam \\
     /path/to/sample.unmap.bam \\
     /path/to/output \\
     /path/to/genes.bed \\
     delly,manta \\
     /custom/scripts/path
EOF
}

# 主处理逻辑
sample_name=$1
sv_dir=$2
hisat2_bam=$3
bwa_bam=$4
output_dir=$5
ref_bed=$6
callers=$7
scripts_dir=${8:-$DEFAULT_SCRIPTS_DIR}  # 如果第8个参数为空，使用默认值

# 参数检查
if [ $# -lt 7 ] || [ $# -gt 8 ]; then
    echo "Error: Invalid number of arguments (got $#, expected 7 or 8)"
    show_usage
    exit 1
fi

# 检查必要目录和文件
if [ ! -d "$sv_dir" ]; then
    echo "Error: SV directory not found: $sv_dir"
    exit 1
fi

if [ ! -f "$hisat2_bam" ]; then
    echo "Error: HISAT2 BAM file not found: $hisat2_bam"
    exit 1
fi

if [ ! -f "$bwa_bam" ]; then
    echo "Error: BWA unmapped BAM file not found: $bwa_bam"
    exit 1
fi

if [ ! -f "$ref_bed" ]; then
    echo "Error: Reference BED file not found: $ref_bed"
    exit 1
fi

if [ ! -d "$scripts_dir" ]; then
    echo "Error: Scripts directory not found: $scripts_dir"
    exit 1
fi

# 检查必要的脚本文件
required_scripts=("filter_delly.py" "preprocess.sh")
missing_scripts=()
for script in "${required_scripts[@]}"; do
    if [ ! -f "$scripts_dir/$script" ]; then
        missing_scripts+=("$script")
    fi
done

if [ ${#missing_scripts[@]} -gt 0 ]; then
    echo "Error: Required scripts not found in $scripts_dir:"
    printf '  %s\n' "${missing_scripts[@]}"
    exit 1
fi

# 验证callers参数
valid_callers=("delly" "manta" "lumpy" "svaba")
IFS=',' read -ra caller_array <<< "$callers"

for caller in "${caller_array[@]}"; do
    if [[ ! " ${valid_callers[@]} " =~ " ${caller} " ]]; then
        echo "Error: Invalid caller '$caller'. Supported callers: ${valid_callers[*]}"
        exit 1
    fi
done

echo "===== Single Sample SV Analysis ====="
echo "Sample: $sample_name"
echo "SV directory: $sv_dir"
echo "HISAT2 BAM: $hisat2_bam"
echo "BWA unmap BAM: $bwa_bam"
echo "Output directory: $output_dir"
echo "Reference BED: $ref_bed"
echo "Scripts directory: $scripts_dir"
echo "Callers to process: $callers"
echo "======================================"

# 从样本名推断WGS样本名（假设相同）
wgs_sample=$sample_name

# 创建样本结果目录
sample_res_dir="$output_dir"/res/$sample_name
if [ ! -d "$sample_res_dir" ]; then
    mkdir -p "$sample_res_dir"
fi

# 记录分析信息
cat > "$sample_res_dir/analysis_info.txt" << EOF
===== Analysis Information =====
Sample: $sample_name
Analysis Date: $(date)
SV Directory: $sv_dir
HISAT2 BAM: $hisat2_bam
BWA Unmap BAM: $bwa_bam
Reference BED: $ref_bed
Scripts Directory: $scripts_dir
Callers: $callers
=================================
EOF

# 处理统计
processed_count=0
failed_count=0
processed_callers=()
failed_callers=()

# 处理每个指定的caller
for caller in "${caller_array[@]}"; do
    echo ""
    echo "--- Processing $caller ---"
    
    case "$caller" in
        "delly")
            delly_vcf="$sv_dir/delly/$sample_name.delly.vcf"
            delly_vcf_f="$sv_dir/delly/$sample_name.delly.f.vcf"
            delly_res_dir="$sample_res_dir"/delly/
            
            if process_sv_caller "delly" "$delly_vcf" "$delly_vcf_f" "$delly_res_dir" "$wgs_sample"; then
                processed_callers+=("delly")
                ((processed_count++))
            else
                failed_callers+=("delly")
                ((failed_count++))
            fi
            ;;
            
        "manta")
            manta_vcf="$sv_dir/manta/$sample_name.manta.invfmt.vcf"
            manta_res_dir="$sample_res_dir"/manta/
            
            if process_sv_caller "manta" "$manta_vcf" "" "$manta_res_dir" "$wgs_sample"; then
                processed_callers+=("manta")
                ((processed_count++))
            else
                failed_callers+=("manta")
                ((failed_count++))
            fi
            ;;
            
        "lumpy")
            lumpy_vcf="$sv_dir/lumpy/$sample_name.lumpy.svtyper.vcf"
            filtered_lumpy_vcf="$sv_dir/lumpy/$sample_name.lumpy.svtyper.f.vcf"
            lumpy_res_dir="$sample_res_dir"/lumpy/
            
            if process_sv_caller "lumpy" "$lumpy_vcf" "$filtered_lumpy_vcf" "$lumpy_res_dir" "$wgs_sample"; then
                processed_callers+=("lumpy")
                ((processed_count++))
            else
                failed_callers+=("lumpy")
                ((failed_count++))
            fi
            ;;
            
        "svaba")
            svaba_vcf="$sv_dir/svaba/$sample_name.svaba.sv.info.vcf"
            svaba_res_dir="$sample_res_dir"/svaba/
            
            if process_sv_caller "svaba" "$svaba_vcf" "" "$svaba_res_dir" "$wgs_sample"; then
                processed_callers+=("svaba")
                ((processed_count++))
            else
                failed_callers+=("svaba")
                ((failed_count++))
            fi
            ;;
    esac
done

# 生成处理报告
report_file="$sample_res_dir/processing_report.txt"
cat > "$report_file" << EOF
===== Processing Summary =====
Sample: $sample_name
Analysis Date: $(date)
Total callers requested: ${#caller_array[@]}
Successfully processed: $processed_count
Failed: $failed_count

Successful callers: ${processed_callers[*]}
Failed callers: ${failed_callers[*]}

Results directory: $sample_res_dir
===============================
EOF

echo ""
echo "===== Processing Summary ====="
echo "Sample: $sample_name"
echo "Total callers requested: ${#caller_array[@]}"
echo "Successfully processed: $processed_count"
echo "Failed: $failed_count"

if [ ${#processed_callers[@]} -gt 0 ]; then
    echo "Successful callers: ${processed_callers[*]}"
fi

if [ ${#failed_callers[@]} -gt 0 ]; then
    echo "Failed callers: ${failed_callers[*]}"
fi

echo "Results directory: $sample_res_dir"
echo "Report saved: $report_file"
echo "==============================="

# 如果所有caller都失败，退出时返回错误码
if [ $processed_count -eq 0 ]; then
    echo "Error: All callers failed to process"
    exit 1
fi

echo "Processing completed for sample $sample_name"