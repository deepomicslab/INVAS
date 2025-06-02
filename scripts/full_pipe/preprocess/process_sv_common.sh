#!/bin/bash

# Simple SV Detection Pipeline
# Optimized version with proper manta workflow

# 安全检查：防止shell层级过深
# if [[ ${SHLVL:-0} -gt 5 ]]; then
#     echo "ERROR: Shell level too high ($SHLVL). Starting fresh shell."
#     exec bash -l "$0" "$@"
# fi

# 防止重复执行
if [[ -n "$SV_PIPELINE_RUNNING" ]]; then
    echo "ERROR: SV pipeline already running in this session"
    exit 1
fi
export SV_PIPELINE_RUNNING=1

set -e

# 清理函数
cleanup() {
    unset SV_PIPELINE_RUNNING
    exit ${1:-0}
}
trap 'cleanup $?' EXIT

# ===== 工具路径配置 =====
# 根据您的环境修改这些路径
SAMTOOLS="${SAMTOOLS:-samtools}"
DELLY="${DELLY:-/home/xuedowang2/app/delly_v1.1.6_linux_x86_64bit}"
# MANTA_CONFIG="/home/xuedowang2/miniconda3/envs/py27/bin/configManta.py"
# MANTA_CONVERT="/home/xuedowang2/miniconda3/envs/py27/bin/convertInversion.py"
MANTA_CONFIG="${MANTA_CONFIG:-configManta.py}"
MANTA_CONVERT="${MANTA_CONVERT:-convertInversion.py}"
SVABA="${SVABA:-svaba}"
LUMPYEXPRESS="${LUMPYEXPRESS:-lumpyexpress}"
LUMPY_EXTRACT="${LUMPY_EXTRACT:-extractSplitReads_BwaMem}"
SVTYPER="${SVTYPER:-svtyper}"
# svaba_converter=/home/xuedowang2/app/svaba_converter.py
SVABA_CONVERTER="${SVABA_CONVERTER:-/home/xuedowang2/app/svaba_converter.py}"

# ===== 默认参数 =====
THREADS=8
CALLERS="delly,manta,svaba,lumpy"
VERBOSE=false

# ===== 帮助信息 =====
show_help() {
    cat << EOF
Usage: $0 -i input.bam -o output_dir -r reference.fa -s sample_name [options]

Required Parameters:
  -i  Input BAM file path
  -o  Output directory  
  -r  Reference genome FASTA file
  -s  Sample name

Optional Parameters:
  -t  Number of threads (default: 8)
  -c  SV callers: delly,manta,svaba,lumpy (default: all)
  -v  Verbose output
  -h  Show this help

Examples:
  $0 -i sample.bam -o ./results -r hg38.fa -s SAMPLE01
  $0 -i sample.bam -o ./results -r hg38.fa -s SAMPLE01 -t 16 -c delly,manta
  $0 -i sample.bam -o ./results -r hg38.fa -s SAMPLE01 -v

Output Structure:
  results/
  ├── delly/
  │   ├── SAMPLE01.delly.inv.vcf
  │   └── SAMPLE01.delly.inv.log
  ├── manta/
  │   ├── SAMPLE01.manta.invfmt.vcf
  │   └── results/variants/diploidSV.vcf.gz
  ├── svaba/
  │   └── SAMPLE01.svaba.sv.vcf
  └── lumpy/
      ├── SAMPLE01.lumpy.vcf
      └── SAMPLE01.lumpy.genotyped.vcf
EOF
}

# ===== 日志函数 =====
log_info() {
    echo "[$(date '+%H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2
}

log_debug() {
    if [[ "$VERBOSE" == true ]]; then
        echo "[$(date '+%H:%M:%S')] DEBUG: $*" >&2
    fi
}

# ===== 参数解析 =====
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i)
                INPUT_BAM="$2"
                shift 2
                ;;
            -o)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -r)
                REFERENCE="$2"
                shift 2
                ;;
            -s)
                SAMPLE_NAME="$2"
                shift 2
                ;;
            -t)
                THREADS="$2"
                shift 2
                ;;
            -c)
                CALLERS="$2"
                shift 2
                ;;
            -v)
                VERBOSE=true
                shift
                ;;
            -h)
                show_help
                cleanup 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_help
                cleanup 1
                ;;
        esac
    done
}

# ===== 参数验证 =====
validate_inputs() {
    local errors=()
    
    # 检查必需参数
    [[ -z "${INPUT_BAM:-}" ]] && errors+=("Input BAM file not specified")
    [[ -z "${OUTPUT_DIR:-}" ]] && errors+=("Output directory not specified")
    [[ -z "${REFERENCE:-}" ]] && errors+=("Reference genome not specified")
    [[ -z "${SAMPLE_NAME:-}" ]] && errors+=("Sample name not specified")
    
    # 检查文件存在性
    [[ -n "${INPUT_BAM:-}" && ! -f "$INPUT_BAM" ]] && errors+=("BAM file not found: $INPUT_BAM")
    [[ -n "${REFERENCE:-}" && ! -f "$REFERENCE" ]] && errors+=("Reference file not found: $REFERENCE")
    
    # 检查数值参数
    [[ ! "$THREADS" =~ ^[0-9]+$ ]] && errors+=("Invalid thread count: $THREADS")
    
    if [[ ${#errors[@]} -gt 0 ]]; then
        log_error "Validation failed:"
        printf '  %s\n' "${errors[@]}"
        cleanup 1
    fi
}

# ===== 工具检查 =====
check_tools() {
    local missing_tools=()
    
    # 检查基础工具
    command -v "$SAMTOOLS" >/dev/null || missing_tools+=("samtools: $SAMTOOLS")
    
    # 检查选择的SV检测工具
    IFS=',' read -ra CALLER_ARRAY <<< "$CALLERS"
    for caller in "${CALLER_ARRAY[@]}"; do
        case $caller in
            delly)
                command -v "$DELLY" >/dev/null || missing_tools+=("delly: $DELLY")
                ;;
            manta)
                [[ -f "$MANTA_CONFIG" ]] || missing_tools+=("manta configManta.py: $MANTA_CONFIG")
                [[ -f "$MANTA_CONVERT" ]] || missing_tools+=("manta convertInversion.py: $MANTA_CONVERT")
                ;;
            svaba)
                command -v "$SVABA" >/dev/null || missing_tools+=("svaba: $SVABA")
                ;;
            lumpy)
                command -v "$LUMPYEXPRESS" >/dev/null || missing_tools+=("lumpyexpress: $LUMPYEXPRESS")
                command -v "$LUMPY_EXTRACT" >/dev/null || missing_tools+=("extractSplitReads_BwaMem: $LUMPY_EXTRACT")
                ;;
        esac
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing tools:"
        printf '  %s\n' "${missing_tools[@]}"
        log_error "Please install missing tools or correct the paths in the script"
        cleanup 1
    fi
    
    log_info "All required tools found"
}

# ===== 主函数 =====
main() {
    local start_time=$(date +%s)
    
    log_info "===== SV Detection Pipeline ====="
    log_info "Sample: $SAMPLE_NAME"
    log_info "Input BAM: $INPUT_BAM"
    log_info "Output Directory: $OUTPUT_DIR"
    log_info "Reference: $REFERENCE"
    log_info "Threads: $THREADS"
    log_info "Callers: $CALLERS"
    log_info "=============================="
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR"
    cd "$OUTPUT_DIR" || cleanup 1
    
    # 检查BAM索引
    check_bam_index
    
    # 运行选择的SV检测工具
    IFS=',' read -ra CALLER_ARRAY <<< "$CALLERS"
    for caller in "${CALLER_ARRAY[@]}"; do
        log_info "Starting $caller..."
        case $caller in
            delly)
                run_delly
                ;;
            manta)
                run_manta
                ;;
            svaba)
                run_svaba
                ;;
            lumpy)
                run_lumpy
                ;;
            *)
                log_error "Unknown caller: $caller"
                ;;
        esac
    done
    
    # 生成汇总报告
    generate_summary
    
    # 计算运行时间
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    log_info "===== SV Detection Complete ====="
    log_info "Total runtime: $(printf '%02d:%02d:%02d' $((runtime/3600)) $((runtime%3600/60)) $((runtime%60)))"
}

# ===== BAM索引检查 =====
check_bam_index() {
    log_info "Checking BAM file index..."
    if [[ ! -f "${INPUT_BAM}.bai" ]]; then
        log_info "Creating BAM index..."
        $SAMTOOLS index -@ "$THREADS" "$INPUT_BAM"
        log_info "BAM index created"
    else
        log_info "BAM index found"
    fi
}

# ===== DELLY =====
run_delly() {
    log_info "Running DELLY (inversions only)..."
    local delly_dir="$OUTPUT_DIR/delly"
    mkdir -p "$delly_dir"
    
    local delly_vcf="$delly_dir/${SAMPLE_NAME}.delly.vcf"
    local delly_log="$delly_dir/${SAMPLE_NAME}.delly.log"
    
    log_debug "DELLY command: $DELLY call -t INV -g $REFERENCE $INPUT_BAM"
    
    if $DELLY call -t INV -g "$REFERENCE" "$INPUT_BAM" \
        > "$delly_vcf" 2>"$delly_log"; then
        
        local inv_count
        inv_count=$(grep -cv '^#' "$delly_vcf" 2>/dev/null || echo "0")
        log_info "DELLY completed: $inv_count inversions detected"
        
    else
        log_error "DELLY failed, check log: $delly_log"
        return 1
    fi
}

# ===== MANTA =====
run_manta() {
    log_info "Running MANTA..."
    local manta_dir="$OUTPUT_DIR/manta"
    mkdir -p "$manta_dir"
    
    # 配置MANTA
    log_info "  Configuring MANTA..."
    log_debug "MANTA config command: $MANTA_CONFIG --bam $INPUT_BAM --referenceFasta $REFERENCE --runDir $manta_dir"
    
    if ! "$MANTA_CONFIG" \
        --bam "$INPUT_BAM" \
        --referenceFasta "$REFERENCE" \
        --runDir "$manta_dir"; then
        log_error "MANTA configuration failed"
        return 1
    fi
    
    # 运行MANTA工作流
    log_info "  Running MANTA workflow..."
    log_debug "MANTA run command: python $manta_dir/runWorkflow.py -j $THREADS"
    
    if ! python "$manta_dir/runWorkflow.py" -j "$THREADS"; then
        log_error "MANTA workflow failed"
        return 1
    fi
    
    # 检查输出文件
    local diploid_vcf="$manta_dir/results/variants/diploidSV.vcf.gz"
    if [[ ! -f "$diploid_vcf" ]]; then
        log_error "MANTA output file not found: $diploid_vcf"
        return 1
    fi
    
    # 转换倒位格式
    log_info "  Converting inversion format..."
    local output_vcf="$manta_dir/${SAMPLE_NAME}.manta.invfmt.vcf"
    log_debug "Conversion command: $MANTA_CONVERT $SAMTOOLS $REFERENCE $diploid_vcf"
    
    if "$MANTA_CONVERT" \
        "$SAMTOOLS" \
        "$REFERENCE" \
        "$diploid_vcf" \
        > "$output_vcf"; then
        
        local sv_count
        sv_count=$(grep -cv '^#' "$output_vcf" 2>/dev/null || echo "0")
        log_info "MANTA completed: $sv_count variants detected"
        
    else
        log_error "MANTA inversion conversion failed"
        return 1
    fi
}

# ===== SVABA =====
run_svaba() {
    log_info "Running SVABA..."
    local svaba_dir="$OUTPUT_DIR/svaba"
    mkdir -p "$svaba_dir"
    
    local current_dir=$(pwd)
    cd "$svaba_dir" || return 1
    
    log_debug "SVABA command: $SVABA run -a $SAMPLE_NAME -G $REFERENCE -t $INPUT_BAM --override-reference-check --read-tracking --germline -p $THREADS -L 6 -I"
    
    if $SVABA run \
        -a "$SAMPLE_NAME" \
        -G "$REFERENCE" \
        -t "$INPUT_BAM" \
        --override-reference-check \
        --read-tracking \
        --germline \
        -p "$THREADS" \
        -L 6 \
        -I; then
        
        cd "$current_dir"
        
        local sv_vcf="$svaba_dir/${SAMPLE_NAME}.svaba.sv.vcf"
        if [[ -f "$sv_vcf" ]]; then
            local sv_count
            sv_count=$(grep -cv '^#' "$sv_vcf" 2>/dev/null || echo "0")
            log_info "SVABA completed: $sv_count variants detected"
        else
            log_error "SVABA output file not found: $sv_vcf"
            return 1
        fi
        # CONVERT breakpoint SV TO INFO VCF
        python2 $SVABA_CONVERTER $sv_vcf "$svaba_dir/${SAMPLE_NAME}.inv.bed" > $svaba_dir/$SAMPLE_NAME.svaba.sv.info.vcf

        
    else
        cd "$current_dir"
        log_error "SVABA execution failed"
        return 1
    fi
}

# ===== LUMPY =====
run_lumpy() {
    log_info "Running LUMPY..."
    local lumpy_dir="$OUTPUT_DIR/lumpy"
    mkdir -p "$lumpy_dir"
    
    local current_dir=$(pwd)
    cd "$lumpy_dir" || return 1
    
    # 提取不一致reads
    log_info "  Extracting discordant reads..."
    $SAMTOOLS view -b -F 1294 -@ "$THREADS" "$INPUT_BAM" \
        > "${SAMPLE_NAME}.discordants.unsorted.bam"
    
    # 提取分裂reads
    log_info "  Extracting split reads..."
    $SAMTOOLS view -h "$INPUT_BAM" \
        | $LUMPY_EXTRACT -i stdin \
        | $SAMTOOLS view -Sb -@ "$THREADS" -o "${SAMPLE_NAME}.splitters.unsorted.bam"
    
    # 排序比对文件
    log_info "  Sorting alignment files..."
    $SAMTOOLS sort -@ "$THREADS" "${SAMPLE_NAME}.discordants.unsorted.bam" \
        -o "${SAMPLE_NAME}.discordants.bam"
    $SAMTOOLS index "${SAMPLE_NAME}.discordants.bam"
    
    $SAMTOOLS sort -@ "$THREADS" "${SAMPLE_NAME}.splitters.unsorted.bam" \
        -o "${SAMPLE_NAME}.splitters.bam"
    $SAMTOOLS index "${SAMPLE_NAME}.splitters.bam"
    
    # 运行LUMPY
    log_info "  Running LumpyExpress..."
    log_debug "LUMPY command: $LUMPYEXPRESS -B $INPUT_BAM -S ${SAMPLE_NAME}.splitters.bam -D ${SAMPLE_NAME}.discordants.bam -o ${SAMPLE_NAME}.lumpy.vcf"
    
    if $LUMPYEXPRESS \
        -B "$INPUT_BAM" \
        -S "${SAMPLE_NAME}.splitters.bam" \
        -D "${SAMPLE_NAME}.discordants.bam" \
        -o "${SAMPLE_NAME}.lumpy.vcf"; then
        
        # SVTyper基因型分析（可选）
        if command -v "$SVTYPER" >/dev/null; then
            log_info "  Running SVTyper for genotyping..."
            if $SVTYPER \
                -B "$INPUT_BAM" \
                -S "${SAMPLE_NAME}.splitters.bam" \
                -i "${SAMPLE_NAME}.lumpy.vcf" \
                -o "${SAMPLE_NAME}.lumpy.genotyped.vcf"; then
                log_debug "SVTyper completed successfully"
            else
                log_error "SVTyper failed, but continuing..."
            fi
        fi
        
        local sv_count
        sv_count=$(grep -cv '^#' "${SAMPLE_NAME}.lumpy.vcf" 2>/dev/null || echo "0")
        log_info "LUMPY completed: $sv_count variants detected"
        
    else
        cd "$current_dir"
        log_error "LUMPY execution failed"
        return 1
    fi
    
    cd "$current_dir"
    
    # 清理临时文件
    rm -f "$lumpy_dir"/*.unsorted.bam
}

# ===== 汇总报告 =====
generate_summary() {
    log_info "Generating summary report..."
    
    local summary_file="$OUTPUT_DIR/${SAMPLE_NAME}.sv_summary.txt"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    cat > "$summary_file" << EOF
===== Structural Variant Detection Summary =====
Sample: $SAMPLE_NAME
Analysis Date: $timestamp
Input BAM: $INPUT_BAM
Reference: $REFERENCE
Threads: $THREADS
Callers: $CALLERS

===== Results Overview =====
EOF
    
    local total_variants=0
    
    IFS=',' read -ra CALLER_ARRAY <<< "$CALLERS"
    for caller in "${CALLER_ARRAY[@]}"; do
        local count=0
        local vcf_file=""
        
        case $caller in
            delly)
                vcf_file="delly/${SAMPLE_NAME}.delly.inv.vcf"
                if [[ -f "$vcf_file" ]]; then
                    count=$(grep -cv '^#' "$vcf_file" 2>/dev/null || echo "0")
                    echo "DELLY (inversions): $count variants" | tee -a "$summary_file"
                fi
                ;;
            manta)
                vcf_file="manta/${SAMPLE_NAME}.manta.invfmt.vcf"
                if [[ -f "$vcf_file" ]]; then
                    count=$(grep -cv '^#' "$vcf_file" 2>/dev/null || echo "0")
                    echo "MANTA: $count variants" | tee -a "$summary_file"
                fi
                ;;
            svaba)
                vcf_file="svaba/${SAMPLE_NAME}.svaba.sv.vcf"
                if [[ -f "$vcf_file" ]]; then
                    count=$(grep -cv '^#' "$vcf_file" 2>/dev/null || echo "0")
                    echo "SVABA: $count variants" | tee -a "$summary_file"
                fi
                ;;
            lumpy)
                vcf_file="lumpy/${SAMPLE_NAME}.lumpy.vcf"
                if [[ -f "$vcf_file" ]]; then
                    count=$(grep -cv '^#' "$vcf_file" 2>/dev/null || echo "0")
                    echo "LUMPY: $count variants" | tee -a "$summary_file"
                fi
                ;;
        esac
        
        total_variants=$((total_variants + count))
    done
    
    echo "" | tee -a "$summary_file"
    echo "Total variants detected: $total_variants" | tee -a "$summary_file"
    echo "" | tee -a "$summary_file"
    echo "===== Output Files =====" >> "$summary_file"
    
    # 列出所有VCF文件
    find . -name "*.vcf" -o -name "*.vcf.gz" | sort >> "$summary_file"
    
    log_info "Summary saved to: $summary_file"
}

# ===== 运行流程 =====
parse_arguments "$@"
validate_inputs
check_tools
main