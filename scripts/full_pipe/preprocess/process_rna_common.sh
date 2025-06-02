#!/bin/bash

# RNA-seq Re-alignment Pipeline (Single Sample)
# Description: Extract fastq from BAM, clean reads, realign with HISAT2, and BWA rescue
# Version: 1.0

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# =====================================================
# Default Parameters
# =====================================================
SAMPLE_NAME=""
INPUT_BAM=""
OUTPUT_DIR=""
REF_FASTA=""
HISAT_INDEX=""

# Tool and environment settings
THREADS=32
CONDA_ENV=""
PICARD_JAR="picard.jar"
TRIMMOMATIC_ADAPTERS="TruSeq3-PE.fa"

# Processing options
FORCE_OVERWRITE=false
SKIP_CLEANUP=false
VERBOSE=false

# Trimmomatic parameters
TRIM_LEADING=3
TRIM_TRAILING=3
TRIM_MINLEN=36
TRIM_ILLUMINACLIP="2:30:10:2:True"

# =====================================================
# Usage Function
# =====================================================
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Parameters:
  -n, --sample-name NAME      Sample name
  -b, --input-bam FILE        Input BAM file
  -o, --output-dir DIR        Output directory
  -r, --reference FILE        Reference FASTA file
  -x, --hisat-index PREFIX    HISAT2 index prefix

Optional Parameters:
  -t, --threads INT           Number of threads (default: $THREADS)
  -e, --conda-env NAME        Conda environment name (optional)
  -p, --picard PATH           Picard JAR file path (default: $PICARD_JAR)
  -a, --adapters FILE         Trimmomatic adapters file (default: $TRIMMOMATIC_ADAPTERS)
  -f, --force                 Force overwrite existing files
  -c, --skip-cleanup          Skip intermediate file cleanup
  -v, --verbose               Verbose output
  -h, --help                  Show this help message

Trimmomatic Parameters:
  --trim-leading INT          Leading quality threshold (default: $TRIM_LEADING)
  --trim-trailing INT         Trailing quality threshold (default: $TRIM_TRAILING)
  --trim-minlen INT           Minimum read length (default: $TRIM_MINLEN)
  --trim-illuminaclip STR     ILLUMINACLIP parameters (default: $TRIM_ILLUMINACLIP)

Example:
  $0 -n sample01 -b input.bam -o output_dir \\
     -r genome.fa -x hisat_index -t 64

  # With conda environment
  $0 -n sample01 -b input.bam -o output_dir \\
     -r genome.fa -x hisat_index -e trimmomatic_env

EOF
    exit 1
}

# =====================================================
# Logging Functions
# =====================================================
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
    exit 1
}

verbose() {
    if [[ "$VERBOSE" == "true" ]]; then
        log "VERBOSE: $*"
    fi
}

# =====================================================
# Parse Command Line Arguments
# =====================================================
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n|--sample-name)
                SAMPLE_NAME="$2"
                shift 2
                ;;
            -b|--input-bam)
                INPUT_BAM="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -r|--reference)
                REF_FASTA="$2"
                shift 2
                ;;
            -x|--hisat-index)
                HISAT_INDEX="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -e|--conda-env)
                CONDA_ENV="$2"
                shift 2
                ;;
            -p|--picard)
                PICARD_JAR="$2"
                shift 2
                ;;
            -a|--adapters)
                TRIMMOMATIC_ADAPTERS="$2"
                shift 2
                ;;
            --trim-leading)
                TRIM_LEADING="$2"
                shift 2
                ;;
            --trim-trailing)
                TRIM_TRAILING="$2"
                shift 2
                ;;
            --trim-minlen)
                TRIM_MINLEN="$2"
                shift 2
                ;;
            --trim-illuminaclip)
                TRIM_ILLUMINACLIP="$2"
                shift 2
                ;;
            -f|--force)
                FORCE_OVERWRITE=true
                shift
                ;;
            -c|--skip-cleanup)
                SKIP_CLEANUP=true
                shift
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -h|--help)
                usage
                ;;
            *)
                error "Unknown option: $1"
                ;;
        esac
    done
}

# =====================================================
# Validation Functions
# =====================================================
validate_inputs() {
    log "Validating inputs..."
    
    # Check required parameters
    [[ -z "$SAMPLE_NAME" ]] && error "Sample name is required (-n/--sample-name)"
    [[ -z "$INPUT_BAM" ]] && error "Input BAM file is required (-b/--input-bam)"
    [[ -z "$OUTPUT_DIR" ]] && error "Output directory is required (-o/--output-dir)"
    [[ -z "$REF_FASTA" ]] && error "Reference FASTA is required (-r/--reference)"
    [[ -z "$HISAT_INDEX" ]] && error "HISAT2 index is required (-x/--hisat-index)"
    
    # Check file existence
    [[ ! -f "$INPUT_BAM" ]] && error "Input BAM file not found: $INPUT_BAM"
    [[ ! -f "$REF_FASTA" ]] && error "Reference FASTA not found: $REF_FASTA"
    
    # Check HISAT2 index files
    if [[ ! -f "${HISAT_INDEX}.1.ht2" ]]; then
        error "HISAT2 index not found: ${HISAT_INDEX}*.ht2"
    fi
    
    # Validate threads
    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
        error "Threads must be a positive integer: $THREADS"
    fi
    
    log "Input validation passed"
}

# =====================================================
# Utility Functions
# =====================================================
setup_conda() {
    if [[ -n "$CONDA_ENV" ]]; then
        verbose "Setting up conda environment: $CONDA_ENV"
        if command -v conda >/dev/null 2>&1; then
            eval "$(conda shell.bash hook)"
        else
            error "Conda not found but conda environment specified"
        fi
    fi
}

activate_conda() {
    if [[ -n "$CONDA_ENV" ]]; then
        verbose "Activating conda environment: $CONDA_ENV"
        conda activate "$CONDA_ENV"
    fi
}

deactivate_conda() {
    if [[ -n "$CONDA_ENV" ]]; then
        verbose "Deactivating conda environment"
        conda deactivate
    fi
}

check_tools() {
    log "Checking required tools..."
    
    local required_tools=("samtools" "hisat2" "bwa" "java")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        error "Missing required tools: ${missing_tools[*]}"
    fi
    
    # Check trimmomatic availability
    if [[ -n "$CONDA_ENV" ]]; then
        activate_conda
        if ! command -v trimmomatic >/dev/null 2>&1; then
            error "Trimmomatic not found in conda environment: $CONDA_ENV"
        fi
        deactivate_conda
    else
        if ! command -v trimmomatic >/dev/null 2>&1; then
            error "Trimmomatic not found. Please specify conda environment with -e option"
        fi
    fi
    
    log "All required tools found"
}

# =====================================================
# File Management Functions
# =====================================================
check_output_exists() {
    local sample_outdir="$OUTPUT_DIR"
    local final_file="$sample_outdir/still_unmap_bwa.s.bam"
    
    if [[ -f "$final_file" && "$FORCE_OVERWRITE" == "false" ]]; then
        log "Output exists for $SAMPLE_NAME, skipping (use -f to force overwrite)"
        return 0
    fi
    return 1
}

cleanup_intermediate_files() {
    local sample_outdir="$1"
    
    if [[ "$SKIP_CLEANUP" == "true" ]]; then
        verbose "Skipping cleanup for $sample_outdir"
        return
    fi
    
    verbose "Cleaning up intermediate files in $sample_outdir"
    
    # Remove intermediate files but keep important ones
    rm -f "$sample_outdir"/*.sam
    rm -f "$sample_outdir"/*.sortn.bam
    rm -f "$sample_outdir"/*unpaired.fastq
    rm -f "$sample_outdir"/*single.fastq
    rm -f "$sample_outdir"/*unpair.fastq.gz
    rm -f "$sample_outdir"/hisat_unmap_*.fastq
    rm -f "$sample_outdir"/hisat.unmap*.stillunmap.fq
    rm -f "$sample_outdir"/hisat.unmap*.stillunmap.*.fq
    rm -f "$sample_outdir"/*.path
}

# =====================================================
# Main Processing Functions
# =====================================================
extract_fastq_from_bam() {
    local sample_outdir="$1"
    
    log "Extracting FASTQ from BAM for $SAMPLE_NAME"
    
    # Sort by name and extract FASTQ
    samtools sort -n -@ "$THREADS" -O BAM -o "$sample_outdir/${SAMPLE_NAME}.sortn.bam" "$INPUT_BAM"
    
    samtools fastq \
        -0 "$sample_outdir/${SAMPLE_NAME}.unpaired.fastq" \
        -1 "$sample_outdir/${SAMPLE_NAME}.1.fastq" \
        -2 "$sample_outdir/${SAMPLE_NAME}.2.fastq" \
        -s "$sample_outdir/${SAMPLE_NAME}.single.fastq" \
        -N "$sample_outdir/${SAMPLE_NAME}.sortn.bam"
}

trim_reads() {
    local sample_outdir="$1"
    
    log "Trimming reads for $SAMPLE_NAME"
    
    # Activate conda environment if specified
    activate_conda
    
    local fq1="$sample_outdir/${SAMPLE_NAME}.1.fastq"
    local fq2="$sample_outdir/${SAMPLE_NAME}.2.fastq"
    local fq1_clean="$sample_outdir/${SAMPLE_NAME}_1.clean.fastq.gz"
    local fq2_clean="$sample_outdir/${SAMPLE_NAME}_2.clean.fastq.gz"
    local fq1_unpair="$sample_outdir/${SAMPLE_NAME}_1.unpair.fastq.gz"
    local fq2_unpair="$sample_outdir/${SAMPLE_NAME}_2.unpair.fastq.gz"
    
    trimmomatic PE -threads "$THREADS" \
        "$fq1" "$fq2" \
        "$fq1_clean" "$fq1_unpair" \
        "$fq2_clean" "$fq2_unpair" \
        "ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:${TRIM_ILLUMINACLIP}" \
        "LEADING:${TRIM_LEADING}" \
        "TRAILING:${TRIM_TRAILING}" \
        "MINLEN:${TRIM_MINLEN}"
    
    deactivate_conda
    
    # Store clean fastq paths for later use
    echo "$fq1_clean" > "$sample_outdir/clean_fq1.path"
    echo "$fq2_clean" > "$sample_outdir/clean_fq2.path"
}

hisat2_alignment() {
    local sample_outdir="$1"
    
    log "Running HISAT2 alignment for $SAMPLE_NAME"
    
    local fq1_clean=$(cat "$sample_outdir/clean_fq1.path")
    local fq2_clean=$(cat "$sample_outdir/clean_fq2.path")
    
    hisat2 -x "$HISAT_INDEX" \
        -1 "$fq1_clean" \
        -2 "$fq2_clean" \
        -S "$sample_outdir/hisat2.sam" \
        -p "$THREADS"
    
    samtools sort -@ "$THREADS" -O BAM -o "$sample_outdir/hisat2.bam" "$sample_outdir/hisat2.sam"
    samtools index -@ "$THREADS" "$sample_outdir/hisat2.bam"
}

process_unmapped_reads() {
    local sample_outdir="$1"
    
    log "Processing unmapped reads for $SAMPLE_NAME"
    
    # Separate mapped and unmapped reads
    local unmap="$sample_outdir/hisat2.unmap.bam"
    local map="$sample_outdir/hisat2.map.bam"
    local unmap_sortn="$sample_outdir/hisat2.unmap.sortn.bam"
    
    samtools view -f4 "$INPUT_BAM" -O BAM -o "$unmap"
    samtools view -F4 "$INPUT_BAM" -O BAM -o "$map"
    samtools sort -n -@ "$THREADS" -O BAM -o "$unmap_sortn" "$unmap"
    
    # Extract unmapped reads to FASTQ
    local unmap_fq_up="$sample_outdir/hisat_unmap_unpaired.fastq"
    local unmap_fq_1="$sample_outdir/hisat_unmap_1.fastq"
    local unmap_fq_2="$sample_outdir/hisat_unmap_2.fastq"
    local unmap_fq_s="$sample_outdir/hisat_unmap_single.fastq"
    
    samtools fastq -0 "$unmap_fq_up" -1 "$unmap_fq_1" -2 "$unmap_fq_2" \
        "$unmap_sortn" -s "$unmap_fq_s" -N -@ "$THREADS"
    
    # Re-align unmapped reads individually
    hisat2 -x "$HISAT_INDEX" -U "$unmap_fq_1" -S "$sample_outdir/hisat.unmap1.remap.sam" -p "$THREADS"
    hisat2 -x "$HISAT_INDEX" -U "$unmap_fq_2" -S "$sample_outdir/hisat.unmap2.remap.sam" -p "$THREADS"
    hisat2 -x "$HISAT_INDEX" -U "$unmap_fq_s" -S "$sample_outdir/hisat.unmapup.remap.sam" -p "$THREADS"
    
    # Merge all re-aligned reads
    samtools merge -f \
        "$sample_outdir/hisat.unmap_remap_all.bam" \
        "$sample_outdir/hisat.unmap1.remap.sam" \
        "$sample_outdir/hisat.unmap2.remap.sam" \
        "$sample_outdir/hisat.unmapup.remap.sam"
    
    samtools sort -@ "$THREADS" -O BAM -o "$sample_outdir/hisat.unmap_remap_all.s.bam" \
        "$sample_outdir/hisat.unmap_remap_all.bam"
    samtools view -F 4 -O BAM -@ "$THREADS" -o "$sample_outdir/hisat.unmap_remap_all.s.map.bam" \
        "$sample_outdir/hisat.unmap_remap_all.s.bam"
    
    # Merge with originally mapped reads
    samtools merge -f "$sample_outdir/hisat.map_unmapremap.bam" \
        "$map" "$sample_outdir/hisat.unmap_remap_all.s.map.bam"
    samtools sort -@ "$THREADS" -O BAM -o "$sample_outdir/hisat.map_unmapremap.s.bam" \
        "$sample_outdir/hisat.map_unmapremap.bam"
}

bwa_rescue_alignment() {
    local sample_outdir="$1"
    
    log "Running BWA rescue alignment for $SAMPLE_NAME"
    
    # Extract still unmapped reads
    samtools view -f4 "$sample_outdir/hisat.unmap1.remap.sam" -O BAM -o "$sample_outdir/hisat.unmap1.stillunmap.bam"
    samtools view -f4 "$sample_outdir/hisat.unmap2.remap.sam" -O BAM -o "$sample_outdir/hisat.unmap2.stillunmap.bam"
    samtools view -f4 "$sample_outdir/hisat.unmapup.remap.sam" -O BAM -o "$sample_outdir/hisat.unmapup.stillunmap.bam"
    
    # Convert to FASTQ using Picard
    java -jar "$PICARD_JAR" SamToFastq \
        -I "$sample_outdir/hisat.unmap1.stillunmap.bam" \
        -F "$sample_outdir/hisat.unmap1.stillunmap.fq"
    
    java -jar "$PICARD_JAR" SamToFastq \
        -I "$sample_outdir/hisat.unmap2.stillunmap.bam" \
        -F "$sample_outdir/hisat.unmap2.stillunmap.fq"
    
    java -jar "$PICARD_JAR" SamToFastq \
        -I "$sample_outdir/hisat.unmapup.stillunmap.bam" \
        -F "$sample_outdir/hisat.unmapup.stillunmap.fq"
    
    # Add read identifiers and combine
    awk '{if(NR%4==1) {print $0 "_r1"} else {print $0}}' \
        "$sample_outdir/hisat.unmap1.stillunmap.fq" > "$sample_outdir/hisat.unmap1.stillunmap.r1.fq"
    awk '{if(NR%4==1) {print $0 "_r2"} else {print $0}}' \
        "$sample_outdir/hisat.unmap2.stillunmap.fq" > "$sample_outdir/hisat.unmap2.stillunmap.r2.fq"
    awk '{if(NR%4==1) {print $0 "_sg"} else {print $0}}' \
        "$sample_outdir/hisat.unmapup.stillunmap.fq" > "$sample_outdir/hisat.unmapup.stillunmap.sg.fq"
    
    cat "$sample_outdir/hisat.unmap1.stillunmap.r1.fq" \
        "$sample_outdir/hisat.unmap2.stillunmap.r2.fq" \
        "$sample_outdir/hisat.unmapup.stillunmap.sg.fq" > \
        "$sample_outdir/still_unmap_r1r2sg.fq"
    
    # BWA alignment
    bwa mem "$REF_FASTA" "$sample_outdir/still_unmap_r1r2sg.fq" -t "$THREADS" > \
        "$sample_outdir/still_unmap_bwa.sam"
    samtools sort "$sample_outdir/still_unmap_bwa.sam" -O BAM -o "$sample_outdir/still_unmap_bwa.s.bam"
    samtools index "$sample_outdir/still_unmap_bwa.s.bam"
}

# =====================================================
# Main Function
# =====================================================
main() {
    local start_time=$(date +%s)
    
    log "Starting RNA-seq Re-alignment Pipeline for sample: $SAMPLE_NAME"
    log "Parameters:"
    log "  Sample name: $SAMPLE_NAME"
    log "  Input BAM: $INPUT_BAM"
    log "  Output dir: $OUTPUT_DIR"
    log "  Reference: $REF_FASTA"
    log "  HISAT2 index: $HISAT_INDEX"
    log "  Threads: $THREADS"
    log "  Conda env: ${CONDA_ENV:-none}"
    log "  Force overwrite: $FORCE_OVERWRITE"
    log "  Skip cleanup: $SKIP_CLEANUP"
    
    # Setup and validation
    validate_inputs
    setup_conda
    check_tools
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Check if output already exists
    if check_output_exists; then
        log "Sample $SAMPLE_NAME already processed, exiting"
        exit 0
    fi
    
    # Process sample
    extract_fastq_from_bam "$OUTPUT_DIR"
    trim_reads "$OUTPUT_DIR"
    hisat2_alignment "$OUTPUT_DIR"
    process_unmapped_reads "$OUTPUT_DIR"
    bwa_rescue_alignment "$OUTPUT_DIR"
    
    # Cleanup intermediate files
    cleanup_intermediate_files "$OUTPUT_DIR"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log "Pipeline completed successfully for sample: $SAMPLE_NAME"
    log "Total time: ${duration}s ($(date -d@${duration} -u +%H:%M:%S))"
    log "Final output: $OUTPUT_DIR/still_unmap_bwa.s.bam"
}

# =====================================================
# Script Entry Point
# =====================================================
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    parse_args "$@"
    main
fi