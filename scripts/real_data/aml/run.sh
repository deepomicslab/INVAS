filename=$1
output_dir=$2
ref_bed="/run/media/wangxuedong/One Touch/virus_rna/hg38_anno/hg38_genes.gencode.bed"
scripts_dir="/run/media/wangxuedong/5C40783348454630/cityu/software/secflow_scripts/scripts/test_inv/preprocess"
# 检查文件是否存在
if [ ! -f "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi
res_dir="$output_dir"/res/
if [ ! -d "$res_dir" ]; then
    mkdir -p "$res_dir"
fi
# 使用循环读取文件的每一行
while read -r line; do
    echo "precessing $line"
    sample_res_dir="$res_dir"/$line
    if [ ! -d "$sample_res_dir" ]; then
        mkdir -p "$sample_res_dir"
    fi
    hisat2_bam="$output_dir"/RNA/$line/hisat.map_unmapremap.s.bam
    bwa_bam="$output_dir"/RNA/$line/still_unmap_bwa.s.bam
    # for delly
    delly_res_dir="$sample_res_dir"/delly/
    if [ ! -d "$delly_res_dir" ]; then
        mkdir -p "$delly_res_dir"
    fi
    delly_vcf_f="$output_dir"/WGS/$line/delly/$line.delly.f.vcf
    delly_vcf="$output_dir"/WGS/$line/delly/$line.delly.vcf
    python $scripts_dir/filter_delly.py --output_vcf "$delly_vcf_f" --input_vcf "$delly_vcf"
    echo "run delly for $line"
    bash $scripts_dir/preprocess.sh "$delly_vcf_f" "$ref_bed" "$bwa_bam" "$hisat2_bam" "$delly_res_dir"
    # for manta
    manta_vcf="$output_dir"/WGS/$line/manta/$line.manta.invfmt.vcf
    manta_res_dir="$sample_res_dir"/manta/
    if [ ! -d "$manta_res_dir" ]; then
        mkdir -p "$manta_res_dir"
    fi
    echo "run manta for $line"
    bash $scripts_dir/preprocess.sh "$manta_vcf" "$ref_bed" "$bwa_bam" "$hisat2_bam" "$manta_res_dir"
done < "$filename"
