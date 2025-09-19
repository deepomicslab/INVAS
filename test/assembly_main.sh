#conda activate invas_assembly
sample_name="test"

# Create output directory
mkdir -p assembly_out/$sample_name
reference_genome="/data6/wangxuedong/invas_nc/wgs_ref/hg19.fa"
wgs_bam="/data6/wangxuedong/invas_nc/INVAS/test/demo/DNA/inv_wgs.hg19.bam"
# Run assembly
python ../scripts/full_pipe/main.py \
    --input_dir candidate_out/res/$sample_name \
    --sample_name $sample_name \
    --extend_region 5000 \
    --output_dir assembly_out/$sample_name \
    --ref_genome $reference_genome \
    --rna_map_with_remap_bam rna_out/$sample_name/hisat.map_unmapremap.s.bam \
    --rna_unmap_bwa_bam rna_out/$sample_name/still_unmap_bwa.s.bam \
    --wgs_bam $wgs_bam \
    --with_normal_haps
