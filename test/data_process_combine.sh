# Create output directory
mkdir -p candidate_out

sample_name="test"
gene_annotation_bed="/data6/wangxuedong/invas_nc/INVAS/scripts/full_pipe/annotation/hg19_genes.gencode.bed"
# Run candidate detection
bash ../scripts/full_pipe/preprocess/combine_sv_rna.sh   
    $sample_name \ 
    sv_out/$sample_name \
    rna_out/$sample_name/hisat.map_unmapremap.s.bam \
    rna_out/$sample_name/still_unmap_bwa.s.bam \
    candidate_out \
    $gene_annotation_bed \
    delly,manta,lumpy,svaba