# process each chromosome separately
inv_vcf=$1
hg38_genes_bed=$2
still_unmap_bwa_bam=$3
hisat2_bam=$4
local_stringtie=/run/media/wangxuedong/5C40783348454630/cityu/software/RNA/stringtie/stringtie


inv_bed=inv.bed
python invvcf2bed_dup.py $inv_vcf $inv_bed 50000

for chr in {1..22} X; do
    echo Processing $chr
    # example: awk '$1 == "12"' hg38_genes.bed > chr12_genes.bed
    awk -v chr=$chr '$1 == chr' $hg38_genes_bed > chr${chr}_genes.bed
    awk 'BEGIN {OFS="\t"} {$1="chr"$1; print}' chr${chr}_genes.bed > temp.bed && mv temp.bed chr${chr}_genes.bed
    gene_bed=chr${chr}_genes.bed
    # extract inversions on this chromosome

    awk -v chr="chr$chr" '$1 == chr' $inv_bed > "chr${chr}_inv.bed"
    
    # bash find_invs_withingenes2.sh inv.bed /run/media/wangxuedong/5C40783348454630/cityu/software/RNA/simulation/RSEM/test/chr12_genes.bed > inv_gene_map.bed
    python find_invs_withingenes2.py chr${chr}_inv.bed chr${chr}_genes.bed chr${chr}_inv_in_gene.bed
    # python inversion_genes_mapping2.py inv_gene_map.bed inv_gene_map.txt
    python inversion_genes_mapping2.py chr${chr}_inv_in_gene.bed chr${chr}_inv_in_gene.txt
    # python filter_dup_inv.py inv_gene_map.txt inv_gene_map_dedup.txt
    python filter_dup_inv.py chr${chr}_inv_in_gene.txt chr${chr}_inv_in_gene_dedup.txt
    # samtools view still_unmap_bwa.s.bam chr12 -O BAM -o still_unmap_bwa.s.chr12.bam
    samtools view $still_unmap_bwa_bam chr${chr} -O BAM -o still_unmap_bwa.s.chr${chr}.bam
    samtools index still_unmap_bwa.s.chr${chr}.bam
    # samtools view $hisat2_bam chr${chr} -O BAM -o hisat2.chr${chr}.bam
    # samtools index hisat2.chr${chr}.bam
    # python get_high_conf_rg.py still_unmap_bwa.s.chr12.bam chr12 5 --gap 0 > high_conf.2.bed
    # python get_high_conf_rg.py still_unmap_bwa.s.chr12.bam chr12 5 --gap 0 > high_conf.2.bed
    python get_high_conf_rg.py still_unmap_bwa.s.chr${chr}.bam chr${chr} 5 --gap 0 > high_conf.${chr}.bed

    # python filter_conf_inv2.py "/run/media/wangxuedong/One Touch/data/colon/colon/vcf/2653e76e-893a-4777-9bea-c4b92568ca8c/rna_seq/genome/9de9a8b1-b42b-4113-bab5-65f00f2c4b5c/pre_rawdata2/high_conf.2.bed" inv_gene_map_dedup.txt inv_gene_map_dedup_highcof.2.txt
    python filter_conf_inv2.py high_conf.${chr}.bed chr${chr}_inv_in_gene_dedup.txt chr${chr}_inv_in_gene_dedup_highcof.txt

    # get gene region bam for each region in chr${chr}_inv_in_gene_dedup_highcof.txt
    region_file=chr${chr}_inv_in_gene_dedup_highcof.txt

    echo "finished chr${chr}"

    # awk 'BEGIN { FS="\t" } {
    # split($1, a, "[:-]");
    # split($2, b, "[:-]");
    # split($3, c, "[:-]");
    # printf("%s:%d-%d\t%s:%d-%d\t%s:%d-%d\n", a[1], a[2], a[3], a[1], b[2], b[3], a[1], c[2], c[3]);
    # }' "$file_path" | while read -r inversion_region gene_region confident_region; do
    #     # echo "Processing $coord1, $coord2, $coord3 with your_program..."
    #     # 调用你的程序，这里用 echo 命令作为示例
    #     # your_program "$coord1" "$coord2" "$coord3"
    #     # 如果 your_program 是个脚本或程序，确保它在 PATH 中或使用完整路径
    #     echo $inversion_region
    # done



 


    




    
    
done
