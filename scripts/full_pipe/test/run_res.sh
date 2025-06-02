sample_name_file=$1




for sample_name in $(cat $sample_name_file); do
    echo $sample_name
    python ../main.py --input_dir "/run/media/wangxuedong/One Touch/virus_rna/ShangHai_GanDan_ICC/res/$sample_name/" --sample_name $sample_name --extend_region 5000 --output_dir ./$sample_name --ref_genome /run/media/wangxuedong/backup1/download/hg19_ref/hg19.fa --rna_map_with_remap_bam "/run/media/wangxuedong/One Touch/virus_rna/ShangHai_GanDan_ICC/RNA/$sample_name/hisat.map_unmapremap.s.bam" --rna_unmap_bwa_bam "/run/media/wangxuedong/One Touch/virus_rna/ShangHai_GanDan_ICC/RNA/$sample_name/still_unmap_bwa.s.bam" --wgs_bam "/run/media/wangxuedong/One Touch/virus_rna/ShangHai_GanDan_ICC/sub_wgs/$sample_name/${sample_name}_focus_sorted.bam"
done

