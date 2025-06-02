


threads=128
ref=/home/xuedowang2/scratch/1Kpg/download/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

HG00512_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/HG00514/HG00512.bam
HG00513_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/HG00514/HG00513.bam
HG00514_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/HG00514/HG00514.bam
HG00731_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/HG00733/HG00731.bam
HG00732_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/HG00733/HG00732.bam
HG00733_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/HG00733/HG00733.bam
NA19240_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/NA19240/NA19240.bam
NA19238_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/NA19240/NA19238.bam
NA19239_tgs_bam=/scratch/project/cs_shuaicli/pgz/trio/tgs/NA19240/NA19239.bam

cute=~/app/spechla_run/SpecHLA/spechla_env/bin/cuteSV

#### for each sample

for sample in HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA19240 NA19238 NA19239; do
    bam_var=${sample}_tgs_bam
    bam=${!bam_var} # 使用间接引用来获取正确的BAM文件路径
    out_sv_vcf=$sample.cute.sv.vcf
    work_dir=${sample}_cute_work
    if [ ! -d "$work_dir" ]; then
        mkdir "$work_dir"
    fi
    "$cute" "$bam" "$ref" "$out_sv_vcf" "$work_dir" \
        --max_cluster_bias_INS 1000 \
        --diff_ratio_merging_INS 0.9 \
        --max_cluster_bias_DEL 1000 \
        --diff_ratio_merging_DEL 0.5  \
        -q 10 \
        --threads "$threads" \
        --sample "$sample" \
        --genotype \
        -mi 500 -md 500 \
        --report_readid

    # 假设sniffles也在PATH中，或者你需要使用完整路径
    sniffles -i "$bam" -v $sample.sniffles.sv.vcf --threads "$threads"
done
